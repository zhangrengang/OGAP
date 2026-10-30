from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import zip
from builtins import object
import sys
import os
import re
import copy
import gzip
from collections import Counter, OrderedDict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# from Bio.SeqUtils import GC
from Bio.SeqUtils import gc_fraction as GC
from Bio import SeqFeature
# from Bio.Alphabet import IUPAC
from lazy_property import LazyWritableProperty as lazyproperty
from .Taxonomy import Taxonomy
from .small_tools import open_file as open

from ete3 import NCBITaxa
ncbi = NCBITaxa()

class GenbankParser(object):
	'''parser of genbank based on SeqIO'''
	def __init__(self, gbfiles):
		if isinstance(gbfiles, str):  # one file
			self.gbfiles = [gbfiles]
		else:
			self.gbfiles = list(gbfiles)
	def __iter__(self):
		return self._parse()

	def _parse(self):
		for gbfile in self.gbfiles:
			if gbfile.endswith(".gz"):
				fp = gzip.open(gbfile, mode="rt")
			else:
				fp = open(gbfile, mode="r")
			try:
				for record in SeqIO.parse(fp, 'genbank'):
					try:
						yield GenbankRecord(record)
					except Exception as e:
						import sys
						print(f"WARNING: skip bad record, error: {e}", file=sys.stderr)
						continue
			finally:
				fp.close()



	def filter_by_taxon(self, taxon=None, exclude=None, taxonomy_dbfile=None):

		i = 0
		if taxonomy_dbfile is None:
			d_species = {}
		else:
			d_species = Taxonomy(jsonfile=taxonomy_dbfile).db
		print('loading genbank records', file=sys.stderr)
		self.tax = []
		for record in self:
			#print >> sys.stderr, record.organism, record.id, record.cds_count, record.rna_count
			if not (record.cds_count + record.rna_count>0):
				continue
			key = record.organism.lower()

			print(f"[DEBUG] rec.org:>{record.organism}< | key:>{key}<", file=sys.stderr)
			dbg_keys = list(d_species.keys())[:10]
			print(f"[DEBUG] d_species keys sample: {dbg_keys}", file=sys.stderr)

			#print >> sys.stderr, record.organism, record.id
			
			if key in d_species:
				taxid, tax_lineage, tax_ranks = d_species[key]
			else:
				print(f'`{record.organism}` ({record.id}) is not found, key={key}', file=sys.stderr)
				continue

			index = record.is_taxon(taxon)
			print(f"[DBG‑is_taxon] taxon={taxon}, index={index}", file=sys.stderr)
			if not index:
				print(f"[DBG‑SKIP] is_taxon False skip {record.id}", file=sys.stderr)
				continue

			if exclude is not None:
				to_exclude = False
				for ex_taxon in exclude:
					if record.is_taxon(ex_taxon):
						to_exclude = True
						break
				if to_exclude:
					print(f"[DBG‑SKIP] excluded skip {record.id}", file=sys.stderr)
					continue

			# 直接写实例字典，绕开 @property 只读限制
			#record.ranks = OrderedDict([(rank, tax) for rank, tax in zip(tax_ranks, tax_lineage)])
			#record.__dict__['taxonomy'] = tax_lineage
			#record.__dict__['ranks'] = tax_ranks
			#record.organism = str(record.ranks['species'])

			# 直接写实例字典，绕开 @property 只读限制
			new_ranks = OrderedDict([(rank, tax) for rank, tax in zip(tax_ranks, tax_lineage)])

			# 1. taxid lineage → 翻译成分类名称字符串列表
			taxid_2_name = ncbi.get_taxid_translator(tax_lineage)
			taxonomy_name_list = [taxid_2_name[tid] for tid in tax_lineage]

			# 2. 获取物种名字
			sp_tid = new_ranks['species']
			sp_name = ncbi.get_taxid_translator([sp_tid])[sp_tid]

			# 全部写入实例 __dict__
			record.__dict__['ranks'] = new_ranks
			record.__dict__['taxonomy'] = taxonomy_name_list
			record.__dict__['organism'] = sp_name
			
			if not record.name_count:	# some taxa without annotation
				continue

			i += 1
			# ✅ 关键修复：self.tax 存入【分类名称字符串切片】，不是taxid数字切片！
			self.tax += [taxonomy_name_list[:index]]
			record.__dict__['taxid'] = taxid
			print(f"[DBG‑KEEP] keep record {record.id}, i={i}", file=sys.stderr)
			yield record


		#print >> sys.stderr, Counter(self.taxa)

	def write_seqs_by_taxon(self, outdir, feat_type='protein', records=None, taxon=None):
		import sys
		features = []
		if records is not None:
			records = list(records)
		else:
			records = list(self.filter_by_taxon(taxon=taxon))
		print(f"[DBG-write_seq] feat_type={feat_type}, total records={len(records)}", file=sys.stderr)
		for rc in records:	# rc is GenbankRecord
			sp_name_safe = rc.species.replace(" ", "_")
			outfile = '{}/{}.fasta'.format(outdir, sp_name_safe)
			print(f"[DBG-write_seq] process species: {rc.species}, seq is None? {rc.seq is None}", file=sys.stderr)
			# =========修复1：protein映射到CDS；rRNA/tRNA元组处理=========
			if feat_type == "protein":
				feat_list = list(rc.get_features("CDS"))
			else:
				feat_list = []
				for ft in feat_type:
					feat_list.extend(rc.get_features(ft))
			print(f"[DBG-write_seq] got {len(feat_list)} features for {rc.species}", file=sys.stderr)
			with open(outfile, 'a') as fout:
				feat_idx = 0
				for feat in feat_list:
					try:
						if feat_type == "protein":
							nuc_seq = feat.location.extract(rc.seq)
							print(f"[DBG‑trans] nuc len={len(nuc_seq)}", file=sys.stderr)
							transl_table = int(feat.qualifiers.get("transl_table", [11])[0])
							prot_seq = nuc_seq.translate(table=transl_table, to_stop=True)
							# 先赋值gene_id，再print、再写fasta
							gene_id = feat.qualifiers.get("gene", feat.qualifiers.get("locus_tag", ["unknown_gene"]))[0]
							gene_id = gene_id.replace(" ", "_")   # fasta标题不能有空格，必加
							unique_gid = f"{gene_id}_{feat_idx}"
							feat_idx += 1
							print(f"[DBG‑trans] prot len={len(prot_seq)} gene_id={unique_gid}", file=sys.stderr)
							fout.write(f">{unique_gid}\n{prot_seq}\n")
						else:
							feat.write(fout, feat_type)
						# if‑else两个分支成功都append
						features.append(feat)
					except Exception as e:
						print(f"[DBG-write_seq] write feat failed: {type(e)} {e}", file=sys.stderr)
			if os.path.exists(outfile) and os.path.getsize(outfile) == 0:
				os.remove(outfile)
		# =========❗❗ return退出，退到for循环外面，全部物种跑完再返回=========
		print(f"[DBG‑write_seq‑RETURN] total collected features = {len(features)}", file=sys.stderr)
		return features

	def get_augustus_train_set(self, genes):
		pass

def format_taxon(taxon):
	return re.compile(r'[^\w\.]').sub('_', taxon)

class GenbankRecord(object):
	'''wrapper of one SeqRecord of genbank'''
	def __init__(self, gb_record):
		self.__dict__ = copy.deepcopy(gb_record.__dict__)
		self.record = gb_record
		self.seq = gb_record.seq
		self.id = gb_record.id
	def prune(self):
		del self.record
		del self.seq
		#del self.features
		return self
	@property
	def species(self):
		#return self.annotations['organism'].replace(' ', '_')
		#return format_taxon(self.annotations['organism'])
		return format_taxon(self.organism)
	@property
	def organism(self):
		return self.annotations['organism']
#		org = self.annotations['organism']
#		if re.compile(r'(subsp|var)').search(org):
#			return re.compile(r'(\S+.?\S+).?(subsp|var)?').match(org).groups()[0]
#		else:
#			return org
	@property
	def genus(self):
		#return self.annotations['organism'].split()[0]
		return self.ranks.get('genus')
	@property
	def family(self):
		return self.ranks.get('family')
	@property
	def order(self):
		return self.ranks.get('order')
#	def get_rank(self, rank):
#		return Taxonomy(spname=self.organism).get_rank(rank)
#	@lazyproperty
#	def ranks(self):	# OrderedDict
#		return Taxonomy(spname=self.organism).get_ranks()
#	@lazyproperty
#	def taxonomy(self):
#		return self.annotations['taxonomy']

	@property
	def ranks(self):	# OrderedDict
		if 'ranks' in self.__dict__:
			return self.__dict__['ranks']
		# 未预先注入，回退原逻辑，依赖 self.Taxonomy
		if hasattr(self, "Taxonomy"):
			return self.Taxonomy.ranks
		else:
			# custom模式没有Taxonomy对象，返回空有序字典
			from collections import OrderedDict
			return OrderedDict()

	@property
	def taxonomy(self):
		if 'taxonomy' in self.__dict__:
			return self.__dict__['taxonomy']
		# 未预先注入，回退原逻辑，依赖 self.Taxonomy
		if hasattr(self, "Taxonomy"):
			return self.Taxonomy.taxonomy
		else:
			# custom模式没有Taxonomy对象，返回空列表
			return []


#	@lazyproperty
#	def Taxonomy(self):
#		return Taxonomy(spname=self.organism)
	@property
	def length(self):
		return len(self.seq)
	@property
	def GC(self):
		return round(GC(self.seq), 1)
	@property
	def cds_count(self):
		return self.count_feature('CDS')
	@property
	def rrn_count(self):
		return self.count_feature('rRNA')
	@property
	def trn_count(self):
		return self.count_feature('tRNA')
	@property
	def rna_count(self):
		return self.rrn_count + self.trn_count
	@property
	def name_count(self):
		return self.count_qualifier('gene')
	
	def count_feature(self, feat_type):
		count = 0
		for feature in self.features:
			if feature.type == feat_type:
				count += 1
		return count
	def count_qualifier(self, qualifier):
		count = 0
		for feature in self.features:
			if feature.type == 'gene' and qualifier in feature.qualifiers \
				and 'pseudo' not in feature.qualifiers and \
				not feature.qualifiers[qualifier][0].startswith('orf'):
				count += 1
		return count

	def is_taxon(self, taxon):
		'''example: 'Viridiplantae' in ['Eukaryota', 'Viridiplantae', 'Streptophyta',]'''
		if taxon is None:
			return True
			try:
				return self.taxonomy.index(taxon) + 1
			except ValueError:
				return False


	def extract_feature_seqs(self, feature_types):
		if isinstance(feature_types, str):
			feature_types = {feature_types}
		elif not isinstance(feature_types, set):
			feature_types = set(feature_types)
		features = set([])
		for i, feature in enumerate(self.features):
			if feature.type in feature_types:
				try:
					feat_id = feature.qualifiers['gene'][0]
				except KeyError:
					try: feat_id = feature.qualifiers['product'][0].replace(' ', '_')
					except KeyError:
						feat_id = 'gene{}'.format(i)
				if feat_id in features:
					feat_id = '{}-{}'.format(feat_id, i)
				try: nucl_seq = feature.extract(self.seq)
				except Exception as e:
					print(e, self.organism, self.id, feat_id, feature.type, feature.location, file=sys.stderr)
					continue
					#raise Exception(e)
				features.add(feat_id)
				yield FeatureRecord(id=feat_id, seq=nucl_seq, index=i, feature=feature)
##						type=feature.type, location=feature.location,
##						qualifiers=feature.qualifiers)
	def write(self, fout):
		SeqIO.write(self.record, fout, 'genbank')

	def write_rrn(self, fout):
		self.write_features(fout, 'rRNA')
	def write_protein(self, fout):
		self.write_features(fout, 'protein')
	def write_cds(self, fout):
		self.write_features(fout, 'CDS')
	def write_features(self, fout, feat_type):
		'''write feature sequences'''
		for feat in self.get_features(feat_type):
			feat.write(fout, feat_type)
	def get_features(self, feat_type):
		if isinstance(feat_type, str):
			feature_types = {feat_type}
		else:
			feature_types = set(feat_type)
		if feat_type == 'protein':
			feature_types.add('CDS')
		for feat in self.extract_feature_seqs(feature_types):
			feat.id = '{}|{}'.format(self.species, feat.id)	# add organism to feat.id
			if feat_type == 'protein':
				feat.seq = feat.pep
			yield feat

class FeatureRecord(object):
	'''wrapper of one SeqFeature'''
	def __init__(self, id, seq, index, feature):
		self.__dict__ = copy.deepcopy(feature.__dict__)
	#	print >>sys.stderr, feature.__dict__
	#	print >>sys.stderr, self.__dict__
		self.id = self.format_id(id)
		self._seq = seq
		self.seq = seq
		self.index = index
	#	self.qualifiers = qualifiers
	#	print >>sys.stderr, vars(self), vars(feature)
		self.description = ''
	def __hash__(self):
		return hash(self.id)
	def __eq__(self, other):
		if self.id == other.id:
			return True
		return False
	def format_id(self, id):
		id = re.compile(r'\W+$').sub('', id)
		id = re.compile(r'\W+').sub('-', id)
		return id
	@lazyproperty
	def gene(self):
		try:
			gene = self.qualifiers['gene'][0]
		except KeyError:
			try:
				gene = self.qualifiers['product'][0]
			except KeyError:
				gene = self.id
		return self.format_id(gene)
	@lazyproperty
	def nexon(self):
		return len(self.location.parts)
	@lazyproperty
	def product(self):
		try:
			product = self.qualifiers['product'][0]
		except KeyError:
			if 'translation' in self.qualifiers:
				product = 'hypothetical protein'
			else:
				product = ''
		return product
	@lazyproperty
	def transl_table(self):
		try:
			table = int(self.qualifiers['transl_table'][0])
		except KeyError:
			table = 1
		return table
	@lazyproperty
	def trans_splicing(self):
		if 'trans_splicing' in self.qualifiers:
			return True
		return self.check_trans_splicing()
	def check_trans_splicing(self):
		strands = [part.strand for part in self.location.parts]
		if len(set(strands)) > 1:
			return True
		starts = [part.start for part in self.location.parts]
		if starts != sorted(starts) and starts != list(reversed(sorted(starts))):
			return True
		return False
	def rearrange_trans_splicing(self, record, flank=200, gap_length=100):
		gap = 'N'*gap_length
		locs, seqs = [], []
		last_end = 0
		for part in self.location.parts:
			start = part.start - flank
			end = part.end + flank
			start = 0 if start < 0 else start
			end = len(record) if end > len(record) else end
			seq = record.seq[start:end]
			if part.strand == -1:
				seq = seq.reverse_complement()
				loc_start = end - part.end
				loc_end  = end - part.start
			else:
				loc_start = part.start - start
				loc_end = part.end - start
			locs += [(loc_start+last_end, loc_end+last_end)]
			seqs += [str(seq)]
			last_end += len(seq) + gap_length
		seq = gap.join(seqs)
		parts = [SeqFeature.FeatureLocation(start, end, strand=1) for start, end in locs]
		if len(parts) > 1:
			location = SeqFeature.CompoundLocation(parts)
			location_operator = 'join'
		else:
			location = parts[0]
			location_operator = None
		feature = SeqFeature.SeqFeature(location, type=self.type, 
						location_operator=location_operator, qualifiers=self.qualifiers)
		source = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(0, len(seq)), type='source', strand=1)
		seq = Seq(seq, IUPAC.ambiguous_dna)
		id = '{}_{}-{}'.format(record.id, self.location.start, self.location.end)
		new_record = SeqRecord(id=id, seq=seq, features=[source, feature])
		return new_record

	@lazyproperty
	def pep(self):
		try:
			return self.qualifiers['translation'][0]
		except KeyError:
			return str(self.seq.translate())
	def write(self, fout, feat_type=None):
		if feat_type == 'protein':
			self.write_pep(fout)
		else:
			self.write_seq(fout)
	def write_seq(self, fout):
		self.write_fasta(fout, id=self.id, seq=self.seq, description=self.description)
	def write_pep(self, fout):
		self.write_fasta(fout, id=self.id, seq=self.pep, description=self.description)
	def write_fasta(self, fout, id, seq, description=''):
		print('>{id} {description}\n{seq}'.format(id=id, description=description, seq=seq), file=fout)
	
