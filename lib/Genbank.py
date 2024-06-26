import sys
import os
import re
import copy
from collections import Counter, OrderedDict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import SeqFeature
from Bio.Alphabet import IUPAC
from lazy_property import LazyWritableProperty as lazyproperty
from Taxonomy import Taxonomy
from small_tools import open_file as open

class GenbankParser():
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
			for record in SeqIO.parse(open(gbfile), 'genbank'):
				yield GenbankRecord(record)
	def filter_by_taxon(self, taxon=None, exclude=None, taxonomy_dbfile=None):
		i = 0
		if taxonomy_dbfile is None:
			d_species = {}
		else:
			d_species = Taxonomy(jsonfile=taxonomy_dbfile).db
		print >>sys.stderr, 'loading genbank records'
		self.tax = []
		for record in self:
			#print >> sys.stderr, record.organism, record.id, record.cds_count, record.rna_count
			if not (record.cds_count + record.rna_count>0):
				continue
			key = record.organism.lower()
			#print >> sys.stderr, record.organism, record.id
			if key in d_species:
				taxid, record.taxonomy, record.ranks = d_species[key]
			else:
				print >>sys.stderr, '`{}` ({}) is not found'.format(record.organism, record.id)
				continue

			index = record.is_taxon(taxon)
			#if record.organism =='Sus scrofa':
			#	print >>sys.stderr, i, record.taxonomy, record.ranks, index
			#print >> sys.stderr, record.organism, record.id
			if not index:
				continue
			if exclude is not None:
				to_exclude = False
				for ex_taxon in exclude:
					if record.is_taxon(taxon):	# exclude
						to_exclude = True
						break
				if to_exclude:
					continue
			record.ranks = OrderedDict([(rank, tax) \
                    for rank, tax in zip(record.ranks, record.taxonomy)])
			record.organism = str(record.ranks['species'])
			if not record.name_count:	# some taxa without annotation
				continue
			#print >> sys.stderr, record.organism, record.id
			i += 1
			self.tax += [record.taxonomy[:index]]
			yield record
		self.ntaxa = i
		#print >> sys.stderr, Counter(self.taxa)
	def write_seqs_by_taxon(self, outdir, feat_type='protein', records=None, taxon=None):
		features = []
		if records is None:
			records = self.filter_by_taxon(taxon=taxon)
		for rc in records:	# rc is GenbankRecord
			outfile = '{}/{}.fasta'.format(outdir, rc.species)
			with open(outfile, 'a') as fout:
				for feat in rc.get_features(feat_type):
					feat.write(fout, feat_type)
					features.append(feat)
			if os.path.getsize(outfile) == 0:
				os.remove(outfile)
		return features
	def get_augustus_train_set(self, genes):
		pass

def format_taxon(taxon):
	return re.compile(r'[^\w\.]').sub('_', taxon)

class GenbankRecord():
	'''wrapper of one SeqRecord of genbank'''
	def __init__(self, gb_record):
		self.__dict__ = copy.deepcopy(gb_record.__dict__)
		self.record = gb_record
		self.seq = gb_record.seq
		self.id = gb_record.id
	def prune(self):
		del self.record
		del self.seq
		del self.features
		return self
	@lazyproperty
	def species(self):
		#return self.annotations['organism'].replace(' ', '_')
		#return format_taxon(self.annotations['organism'])
		return format_taxon(self.organism)
	@lazyproperty
	def organism(self):
		return self.annotations['organism']
#		org = self.annotations['organism']
#		if re.compile(r'(subsp|var)').search(org):
#			return re.compile(r'(\S+.?\S+).?(subsp|var)?').match(org).groups()[0]
#		else:
#			return org
	@lazyproperty
	def genus(self):
		#return self.annotations['organism'].split()[0]
		return self.ranks.get('genus')
	@lazyproperty
	def family(self):
		return self.ranks.get('family')
	@lazyproperty
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
	@lazyproperty
	def ranks(self):	# OrderedDict
		return self.Taxonomy.ranks
	@lazyproperty
	def taxonomy(self):
		return self.Taxonomy.taxonomy
#	@lazyproperty
#	def Taxonomy(self):
#		return Taxonomy(spname=self.organism)
	@lazyproperty
	def length(self):
		return len(self.seq)
	@lazyproperty
	def GC(self):
		return round(GC(self.seq), 1)
	@lazyproperty
	def cds_count(self):
		return self.count_feature('CDS')
	@lazyproperty
	def rrn_count(self):
		return self.count_feature('rRNA')
	@lazyproperty
	def trn_count(self):
		return self.count_feature('tRNA')
	@lazyproperty
	def rna_count(self):
		return self.rrn_count + self.trn_count
	@lazyproperty
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
					feat_id = 'gene{}'.format(i)
				if feat_id in features:
					feat_id = '{}-{}'.format(feat_id, i)
				try: nucl_seq = feature.extract(self.seq)
				except Exception as e:
					print >> sys.stderr, e, self.organism, self.id, feat_id, feature.type, feature.location
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

class FeatureRecord():
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
		print >> fout, '>{id} {description}\n{seq}'.format(id=id, description=description, seq=seq)
	
