from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import map
from builtins import object
from past.utils import old_div
import sys, os
import csv
import gzip
import json
from ete3 import NCBITaxa
import re
import time
import shutil
import glob
import argparse
import numpy as np
from collections import Counter, OrderedDict
from Bio import SeqIO
from Bio import Phylo
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from lazy_property import LazyWritableProperty as lazyproperty
import subprocess
import logging
logger = logging.getLogger(__name__)
import tempfile
from io import StringIO
import tarfile

from .Genbank import GenbankParser, format_taxon
from .OrthoFinder import OrthoFinder, OrthoMCLGroupRecord
from .tRNA import tRNAscanRecord
from .RunCmdsMP import logger, run_cmd, run_job
from .small_tools import mkdirs, rmdirs

AUGUSTUS_CONFIG_PATH=os.environ['AUGUSTUS_CONFIG_PATH']

d_lineage = {}
class Database(object):
	def __init__(self, organ=None, taxon=None, dbrootdir=None,
				gbfiles=None, custom=None, version=None,
				include_orf=False,
				tmpdir='/dev/shm/tmp', 
				ncpu=24,
				upper_limit=500,
				lower_limit=2,
				exclude_taxa=None,
				allow_single_taxon=False,   # 新增，默认False，不允许单taxon
				cleanup=True,
				**kargs):
		if dbrootdir is None:
			dbrootdir = self.get_dbdir()
		self.taxon = taxon
		self.custom = custom

		if custom is None:
			if taxon is None:
				taxon = 'custom'
				#ogtype = '{}/custom'.format(organ, )	# db name and prefix
			else:
				taxon = format_taxon(taxon)
				#ogtype = '{}/{}'.format(organ, format_taxon(taxon))  # e.g. 'mt-Lamiales'
		else:
			taxon = format_taxon(custom)
			#ogtype = '{}/{}'.format(organ, format_taxon(custom))
		ogtype = '{}-{}'.format(organ, taxon)
		self.ogtype = ogtype
		self.organ = organ
		self.dbrootdir = dbrootdir
		self.dborgndir = '{}/{}'.format(dbrootdir, organ)
		#self.taxon = taxon
		self.ncpu = ncpu
		self.gbfiles = gbfiles
		self.version = version
		self.include_orf = include_orf
		self.exclude_taxa = exclude_taxa
		self.cleanup = cleanup
		if self.exclude_taxa is not None:
			ogtype += '_ex_{}'.format(
					format_taxon( 
						'_'.join(self.exclude_taxa)))

		self.upper_limit = upper_limit
		self.lower_limit = lower_limit
		# folders
		#dbdir = '{}/{}'.format(dbrootdir, ogtype)
		dbdir = '{}/{}'.format(self.dborgndir, taxon)
		tmpdir = '{}/{}'.format(tmpdir, ogtype)
		self.dbdir = dbdir = os.path.realpath(dbdir)
		self.tmpdir = tmpdir = os.path.realpath(tmpdir)
		self.pepdir = '{}/prot'.format(tmpdir)	# tmp dir of CDS sequence for OrthoFinder
		self.rnadir = '{}/rrna'.format(tmpdir)	# tmp dir of RNA sequence for OrthoFinder
		self.alndir = '{}/alns'.format(tmpdir)	# tmp dir for alignment
		self.trndir = '{}/trng'.format(tmpdir)	# tmp dir for training augustus
		self.seqdir = '{}/seqs'.format(dbdir)	# sequence for exonerate
		self.msadir = '{}/msas'.format(dbdir)	# multiple sequence alignment
		self.hmmdir = '{}/hmms'.format(dbdir)	# hmm-profile for HMMER
		self.pfldir = '{}/prfl'.format(dbdir)	# protein-profile for augustus
		self.tnddir = '{}/agst'.format(dbdir)	# species for augustus

		# files
		self.version_file = '{}/version'.format(self.dbdir)
		self.cds_file = '{}/cds.info'.format(self.dbdir)
		self.rna_file = '{}/rna.info'.format(self.dbdir)
		self.species_file = '{}/genome.info'.format(self.dbdir)
		self.taxonomy_file = '{}/taxonomy.info'.format(self.dbdir)
		self.name_mapping = '{}/{}.name_mapping'.format(self.dbrootdir, self.organ, )
		self.name_banning = '{}/{}.name_banning'.format(self.dbrootdir, self.organ, )
		self.taxonomy_dbfile = '{}/taxonomy.json.gz'.format(self.dbrootdir)
	def check_augustus(self):
		if not os.access(AUGUSTUS_CONFIG_PATH, os.W_OK):
			raise ValueError('{} is not writable'.format(AUGUSTUS_CONFIG_PATH))
	def listdb(self):
		for db in self.getdb():
			print(db.organ, db.taxon, file=sys.stdout)
	def getdb(self):
		for name in sorted(os.listdir(self.dborgndir)):
			filename = os.path.join(self.dborgndir, name)
			if os.path.isdir(filename):
				try:
					taxon = name.split('-')[0]
				except ValueError:
					continue
				yield Database(organ=self.organ, taxon=taxon)
	def select_db(self, organism):
		lineage = self.get_taxonomy(organism)
		dbs = []
		for db in self.getdb():
			#if not db.organ == self.organ:
			#	continue
			try: db.index = lineage.index(db.taxon)
			except ValueError: continue
			dbs += [db]
		return max(dbs, key=lambda x: x.index)
	
	def get_taxonomy(self, organism):
		organism = organism.split()[0]	# use genus, which is more common
		try:
			lineage = d_lineage[organism] 
			print('Using existed', organism, file=sys.stderr)
			return lineage
		except KeyError: pass
		cmd = 'ete3 ncbiquery --search "{}" --info'.format(organism)
		stdout, stderr, status = run_cmd(cmd, log=True)
		for info in Ete3TaxonomyInfo(stdout):
			lineage = info.named_lineage
			d_lineage[organism] = lineage
			return lineage
			
	def checkdb(self, untar=True):
		# dbdir
		logger.info('using database `{}` in `{}`'.format(self.ogtype, self.dbdir))
		# version
		self.version = open(self.version_file).read().strip()
		logger.info('database version: `{}`'.format(self.version))
		# taxa
		self.species_info = SpeciesInfo(self.species_file)
		self.ntaxa = self.species_info.ntaxa
		logger.info('database based on `{}` taxa'.format(self.ntaxa))
		# taxonomy
		self.taxonomy_info = TaxonomyInfo(self.taxonomy_file)
		self.transl_table = self.taxonomy_info.transl_table
		logger.info('using transl_table `{}`'.format(self.transl_table))
		# gene
		self.cds_info = GeneInfo(self.cds_file, include_orf=self.include_orf)
		self.rna_info = GeneInfo(self.rna_file)
		self.cds_genes = self.cds_info.genes
		self.rna_genes = self.rna_info.genes
		self.name_info = NameInfo(self.name_mapping)
		name_dict = self.name_info.dict
		self.ban_info  = NameInfo(self.name_banning)
		ban_dict = self.ban_info.dict
		self.gene_info = OrderedDict()
		ban_names = set([])
		# update unified gene name and product
		for gene in self.cds_genes + self.rna_genes:
			name, product = gene.name, gene.product
			key = NameInfo.convert_name(name)
			self.gene_info[gene.id] = gene
			if set([name, key]) & set(ban_dict):
				ban_names.add(name)
				logger.info('gene {} is banned from file {}'.format(name, self.name_banning))
				continue
			if not name.startswith('orf') and key not in name_dict:
				logger.warn('{} is not found in name mapping'.format(name))
				continue
			elif name.startswith('orf'):
				continue
			mapped_gene = name_dict[key]
			mapped_name, mapped_product = mapped_gene.name, mapped_gene.product
			#print >> sys.stderr, [name, key, mapped_name], set(ban_dict), set([name, key, mapped_name]) & set(ban_dict)
			if set([name, key, mapped_name]) & set(ban_dict):
				ban_names.add(name)
				logger.info('gene {} is banned from file {}'.format(name, self.name_banning))
				continue
			if name != mapped_name:
				logger.info('mapping gene name {} -> {}'.format(name, mapped_name))
			if product != mapped_product:
				logger.info('mapping gene product {} -> {}'.format(product, mapped_product))
			gene.name, gene.product = mapped_name, mapped_product
#			self.gene_info[gene.id] = gene
		#logger.info('{} CDS, {} RNA'.format(len(self.cds_genes), len(self.rna_genes)))
		# ban
		self.cds_genes = [gene for gene in self.cds_genes if gene.name not in ban_names]
		self.rna_genes = [gene for gene in self.rna_genes if gene.name not in ban_names]
		#print >> sys.stderr, sorted([gene.name for gene in self.cds_genes])
		logger.info('{} CDS, {} RNA'.format(len(self.cds_genes), len(self.rna_genes)))
		self.cds_hmmfiles = [self.get_hmmfile(gene) for gene in self.cds_genes]
		self.rna_hmmfiles = [self.get_hmmfile(gene) for gene in self.rna_genes]
		self.cds_pflfiles = [self.get_pflfile(gene) for gene in self.cds_genes]
		self.augustus_species = [self.get_augustus_species(gene) for gene in self.cds_genes]
		# untar
		if untar:
			self.check_augustus()
			self.untgz_dirs(self.seqdir, self.hmmdir, self.pfldir, self.tnddir)
			self.check_exists(*(self.cds_hmmfiles+self.rna_hmmfiles+self.cds_pflfiles))

	def check_exists(self, *files):
		for file in files:
			if not os.path.exists(file):
				logger.warn('{} do not exists'.format(file))
	def prune_by_rank(self, records, rank, topn):
		d_records = {}  # by rank
		for record in records:
			taxon = getattr(record, rank)
			if taxon is None:
				continue
			try: d_records[taxon] += [record]
			except KeyError: d_records[taxon] = [record]
		reserved_records = []
		for taxon, bin_records in list(d_records.items()):
			bin_records = sorted(bin_records, key=lambda x:-x.name_count)
			reserved_records += bin_records[:topn]
		return reserved_records
		
	def prune(self, records):
		organisms = {record.organism for record in records}
		norganisms = len(organisms)
		#logger.info('total {} organisms'.format(norganisms))
		if norganisms <= self.upper_limit:
			return records
		elif norganisms <= self.lower_limit:
			logger.error('too less records to build database')
		# prune by genus evenly
		rank = 'genus'
		genera = {record.genus for record in records}
		ngenera = len(genera)
		if ngenera <= self.upper_limit:
			n_per_genus = int(round(old_div(1.0*self.upper_limit, ngenera), 0))
			logger.info('pruning by {} with > {} organisms'.format(rank, n_per_genus))
			return self.prune_by_rank(records, rank, n_per_genus)
		else:
			records = self.prune_by_rank(records, rank, 1)
		
		# prune by family evenly
		rank = 'family'
		families = {record.family for record in records if record.family is not None}
		nfamily = len(families)
		if nfamily <= self.upper_limit:
			n_per_family = int(round(old_div(1.0*self.upper_limit, nfamily), 0))
			logger.info('pruning by {} with > {} organisms'.format(rank, n_per_family))
			d_records = {}  # by family
			return self.prune_by_rank(records, rank, n_per_family)
		else:
			records = self.prune_by_rank(records, rank, 1)
		
		# prune by order evenly, regardless of upper limit
		rank = 'order'
		orders = {record.order for record in records if record.order is not None}
		norder = len(orders)
		n_per_order = int(round(old_div(1.0*self.upper_limit, norder), 0))
		if n_per_order < 1:
			n_per_order = 1
		logger.info('pruning by {} with > {} organisms'.format(rank, n_per_order))
		return self.prune_by_rank(records, rank, n_per_order)
	
	def get_species_tree(self, taxa, name_to_taxid=None):
		from ete3 import NCBITaxa
		ncbi = NCBITaxa()
		taxa = set(taxa)
		# 注意：taxa 现在是物种名字符串集合，不是taxid，单物种分支暂时先禁用，测试库要求>=2物种
		if len(taxa) == 1:
			raise RuntimeError(
				f"Only {len(taxa)} valid taxon left. "
				"Building database requires ≥2 valid taxa. "
				"Set allow_single_taxon=True only for test purpose."
			)

		species_tree = '{}/species.treefile'.format(self.dbdir)

		# 完全绕过 ete3 ncbitree 命令行路径，直接使用API构建树
		if name_to_taxid is not None:
			taxid_list = []
			for sp_name in taxa:
				key = sp_name.lower()
				tid = name_to_taxid.get(key)
				if tid is None:
					raise KeyError(f"Cannot find taxid for species: {sp_name}")
				taxid_list.append(tid)
			tree = ncbi.get_topology(taxid_list, intermediate_nodes=True)
			# 将叶子taxid数字替换为原始物种拉丁名
			tid2name = {v: k for k, v in name_to_taxid.items()}
			for leaf in tree.get_leaves():
				try:
					tid = int(leaf.name)
				except ValueError:
					continue
				if tid in tid2name:
					leaf.name = tid2name[tid]
		else:
			tree = ncbi.get_topology(list(taxa), intermediate_nodes=True)

		tree.write(outfile=species_tree, format=0)
		return species_tree
		
	def count_taxa(self):
		logger.info(f"DEBUG count_taxa: self.taxonomy_dbfile = {repr(self.taxonomy_dbfile)}")
		# 第一轮：仅统计，不prune，新建独立parser
		gb_stat = GenbankParser(self.gbfiles)
		# custom模式 或者 taxonomy_dbfile=None，直接全部读取，不走filter_by_taxon
		if self.taxonomy_dbfile is None or self.custom is not None:
			records_stat = [rec for rec in gb_stat]
		else:
			records_stat = [rec for rec in gb_stat.filter_by_taxon(taxon=self.taxon)]
		
		logger.info(f"DEBUG: self.custom = {self.custom}")
		if self.custom is None:
			self.count_rank(records_stat)
		else:
			logger.info("custom mode: skip count_rank, taxonomy stats not required; lineage info saved in taxonomy.json.gz for tree/augustus")

		# read into MEM：重新新建parser，拿到完整带seq的记录
		gb_data = GenbankParser(self.gbfiles)
		# custom模式强制跳过filter_by_taxon，规避species rank KeyError
		if self.taxonomy_dbfile is None or self.custom is not None:
			logger.warning("taxonomy_dbfile is None OR custom mode enabled, skip filter_by_taxon, load all gb records")
			records = list(gb_data)
		else:
			records = list(gb_data.filter_by_taxon(taxon=self.taxon,
								exclude=self.exclude_taxa,
								taxonomy_dbfile=self.taxonomy_dbfile))
		# ==========新增：初始化并填充 self.cds_dict==========
		self.cds_dict = dict()
		self.prot_id2gene = dict()
		from Bio.SeqRecord import SeqRecord
		for rec in records:
			for feat in rec.features:
				if feat.type == "CDS":
					if "locus_tag" in feat.qualifiers:
						gene_id = feat.qualifiers["locus_tag"][0]
					elif "gene" in feat.qualifiers:
						gene_id = feat.qualifiers["gene"][0]
					else:
						continue
					try:
						cds_seq = feat.location.extract(rec.seq)
					except Exception:
						continue
					cds_rec = SeqRecord(cds_seq, id=gene_id, annotations=feat.qualifiers)
					self.cds_dict[gene_id] = cds_rec
					# =========新增：CDS里面直接拿protein_id，补充prot_id2gene =========
					pid = feat.qualifiers.get("protein_id", [None])[0]
					gname = feat.qualifiers.get("gene", [None])[0]
					if pid and gname:
						self.prot_id2gene[pid] = gname
			#	if feat.type == "protein":
			#		pid = feat.qualifiers.get("protein_id", [None])[0]
			#		gname = feat.qualifiers.get("gene", [None])[0]
			#		if pid and gname:
			#			self.prot_id2gene[pid] = gname
		logger.info(f"count_taxa: collected {len(self.cds_dict)} CDS records into self.cds_dict")
		logger.info(f"count_taxa: collected {len(self.prot_id2gene)} protein_id -> gene mapping")
		return records
	
	def count_rank(self, records):
		organisms = {record.organism for record in records}
		norganisms = len(organisms)
		genera = {record.genus for record in records}
		ngenera = len(genera)
		families = {record.family for record in records if record.family is not None}
		nfamily = len(families)
		orders = {record.order for record in records if record.order is not None}
		norder = len(orders)
		info = '{} records: {} orders, {} families, {} genera, {} organism'.format(
				len(records), norder, nfamily, ngenera, norganisms)
		logger.info(info)
		return len(records)
		
	def tgz_dirs(self, *dirs):
		ckp = '{basename}.x'
		cmd = 'cd {dirname} && tar czf {basename}.tgz {basename} --remove-files'
		self.tar_dirs(cmd, ckp, *dirs)
	def untgz_dirs(self, *dirs):
		ckp = '{basename}'
		cmd = 'cd {dirname} && tar xzf {basename}.tgz'
		self.tar_dirs(cmd, ckp, *dirs)
	def tar_dirs(self, cmd, ckp, *dirs):
		for _dir in dirs:
			_dir = _dir.rstrip('/')
			dirname, basename = os.path.dirname(_dir), os.path.basename(_dir)
			_ckp = ckp.format(basename=_dir)
			print('check point: ' + _ckp, file=sys.stderr)
			if not os.path.exists(_ckp):
				_cmd = cmd.format(dirname=dirname, basename=basename)
				run_cmd(_cmd, log=True)

	def export_ref_sequences(self, records):
		"""
		导出数据库参考序列与gene信息
		输出：
			seqdir/: protein.fasta, cds.fasta
			数据库根目录: cds.info, protein.info, rna.info
		"""
		prot_out = os.path.join(self.seqdir, "protein.fasta")
		cds_out = os.path.join(self.seqdir, "cds.fasta")
		# info文件输出数据库根目录 self.dir
		cds_info_out = os.path.join(self.dbdir, "cds.info")
		prot_info_out = os.path.join(self.dbdir, "protein.info")
		rna_info_out = os.path.join(self.dbdir, "rna.info")
		prot_fh = open(prot_out, "w")
		cds_fh = open(cds_out, "w")
		# =========注释下面这一行=========
		# cdsinfo_fh = open(cds_info_out, "w")
		protinfo_fh = open(prot_info_out, "w")
		# =====改动点1：注释掉rnainfo_fh打开句柄=====
		# rnainfo_fh = open(rna_info_out, "w")
		for rc in records:
			spec = rc.species.replace(" ","_")
			# CDS部分
			for feat in rc.get_features("CDS"):
				try:
					locus = feat.qualifiers.get("gene", ["unknown_gene"])[0]
					protein_id = feat.qualifiers.get("protein_id", ["no_protid"])[0]
					cds_seq = feat.location.extract(rc.seq)
					cds_id = f">{spec}|{locus}|{protein_id}"
					cds_fh.write(f"{cds_id}\n{cds_seq}\n")
					# =========注释掉cds.info原始写入行=========
					# cdsinfo_fh.write(f"{locus}\t{spec}\t{protein_id}\t{len(cds_seq)}\n")
					tt = int(feat.qualifiers.get("transl_table", [11])[0])
					prot_seq = cds_seq.translate(table=tt, to_stop=True)
					prot_id = f">{spec}|{locus}|{protein_id}"
					prot_fh.write(f"{prot_id}\n{prot_seq}\n")
					protinfo_fh.write(f"{locus}\t{spec}\t{protein_id}\t{len(prot_seq)}\n")
				except Exception as e:
					print(f"[DBG‑export_ref] skip CDS feature {e}", file=sys.stderr)
					continue
			# =====改动点2：注释整个RNA写入循环=====
			# # rRNA + tRNA，写入 rna.info
			# for feat in list(rc.get_features("rRNA")) + list(rc.get_features("tRNA")):
			# 	try:
			# 		gene_name = feat.qualifiers.get("gene", ["unknown_rna"])[0]
			# 		product = feat.qualifiers.get("product", [""])[0]
			# 		rna_seq = feat.location.extract(rc.seq)
			# 		rnainfo_fh.write(f"{gene_name}\t{spec}\t{product}\t{len(rna_seq)}\n")
			# 	except Exception as e:
			# 		print(f"[DBG‑export_ref] skip RNA feature {e}", file=sys.stderr)
			# 		continue
		prot_fh.close()
		cds_fh.close()
		# =========注释这一行=========
		# cdsinfo_fh.close()
		protinfo_fh.close()
		# =====改动点3：注释rnainfo_fh.close()=====
		# rnainfo_fh.close()
		logger.info(f"export ref seqs done: {prot_out}, {cds_out}")
		# =====改动点4：保留日志，但rna.info还在，只是不再由这个函数生成=====
		logger.info(f"export info files: {cds_info_out}, {prot_info_out}, {rna_info_out}")

	# def makedb(self):
	def makedb(self, ncpu=None, orthofinder_exe=None):
		rmdirs(self.pepdir, self.rnadir)  # for OrthoFinder
		mkdirs(self.pepdir, self.rnadir, self.alndir, self.msadir, self.trndir) 
		mkdirs(self.seqdir, self.hmmdir, self.pfldir, self.tnddir)

		self.cds_dict = dict()

		# ==========【在这里整块插入】==========
		# 优先使用传入ncpu，更新实例属性
		if ncpu is not None:
			self.ncpu = ncpu
		# 解析orthofinder路径
		if orthofinder_exe is None:
			orthofinder_exe = shutil.which("orthofinder")
		if orthofinder_exe is None or not os.path.isfile(orthofinder_exe):
			raise RuntimeError(
				"OrthoFinder‑3.x not found!\n"
				"OrthoFinder is NOT bundled with OGAP.\n"
				"Please install OrthoFinder‑3.x, either put OrthoFinder.py into PATH,\n"
				"or specify via --orthofinder-exe /path/to/OrthoFinder.py"
			)
		logger.info(f"Using OrthoFinder script: {orthofinder_exe}")
		# ======================================
		logger.info('parsing {}'.format(self.gbfiles))
		# ===== custom模式：第一轮扫描gb，构建临时taxonomy映射json =====
		if self.custom is not None:
			ncbi = NCBITaxa()
			temp_species_map = {}
			gb_scan = GenbankParser(self.gbfiles)
			# 物种名称清洗函数：去掉 subsp./var./f./cv./Group 亚种栽培后缀
			def clean_species_name(s):
				s = re.sub(r"\s+(subsp\.|var\.|f\.|cv\.).*", "", s)
				s = re.sub(r"\s+[A-Za-z ]+Group.*", "", s)
				return s.strip()
			for rec in gb_scan:
				org = rec.organism
				key = org.lower()
				if key in temp_species_map:
					continue
				org_clean = clean_species_name(org)
				# 使用清洗后的简化学名查询NCBI
				name_map = ncbi.get_name_translator([org_clean])
				if org_clean not in name_map or len(name_map[org_clean]) == 0:
					logger.warning(f"cannot find taxid for organism: {org} (clean name:{org_clean}), skip")
					continue
				taxid = name_map[org_clean][0]
				lineage = ncbi.get_lineage(taxid)
				rank_dict = ncbi.get_rank(lineage)
				tax_ranks = [rank_dict[tid] for tid in lineage]
				# key仍然用原始完整物种名小写
				temp_species_map[key] = [taxid, lineage, tax_ranks]
			
			# =========新增判断：空map不生成损坏json.gz=========
			if len(temp_species_map) == 0:
				logger.warning("custom mode: NO species taxid found from NCBI, skip taxonomy filter")
				taxonomy_dbfile = None
			else:
				temp_tax_json = os.path.join(self.dbdir, "taxonomy.json.gz")
				with gzip.open(temp_tax_json, "wt") as fp_out:
					json.dump(temp_species_map, fp_out)
				taxonomy_dbfile = temp_tax_json
			# =========【极其关键！！把局部变量写回实例属性】=========
			self.taxonomy_dbfile = taxonomy_dbfile
		else:
			taxonomy_dbfile = self.taxonomy_dbfile
		# ✅ 现在taxonomy_dbfile已经就绪，再调用count_taxa
		records = self.count_taxa()
		#gb = GenbankParser(self.gbfiles)
		if self.version is None:
			version = time.strftime("%Y-%m-%d", time.gmtime(os.path.getmtime(self.gbfiles[0])))
			today = time.strftime("%Y-%m-%d", time.gmtime(time.time()))
			if re.compile(r'\S+.genomic.gbff').search(self.gbfiles[0]):
				source = 'RefSeq'
			else:
				source = 'Unknown'
			self.version = '{} {}: build {}'.format(source, version, today)
			logger.info('version: {}'.format(self.version))
			with open(self.version_file, 'w') as fout:
				print(self.version, file=fout)


		# 新增：主动触发CDS翻译，生成内存protein虚拟feature
		#for rec in records:
    		# 遍历CDS，触发内部翻译逻辑，把protein feature挂载到rec对象
			#_ = list(rec.get_features("CDS"))


		# count
		#self.count_rank(records)
		# count
		if self.custom is None:
			self.count_rank(records)
		else:
			logger.info("custom mode: skip count_rank in makedb, taxonomy stats are not required")

		# prune
		#logger.info('pruning taxa')
		#records = self.prune(records)
		#self.count_rank(records)

		#self.count_rank(records)

		# prune

		#logger.info('pruning taxa')

		#records = self.prune(records)

		#self.count_rank(records)

		# --- 新增：prune之后构建 name -> taxid 映射 ---
		name_to_taxid = {}
		if self.custom is None:
			# 原生模式：从record对象读取taxid
			for rec in records:
				sp_name = rec.organism.strip()
				tid = rec.__dict__['taxid']
				name_to_taxid[sp_name.lower()] = tid
		else:
			# custom模式：从提前生成的taxonomy.json.gz加载taxid，不从record读取
			logger.info("custom mode: load name_to_taxid from taxonomy.json.gz")
			with gzip.open(self.taxonomy_dbfile, "rt") as fp:
				temp_species_map = json.load(fp)
			# temp_species_map key:物种名小写，value=[taxid, lineage, tax_ranks]
			for sp_lower, val in temp_species_map.items():
				tid = val[0]
				name_to_taxid[sp_lower] = tid


		taxa = [record.organism for record in records]

		logger.info('building taxonomy tree')
		# 传入映射字典给 get_species_tree
		sp_tree = self.get_species_tree(set(taxa), name_to_taxid=name_to_taxid)


		#taxa = [record.organism for record in records]
		#logger.info('building taxonomy tree')
		#sp_tree = self.get_species_tree(set(taxa))
		#logger.info('using speceis tree `{}`'.format(sp_tree))
		logger.info('using species tree `{}`'.format(sp_tree))

		prfl_entries = []   # 全局：存储蛋白OGxxx + RNA基因名，统一输出model.lst

		for feat_type, seqdir, prefix, of_opts, build_pfl, min_ratio, cdhit, out_info in zip(
				['protein', ('rRNA','tRNA')],
				[self.pepdir, self.rnadir],
				[self.cds_file, self.rna_file],
				['', '-S blastn'],
				[True, False],
				[0.08, 0.02],
				['cd-hit', 'cd-hit-est'],
				[True, False],
				):
			self.cdhit = cdhit
			logger.info('writing {} sequences into {}'.format(feat_type, seqdir))
			mkdirs(seqdir)

			# 修复Bug1：使用 self.write_seqs_by_taxon，不是 gb.write_seqs_by_taxon
			gb = GenbankParser(self.gbfiles)
			features = gb.write_seqs_by_taxon(seqdir, feat_type=feat_type, records=records, taxon=self.taxon)

			if not features:
				logger.warning(f"No features for {feat_type}, skip this feature group")
				continue

			# ============ OF 调用块 START ============
			if feat_type == "protein":
				#from ogap_build_db import OF_SCRIPT, get_of_subenv
				from ogap_build_db import get_of_subenv

				# =========【这里插入AGST开关判断代码】==========
				# AGST执行开关：判断是否存在gb/gbff注释文件
				run_agst_augustus = False
				gb_exts = {".gb", ".gbff", ".gb.gz", ".gbff.gz"}
				if hasattr(self, "gbfiles") and self.gbfiles:
					for gb_file in self.gbfiles:
						low_name = os.path.basename(gb_file).lower()
						if any(low_name.endswith(ext) for ext in gb_exts):
							run_agst_augustus = True
							logger.info("Detect GB/GBFF annotation file, enable AGST augustus training")
							break
				else:
					logger.info("No GB/GBFF input found, skip AGST augustus training")
				logger.debug(f"run_agst_augustus = {run_agst_augustus}")
				# ==============================================

				# =====================【✅步骤1代码，每行前置4个Tab】=====================
				# ============ 新增：预解析全部GB，构建CDS/基因索引，用于后续AGST拼装gb ============
				#from Bio import SeqIO
				gb_gene_index = {}  # key: (organism, gene_name), value: list of (SeqRecord, cds_feature)
				gb_cds_id_map = {}  # key: protein_id/locus_tag, value: (organism, gene_name, SeqRecord, cds_feature)

				if run_agst_augustus and hasattr(self, "gbfiles") and self.gbfiles:
					logger.info("Pre-parse GenBank files for AGST GB reconstruction ...")
					for gb_path in self.gbfiles:
						try:
							for rec in SeqIO.parse(gb_path, "genbank"):
								org_name = rec.annotations.get("organism", "unknown_species").strip()
								for feat in rec.features:
									if feat.type != "CDS":
										continue
									gene_name = feat.qualifiers.get("gene", [None])[0]
									prot_id = feat.qualifiers.get("protein_id", [None])[0]
									locus_tag = feat.qualifiers.get("locus_tag", [None])[0]
									if not gene_name:
										continue
									# 存入索引
									key_org_gene = (org_name.lower(), gene_name.lower())
									if key_org_gene not in gb_gene_index:
										gb_gene_index[key_org_gene] = []
									gb_gene_index[key_org_gene].append( (rec, feat) )
									# protein_id / locus_tag 映射，用于和OrthoFinder序列匹配
									if prot_id:
										gb_cds_id_map[prot_id] = (org_name, gene_name, rec, feat)
									if locus_tag:
										gb_cds_id_map[locus_tag] = (org_name, gene_name, rec, feat)
						except Exception as e:
							logger.error(f"Parse GB file {gb_path} failed: {str(e)}")
				# ===================================================================

				logger.info('Running OrthoFinder on protein sequences: %s', seqdir)
				of_out_dir = os.path.join(self.dbdir, "of_results")

				# 删除旧目录：if只做删除，cmd不能嵌套在这个if内部！！【修复Bug2缩进】
				if os.path.exists(of_out_dir):
					rmdirs(of_out_dir)

				# 👇 cmd、subprocess 挪到if外面，无论目录存在与否都执行
				cmd = [
					sys.executable,
					orthofinder_exe,
					"-f", seqdir,
					"-o", of_out_dir,
					"-S", "diamond",
					"-t", str(self.ncpu)
					]
				sub_env = get_of_subenv()
				proc = None
				try:
					proc = subprocess.run(
						cmd,
						env=sub_env,
						cwd=self.dbdir,
						stdout=sys.stdout,
						stderr=sys.stderr,
						text=True
						)
				except Exception as e:
					raise RuntimeError(f"OrthoFinder execute error: {e}")
				if proc is not None:
					if proc.returncode != 0:
						raise RuntimeError(f"OrthoFinder failed, ret={proc.returncode}")
				logger.info("OrthoFinder run finished.")
				logger.info(f"[DEBUG‑OF] of_out_dir = {of_out_dir}")

				res_dirs = []
				res_dirs.extend(glob.glob(os.path.join(of_out_dir, "Results_*")))
				res_dirs.extend(glob.glob(os.path.join(of_out_dir, "OrthoFinder", "Results_*")))
				if not res_dirs:
					logger.warning("OrthoFinder output missing, skip ortholog/hmm/prfl for this feature type.")
					continue
				result_dir = res_dirs[0]

				og_tsv_v3 = os.path.join(result_dir, "Orthogroups", "Orthogroups.tsv")
				og_tsv_v2 = os.path.join(result_dir, "Orthogroups.tsv")
				if os.path.exists(og_tsv_v3):
					og_tsv = og_tsv_v3
				elif os.path.exists(og_tsv_v2):
					og_tsv = og_tsv_v2
				else:
					logger.warning("Orthogroups.tsv not found (OF2/OF3 checked), skip ortholog/hmm/prfl")
					continue
				of_msa_dir = os.path.join(result_dir, "MultipleSequenceAlignments")
				logger.info("Parsing orthogroups file: %s", og_tsv)
			# ============ OF调用块 END ============

				#import csv
				#from Bio import SeqIO
				#from io import StringIO
				# -------------------------- 配置：线粒体核心基因白名单 --------------------------
				mt_core_genes = {
					"atp1","atp4","atp6","atp8","atp9",
					"nad1","nad2","nad3","nad4","nad4L","nad5","nad6","nad7","nad9",
					"cob","cox1","cox2","cox3",
					"ccmB","ccmC","ccmFc","ccmFn","matR","mttB",
					"rpl2","rpl5","rpl10","rpl16",
					"rps1","rps3","rps4","rps7","rps10","rps12","rps13","rps14","rps19",
					"sdh3","sdh4"
				}
				# 次优先：其他标准基因，最低优先级：orf开头
				def pick_best_gene_name(name_set):
					# 优先匹配白名单核心基因
					hit = name_set & mt_core_genes
					if hit:
						return sorted(list(hit))[0]
					# 排除orf类，剩下的标准基因
					non_orf = [n for n in name_set if not n.startswith("orf")]
					if non_orf:
						return sorted(non_orf)[0]
					# 兜底，只能选orf
					return sorted(list(name_set))[0]

				# -------------------------- Step1：解析Orthogroups.tsv，构建orthogroups字典 --------------------------
				orthogroups = dict()   # key: OG编号, value: list[gene_id]
				og_gene_names = dict() # key: OG编号, value: set(所有该OG内基因名)

				with open(og_tsv, "r", encoding="utf-8") as f:
					reader = csv.reader(f, delimiter="\t")
					header = next(reader)
					# header[0]是Orthogroup，后面每一列是物种名
					for row in reader:
						og_id = row[0]
						gene_ids = []
						for cell in row[1:]:
							if cell.strip() == "":
								continue
							gene_ids.extend([g.strip() for g in cell.split(",")])
						orthogroups[og_id] = gene_ids
						# 从全局cds_dict拿到每个gene_id对应的gene_name，收集到集合
						name_set = set()
						for gid in gene_ids:
							#if gid in self.cds_dict:
								#gname = self.cds_dict[gid].annotations.get("gene", None)
							if gid in self.prot_id2gene:
								gname = self.prot_id2gene[gid]
								if gname:
									name_set.add(gname)
						og_gene_names[og_id] = name_set

				# -------------------------- Step2：构建 og2gene 映射（白名单选基因名） --------------------------
				og2gene = dict()
				for og_name, name_set in og_gene_names.items():
					if len(name_set) ==0:
						logger.warning(f"{og_name}: no gene name found, skip og2gene mapping")
						continue
					best_gene = pick_best_gene_name(name_set)
					og2gene[og_name] = best_gene
					logger.debug(f"{og_name} mapped to gene {best_gene}, raw names: {name_set}")

				# ----------------在这里粘贴昨天给到的sto→hmmbuild代码----------------

				of_msa_dir = os.path.join(result_dir, "MultipleSequenceAlignments")
				logger.info(f"Reading OF multiple alignments from: {of_msa_dir}")

				# =========【在这里整块插入 og2gene 映射代码】=========
				og_tsv = os.path.join(result_dir, "Orthogroups", "Orthogroups.tsv")
				#og2gene = dict()
				if not os.path.exists(og_tsv):
					raise FileNotFoundError(
						f"Orthogroups.tsv not found: {og_tsv}\n"
						"OrthoFinder output incomplete, cannot build OG‑gene mapping."
					)
				with open(og_tsv, newline='') as fh:
					reader = csv.reader(fh, delimiter='\t')
					next(reader)  # skip header line
					for row in reader:
						og_name = row[0]
						seq_ids = [s.strip() for s in row[1:] if s.strip()]
						if not seq_ids:
							continue
						# 遍历本组所有序列ID，提取基因名，校验是否统一
						gene_candidates = set()
						for sid in seq_ids:
							gene_candidates.add(sid.split("_")[0])
						real_gene = list(gene_candidates)[0]
						if len(gene_candidates) > 1:
							logger.warning(f"OG {og_name} 存在多个基因名称 {gene_candidates}, 选用 {real_gene}")
						og2gene[og_name] = real_gene
				logger.info(f"Build OG-to-gene mapping done, total entries: {len(og2gene)}")

				# 新增：初始化agst输出目录
				agst_dir = os.path.join(self.trndir, "agst")
				mkdirs(agst_dir)
				if run_agst_augustus:
					for og_name, real_gene in og2gene.items():
						msa_fa = os.path.join(self.msadir, f"{og_name}.fa")
						if not os.path.exists(msa_fa):
							logger.warning(f"OG {og_name} msa file missing, skip AGST for {real_gene}")
							continue
						# 改用 self.dbname，增加 mt- 前缀对齐旧库命名规范
						# 根据细胞器ogtype自动切换 mt-/pt- 前缀，对齐旧库命名规范
						logger.info(f"DEBUG ogtype value: [{self.ogtype}], db_prefix: [{db_prefix}]")

						if self.ogtype == "mt":
							species_name = f"mt-{self.dbname}-{real_gene}"
						elif self.ogtype == "pt":
							species_name = f"pt-{self.dbname}-{real_gene}"
						else:
							raise ValueError(f"ogtype [{self.ogtype}] 仅支持 mt / pt")
						#species_name = f"mt-{self.dbname}-{real_gene}"
						species_cfg_dir = os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"], "species", species_name)

						# 【新增】清理旧的species目录，避免缓存干扰
						if os.path.exists(species_cfg_dir):
							logger.info(f"Remove old augustus species dir: {species_cfg_dir}")
							shutil.rmtree(species_cfg_dir)
						train_tmp = os.path.join(agst_dir, species_name)
						mkdirs(train_tmp)
						# step1 new_species.pl
						cmd_newsp = [
							"new_species.pl",
							"--species", species_name,
							"--AUGUSTUS_CONFIG_PATH", os.environ["AUGUSTUS_CONFIG_PATH"]
						]
						logger.info(f"AGST new_species.pl gene={real_gene}, species={species_name}")
						ret1 = subprocess.run(cmd_newsp, capture_output=True, text=True)
						if ret1.returncode != 0:
							logger.error(f"new_species.pl failed for {real_gene}\nSTDERR:\n{ret1.stderr}")
							shutil.rmtree(train_tmp)
							continue
						# ========= 拼装训练用Genbank =========
						query_gene = real_gene.lower()
						collect_records = []
						for (org_lower, gname_lower), rec_list in gb_gene_index.items():
							if gname_lower != query_gene:
								continue
							for parent_full_rec, cds_feat in rec_list:
								cds_seq = cds_feat.location.extract(parent_full_rec.seq)
								# 新增：过滤空序列
								if len(cds_seq) == 0:
									logger.warning(f"Skipping {real_gene}: extracted cds_seq length is zero!")
									continue
								cds_len = len(cds_seq)
								new_rec = SeqRecord(
									cds_seq,
									id=f"{org_lower}_{real_gene}",
									name=real_gene,
									description=f"{org_lower} {real_gene} CDS for AGST training"
								)
								new_cds_loc = FeatureLocation(0, cds_len, strand=cds_feat.location.strand)
								new_cds_feat = SeqFeature(new_cds_loc, type="CDS", qualifiers=cds_feat.qualifiers.copy())
								new_rec.features = [new_cds_feat]
								# 被子植物线粒体使用 table1 标准密码子
								new_rec.annotations["codon_table"] = 1
								collect_records.append(new_rec)
						if len(collect_records) == 0:
							logger.warning(f"AGST: no CDS found for gene {real_gene}, skip")
							shutil.rmtree(train_tmp)
							continue
						train_gb = os.path.join(self.msadir, f"{og_name}.gb")
						with open(train_gb, "w") as fgb:
							SeqIO.write(collect_records, fgb, "genbank")
						logger.info(f"AGST build training gb for {real_gene}, records={len(collect_records)} file={train_gb}")
						# step2 etraining
						cmd_etrain = [
							"etraining",
							"--species", species_name,
							train_gb
						]
						logger.info(f"AGST run etraining for gene {real_gene}")
						ret2 = subprocess.run(cmd_etrain, capture_output=True, text=True)
						if ret2.returncode != 0:
							logger.error(f"etraining failed for {real_gene}\nSTDERR:\n{ret2.stderr}")
							if os.path.exists(species_cfg_dir):
								shutil.rmtree(species_cfg_dir)
							shutil.rmtree(train_tmp)
							continue
						# 打包模型
						shutil.copytree(species_cfg_dir, os.path.join(train_tmp, species_name))
						tgz_file = os.path.join(agst_dir, f"{species_name}.tgz")
						with tarfile.open(tgz_file, "w:gz") as tar:
							tar.add(os.path.join(train_tmp, species_name), arcname=species_name)
						logger.info(f"AGST successfully generated gene model: {tgz_file}")

						# 清理临时
						if os.path.exists(species_cfg_dir):
							shutil.rmtree(species_cfg_dir)
						shutil.rmtree(train_tmp)

				# ========== 复制OrthoFinder MSA文件到msadir ==========
				# 自动匹配 Results_* 目录，避免写死日期
				result_subdirs = glob.glob(os.path.join(of_out_dir, "Results_*"))
				if len(result_subdirs) != 1:
					raise FileNotFoundError(f"Found {len(result_subdirs)} Results_* dir in {of_out_dir}, expected 1")
				of_result_dir = result_subdirs[0]
				of_msa_root = os.path.join(of_result_dir, "MultipleSequenceAlignments")
				logger.info(f"Copying OF MSA files from {of_msa_root} to {self.msadir}")

				copy_count = 0
				for og_id in og2gene.keys():
					src_msa = os.path.join(of_msa_root, f"{og_id}.fa")
					dst_msa = os.path.join(self.msadir, f"{og_id}.fa")
					if os.path.exists(src_msa):
						try:
							shutil.copyfile(src_msa, dst_msa)
							copy_count += 1
						except Exception as e:
							logger.warning(f"Copy MSA failed for {og_id}: {e}")
					else:
						logger.warning(f"MSA file not found for {og_id}: {src_msa}")
				logger.info(f"MSA copy done: {copy_count}/{len(og2gene)} files copied to msas/")
				# ========== MSA拷贝结束 ==========

				# ========== 原生AGST：new_species.pl + etraining 训练augustus物种模型（复刻老师库） ==========
				agst_dir = os.path.join(self.dbdir, "agst")
				mkdirs(agst_dir)
				db_prefix = os.path.basename(self.dbdir)
				msa_dir = self.msadir

				for og_name, real_gene in og2gene.items():
					#msa_fa = os.path.join(msa_dir, f"{og_name}.fasta")
					msa_fa = os.path.join(msa_dir, f"{og_name}.fa")
					if not os.path.exists(msa_fa):
						msa_fa = os.path.join(msa_dir, f"{og_name}.fasta")
					if not os.path.exists(msa_fa):
						logger.warning(f"OG {og_name} -> gene {real_gene}: msa file {msa_fa} missing, skip agst train")
						continue

				#	species_name = f"{db_prefix}-{real_gene}"
					# 新增：从拼接字符串里提取细胞器标记，仅用于species_name拼接
					work_ogtype = self.ogtype.split("-")[0]
					
					if work_ogtype == "mt":
						species_name = f"mt-{db_prefix}-{real_gene}"
					elif work_ogtype == "pt":
						species_name = f"pt-{db_prefix}-{real_gene}"
					else:
						raise ValueError(f"ogtype [{self.ogtype}] 仅支持 mt / pt")


					# 物种模型目录，放在 augustus config 下面！！！重点！
					species_cfg_dir = os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"], "species", species_name)
					# 临时工作目录
					train_tmp = os.path.join(agst_dir, species_name)
					mkdirs(train_tmp)

    				# step1: new_species.pl 新建物种模板
					cmd_newsp = [
						"new_species.pl",
						"--species=" + species_name,
						f"--AUGUSTUS_CONFIG_PATH={os.environ['AUGUSTUS_CONFIG_PATH']}"
					]
					logger.info(f"AGST step1 new_species.pl for {real_gene}: {' '.join(cmd_newsp)}")
					ret1 = subprocess.run(cmd_newsp, capture_output=True, text=True)
					if ret1.returncode != 0:
						logger.error(f"new_species.pl failed {real_gene}, stderr:\n{ret1.stderr[:800]}")
						shutil.rmtree(train_tmp)
						continue

					# ========== 拼装单基因GenBank文件 ==========
					# ========= 拼装训练用Genbank =========
					query_gene = real_gene.lower()
					collect_records = []
					for (org_lower, gname_lower), rec_list in gb_gene_index.items():
						if gname_lower != query_gene:
							continue
					# ========= 拼装训练用Genbank =========
					query_gene = real_gene.lower()
					collect_records = []
					for (org_lower, gname_lower), rec_list in gb_gene_index.items():
						if gname_lower != query_gene:
							continue
						for parent_full_rec, cds_feat in rec_list:
							cds_seq = cds_feat.location.extract(parent_full_rec.seq)
							cds_len = len(cds_seq)
							new_rec = SeqRecord(
								cds_seq,
								id=f"{org_lower}_{real_gene}",
								name=real_gene,
								description=f"{org_lower} {real_gene} CDS for AGST training"
							)
							# 新本地坐标：0到cds序列全长
							new_loc = FeatureLocation(0, cds_len, strand=cds_feat.location.strand)
							# source feature（Augustus etraining 强制需要）
							source_feat = SeqFeature(new_loc, type="source", qualifiers={"organism":[org_lower]})
							# gene feature
							new_gene_feat = SeqFeature(new_loc, type="gene", qualifiers={"gene": [real_gene]})
							# CDS feature
							new_cds_feat = SeqFeature(new_loc, type="CDS", qualifiers={"gene":[real_gene]})
							new_rec.features = [source_feat, new_gene_feat, new_cds_feat]
							new_rec.annotations["codon_table"] = 1
							new_rec.annotations["molecule_type"] = "DNA"
							collect_records.append(new_rec)

					if len(collect_records) == 0:
						logger.warning(f"AGST: no CDS found for gene {real_gene}, skip")
						shutil.rmtree(train_tmp)
						continue
					train_gb = os.path.join(self.msadir, f"{og_name}.gb")
					with open(train_gb, "w") as fgb:
						SeqIO.write(collect_records, fgb, "genbank")
					logger.info(f"AGST build training gb for {real_gene}, records={len(collect_records)} file={train_gb}")
					# =========新增：少于3条序列直接跳过etraining=========
					if len(collect_records) < 3:
						logger.warning(f"AGST skip {real_gene}: only {len(collect_records)} records, need >=3 for etraining")
						shutil.rmtree(train_tmp)
						continue
					# ========== 拼装代码结束 ==========

					# step2: etraining 训练生成pbl概率文件
					# step2 etraining
					cmd_etrain = [
						"etraining",
						f"--species={species_name}",
						train_gb
					]
					logger.info(f"AGST run etraining for gene {real_gene}, training_gb={train_gb}")
					ret2 = subprocess.run(cmd_etrain, capture_output=True, text=True)
					if ret2.returncode != 0:
						logger.error(f"==== etraining FAILED for {real_gene} ====")
						logger.error(f"STDOUT:\n{ret2.stdout}")
						logger.error(f"STDERR:\n{ret2.stderr}")
						if os.path.exists(species_cfg_dir):
							shutil.rmtree(species_cfg_dir)
						shutil.rmtree(train_tmp)
						continue


					# step3: 把训练产物从config/species拷贝到agst临时目录，打包
					shutil.copytree(species_cfg_dir, os.path.join(train_tmp, species_name))
					tgz_out = os.path.join(agst_dir, f"{species_name}.tgz")
					cmd_tar = ["tar", "-zcf", tgz_out, "-C", train_tmp, species_name]
					ret_tar = subprocess.run(cmd_tar, capture_output=True, text=True)
					if ret_tar.returncode !=0:
						logger.warning(f"tar pack agst model failed {species_name}, stderr:{ret_tar.stderr[:300]}")
					# 清理临时文件
					shutil.rmtree(train_tmp)
					# 【重要】删除config目录下的物种文件，避免config堆积大量物种
					if os.path.exists(species_cfg_dir):
						shutil.rmtree(species_cfg_dir)
					logger.info(f"AGST model success: {species_name}.tgz")

				logger.info("All augustus species models packed into agst/")

				# ======================================================

				fa_files = glob.glob(os.path.join(of_msa_dir, "*.fa"))
				logger.info(f"Found total orthogroup fa alignments: {len(fa_files)}")
				# prfl_entries = []

				# 循环外：初始化prfl目录、perl、定位msa2prfl.pl
				prfl_dir = os.path.join(self.dbdir, "prfl")
				if not os.path.exists(prfl_dir):
					os.makedirs(prfl_dir)
				perl_bin = shutil.which("perl")
				run_env = os.environ.copy()

				aug_bin_path = shutil.which("augustus")
				if aug_bin_path:
					aug_base = os.path.dirname(os.path.dirname(aug_bin_path))
					msa2prfl = os.path.join(aug_base, "scripts", "msa2prfl.pl")
				else:
					msa2prfl = None
					logger.warning("augustus not found, skip generating *.prfl files, only model.lst written")

				for fa_path in fa_files:
					base = os.path.basename(fa_path)
					og_name = os.path.splitext(base)[0]
					hmm_out = os.path.join(self.hmmdir, f"{og_name}.hmm")
					fa_dst = os.path.join(self.msadir, base)
					shutil.copyfile(fa_path, fa_dst)
					cmd_hmmbuild = ["hmmbuild", hmm_out, fa_path]
					ret = subprocess.run(cmd_hmmbuild, capture_output=True, text=True)
					if ret.returncode != 0:
						logger.warning(f"hmmbuild failed for {og_name}, skip. stderr:{ret.stderr[:200]}")
						continue

					if og_name in og2gene:
						real_gene = og2gene[og_name]
						gene_named_hmm = os.path.join(self.hmmdir, f"{real_gene}.hmm")
						try:
							os.link(hmm_out, gene_named_hmm)
						except OSError:
							shutil.copyfile(hmm_out, gene_named_hmm)
						logger.debug(f"Map {og_name} → {real_gene}, create {real_gene}.hmm")
						# 写入真实基因名到prfl_entries，最终进model.lst
						prfl_entries.append(real_gene)
					else:
						logger.warning(f"OG {og_name} missing in og2gene mapping, skip gene‑named hmm")

					# 生成原始og_name.prfl
					if msa2prfl and os.path.exists(msa2prfl) and perl_bin:
						prfl_out = os.path.join(prfl_dir, f"{og_name}.prfl")
						cmd_prfl = [
							perl_bin,
							msa2prfl,
							"--setname", og_name,
							fa_path
						]
						logger.info(f"DEBUG msa2prfl.pl cmd: {cmd_prfl}")
						ret_p = subprocess.run(cmd_prfl, env=run_env, capture_output=True, text=True)
						out_content = ret_p.stdout if ret_p.stdout is not None else ""
						#with open(prfl_out, "w") as fh:
						#	fh.write(out_content)
						# 空输出直接跳过，不生成0字节prfl文件
						if len(out_content.strip()) == 0:
							logger.warning(f"msa2prfl.pl output empty for {og_name}, skip write prfl")
							continue
						with open(prfl_out, "w") as fh:
							fh.write(out_content)

						stderr_text = ret_p.stderr
						if stderr_text:
							logger.warning(f"msa2prfl.pl failed for {og_name}, ret={ret_p.returncode}, stderr:{stderr_text}")
						if not os.path.isfile(prfl_out):
							logger.warning(f"prfl file not generated: {prfl_out}")
						elif os.path.getsize(prfl_out) == 0:
							logger.warning(f"prfl file is empty(No blocks): {prfl_out}")
						else:
							# 生成 real_gene.prfl 软链接
							if og_name in og2gene:
								real_gene = og2gene[og_name]
								gene_prfl = os.path.join(prfl_dir, f"{real_gene}.prfl")
								try:
									os.link(prfl_out, gene_prfl)
								except OSError:
									shutil.copyfile(prfl_out, gene_prfl)
								logger.debug(f"Map prfl {og_name} → {real_gene}: create {real_gene}.prfl")

								# ========== 新增：修改prfl内部[name]字段为真实基因名 ==========
								try:
									with open(gene_prfl, "r") as f:
										lines = f.readlines()
										in_name_block = False
									for idx, line in enumerate(lines):
										stripped = line.strip()
										if stripped == "[name]":
											in_name_block = True
											continue
										if in_name_block:
											# 只替换[name]紧跟的下一行，其余内容不动
											lines[idx] = f"{real_gene}\n"
											break
									with open(gene_prfl, "w") as f:
										f.writelines(lines)
									logger.debug(f"Fix prfl inner name: {gene_prfl}, set inner name to {real_gene}")
								except Exception as e:
									logger.warning(f"Failed to modify inner name of {gene_prfl}: {str(e)}")
								# ================================================================
							# =================================================================

							# ===================== 追加 seqs fa 映射 =====================
							if og_name in og2gene:
								real_gene = og2gene[og_name]
								gene_fa = os.path.join(self.seqdir, f"{real_gene}.fa")
								# 去除MSA比对gap，输出exonerate用的fa
								with open(gene_fa, "w") as fw:
									for rec in SeqIO.parse(fa_path, "fasta"):
										s_clean = str(rec.seq).replace("-", "")
										fw.write(f">{rec.id}\n{s_clean}\n")
								logger.debug(f"Map seq fa {og_name} → {real_gene}: create {real_gene}.fa (stripped gap)")

								# ==========================================================

				# for循环结束后写model.lst-暂时不再使用
				#prfl_file = os.path.join(self.pfldir, "model.lst")
				#with open(prfl_file, "w") as f:
				#	for name in prfl_entries:
				#		f.write(f"{name}\n")

				#logger.info(f"Success build {len(prfl_entries)} HMM profiles.")
				logger.info('==== LEAVE protein OF branch SUCCESS ====')
				# ----------------粘贴结束----------------

			# ==========【在这里整块粘贴RNA代码】==========
			#elif feat_type != "protein":
			else:
				rna_groups = dict()
				rna_collection = dict()
				rna_built = 0
				for rec in records:
					sp = rec.organism
					for ftype in ("rRNA", "tRNA"):
						for ft in rec.get_features(ftype):
							qual = ft.qualifiers
							if "pseudo" in qual or "pseudogene" in qual:
								continue
							gene_raw = qual.get("gene", [None])[0]
							if not gene_raw:
								continue
							product_raw = qual.get("product", [""])[0]
							gene_norm = re.sub(r"\((.+?)\)", r"-\1", gene_raw).strip()
							loc = ft.location
							s = int(loc.start)
							e = int(loc.end)
							seq_frag = rec.seq[s:e]
							if loc.strand == -1:
								seq_frag = seq_frag.reverse_complement()
							seq_str = str(seq_frag)
							if len(seq_str) < 30:
								continue
							if gene_norm not in rna_collection:
								rna_collection[gene_norm] = {
									"name": gene_norm,
									"product": product_raw,
									"seq_type": ftype,
									"exon_count": 1,
									"trans_splicing": "false",
									"gene_count": 0,
									"species_set": set(),
									"names": set(),
									"products": set(),
								}
							entry = rna_collection[gene_norm]
							entry["gene_count"] += 1
							entry["species_set"].add(sp)
							entry["names"].add(gene_raw)
							entry["products"].add(product_raw)
							if gene_norm not in rna_groups:
								rna_groups[gene_norm] = []
							rna_groups[gene_norm].append((sp, seq_str))
				if len(rna_groups) == 0:
					logger.warning("No valid rRNA/tRNA features for building RNA HMM, skip RNA profile")
					rna_info_path = os.path.join(self.dbdir, "rna.info")
					with open(rna_info_path, "w", encoding="utf-8") as fout:
						fout.write("#id\tname\tproduct\tseq_type\texon_count\ttrans_splicing\tgene_count\tspecies_count\tnames\tproducts\tspecies\n")
					continue
				mafft_bin = shutil.which("mafft")
				hmmbuild_bin = shutil.which("hmmbuild")
				if not (mafft_bin and hmmbuild_bin):
					logger.warning("Missing mafft/hmmbuild, skip RNA profile build")
				run_env = os.environ.copy()
				tmp_rna_root = tempfile.TemporaryDirectory(dir=self.tmpdir)
				tmp_rna_dir = tmp_rna_root.name

				for rna_gene, seq_list in rna_groups.items():
					gene_info = rna_collection[rna_gene]
					hmm_out = os.path.join(self.hmmdir, f"{rna_gene}.hmm")
					# =====注释掉rRNA占位分支=====
					# if gene_info["seq_type"] == "rRNA":
					# 	logger.debug(f"rRNA {rna_gene}: create dummy placeholder hmm")
					# 	with open(hmm_out, "w", encoding="utf-8") as fh:
					# 		fh.write("HMMER3/f [dummy placeholder for rRNA, not used for hmmsearch]\n")
					# 	continue
					# ===========================
					fa_tmp = os.path.join(tmp_rna_dir, f"{rna_gene}.fa")
					with open(fa_tmp, "w") as fout:
						for sp_name, seqtext in seq_list:
							fout.write(f">{sp_name}\n{seqtext}\n")
					aln_tmp = os.path.join(tmp_rna_dir, f"{rna_gene}.aln")
					cmd_mafft = [mafft_bin, "--auto", fa_tmp]
					ret_mafft = subprocess.run(cmd_mafft, env=run_env, stdout=open(aln_tmp,"w"), stderr=subprocess.PIPE, text=True)
					if ret_mafft.returncode != 0:
						logger.warning(f"mafft fail for {rna_gene}, skip")
						continue
					cmd_hmmbuild = [hmmbuild_bin, "--dna", hmm_out, aln_tmp]
					ret_hmm = subprocess.run(cmd_hmmbuild, env=run_env, capture_output=True, text=True)
					if ret_hmm.returncode !=0:
						logger.warning(f"hmmbuild fail {rna_gene}")
						continue
					rna_built += 1
					# ==========新增代码 开始==========
					shutil.copy(fa_tmp, os.path.join(self.seqdir, f"{rna_gene}.fa"))
					# ==========新增代码 结束==========
					logger.debug(f"RNA DNA‑HMM built: {rna_gene}, type:{gene_info['seq_type']}")

				tmp_rna_root.cleanup()
				logger.info(f"RNA profile done, RNA models built: {rna_built}")
				def clean_text(s):
					return str(s).replace(",", ";").replace("\t", " ").strip()
				rna_info_path = os.path.join(self.dbdir, "rna.info")
				with open(rna_info_path, "w", encoding="utf-8") as fout:
					header = "#id\tname\tproduct\tseq_type\texon_count\ttrans_splicing\tgene_count\tspecies_count\tnames\tproducts\tspecies\n"
					fout.write(header)
					for gid, item in rna_collection.items():
						species_list = clean_text(";".join(sorted(item["species_set"])))
						names_str = clean_text(";".join(sorted(item["names"])))
						products_str = clean_text(";".join(sorted(item["products"])))
						product_clean = clean_text(item["product"])
						line = (
							f"{gid}\t{item['name']}\t{product_clean}\t{item['seq_type']}\t"
							f"{item['exon_count']}\t{item['trans_splicing']}\t{item['gene_count']}\t"
							f"{len(item['species_set'])}\t{names_str}\t{products_str}\t{species_list}\n"
						)
						fout.write(line)


			# ==========RNA代码块结束==========


			ntaxa = len(set(taxa))	# some taxa have >1 records
			if out_info:
				logger.info('{} taxa: {}'.format(ntaxa, sorted(set(taxa))))
				# species_info
				with open(self.species_file, 'w') as fout:
					SpeciesInfoLine().write(fout)
					for rc in sorted(records, key=lambda x: x.organism):
						taxonomy = ';'.join(rc.taxonomy)
						line = [rc.id, rc.organism, rc.name_count, rc.cds_count, rc.trn_count,
								rc.rrn_count, rc.length, rc.GC, taxonomy]
						SpeciesInfoLine(line).write(fout)
				# taxonomy_info
				taxonomy = [';'.join(rc.taxonomy) for rc in records]
				taxon, taxon_count = self.get_most_common(taxonomy)
				logger.info('get taxon `{}` from {}'.format(taxon, taxon_count))

				tables = []
				for rc in records:
					for ft in rc.get_features("CDS"):
						tables.append(ft.transl_table)
				transl_table, table_count = self.get_most_common(tables)
				logger.info('get transl_table `{}` from {}'.format(transl_table, table_count))
				line = [self.organ, self.ogtype, transl_table, taxon]
				with open(self.taxonomy_file, 'w') as fout:
					TaxonomyInfoLine().write(fout)
					TaxonomyInfoLine(line).write(fout)


			#orf_candidates = glob.glob('{}/OrthoFinder/*/'.format(seqdir))
			#if not orf_candidates:
				#logger.warning("OrthoFinder output missing, skip ortholog/hmm/prfl for this feature type.")
				#continue

			# ==========【CDS 聚合：生成 11 列 cds.info】==========
			def clean_text(s):
				return str(s).replace(",", ";").replace("\t", " ").strip()
			cds_collection = dict()
			for rc in records:
				spec = rc.organism.replace(" ", "_")
				for feat in rc.get_features("CDS"):
					try:
						gene_raw = feat.qualifiers.get("gene", [None])[0]
						if not gene_raw:
							continue
						product_raw = feat.qualifiers.get("product", [""])[0]
						gene_norm = gene_raw.strip()
						if gene_norm not in cds_collection:
							cds_collection[gene_norm] = {
								"name": gene_norm,
								"product": product_raw,
								"seq_type": "CDS",
								"exon_count": 1,
								"trans_splicing": "false",
								"gene_count": 0,
								"species_set": set(),
								"names": set(),
								"products": set(),
							}
						entry = cds_collection[gene_norm]
						entry["gene_count"] += 1
						entry["species_set"].add(spec)
						entry["names"].add(gene_raw)
						entry["products"].add(product_raw)
					except Exception as e:
						print(f"[DBG-cds-collect] skip CDS feature {e}", file=sys.stderr)
						continue
			cds_info_path = os.path.join(self.dbdir, "cds.info")
			with open(cds_info_path, "w", encoding="utf-8") as fout:
				header = "#id\tname\tproduct\tseq_type\texon_count\ttrans_splicing\tgene_count\tspecies_count\tnames\tproducts\tspecies\n"
				fout.write(header)
				for gid, item in cds_collection.items():
					species_list = clean_text(";".join(sorted(item["species_set"])))
					names_str = clean_text(";".join(sorted(item["names"])))
					products_str = clean_text(";".join(sorted(item["products"])))
					product_clean = clean_text(item["product"])
					line = (
						f"{gid}\t{item['name']}\t{product_clean}\t{item['seq_type']}\t"
						f"{item['exon_count']}\t{item['trans_splicing']}\t{item['gene_count']}\t"
						f"{len(item['species_set'])}\t{names_str}\t{products_str}\t{species_list}\n"
					)
					fout.write(line)
			logger.info(f"CDS info written: {cds_info_path}, total genes: {len(cds_collection)}")
			# ==========CDS 聚合结束==========
			self.export_ref_sequences(records)

		# =========退出feat_type循环，往下==========

		# 统一输出model.lst（蛋白 + RNA全部）
		prfl_file = os.path.join(self.pfldir, "model.lst")
		with open(prfl_file, "w") as f:
			for name in prfl_entries:
				f.write(f"{name}\n")
		logger.info(f"Total built HMM profiles(protein+RNA): {len(prfl_entries)}")

		logger.info('pruning taxa')

		records = self.prune(records)

		self.count_rank(records)

		# tar
		self.tgz_dirs(self.seqdir, self.msadir, self.hmmdir, self.pfldir, self.tnddir)
		if self.cleanup:
			#logger.info('cleaning —— DEBUG SKIP REMOVE tmpdir')
			logger.info("cleaning temporary directory")
			rmdirs(self.tmpdir)
		logger.info('build completed')
	def check_unique(self, files):
		if len(files) != len(set(files)):
			logger.error('{}..are not unique'.format(files[0]))

	def parse_orthologs(self, orf, features, outdir, fout):
		d_seqs = {feat.id: feat.seq for feat in features}
		d_names = {feat.id: feat.gene for feat in features}
		#print >>sys.stderr, d_names
		d_products = {feat.id: feat.product for feat in features}
		d_types = {feat.id: feat.type for feat in features}
		d_nexon = {feat.id: feat.nexon for feat in features}
		d_trans_splicing = {feat.id: feat.trans_splicing for feat in features}

		logger.info('using min_count={} to get gene name'.format(self.min_count))
		d_genes = {}
		lines = []
		gene_names, outseqs, groups = [], [], []
		if len(list(orf.get_orthologs_cluster())) == 0:
			orthologs = orf.get_orthogroups()
			logger.info('using orthogroups')
		else:
			orthologs = orf.get_orthologs_cluster()
			logger.info('using orthologues groups')
#		for i, group in enumerate(orf.get_orthologs_cluster()):		# group by feat.id
#		for i, group in enumerate(orf.get_orthogroups()):
		for i, group in enumerate(orthologs):
		#	print group, d_names
			gene_name, name_count = self.get_gene_name(group, d_names,)
			tmpfix = '{}/{}.para'.format(outdir, gene_name)
			group = self.filter_paralogs(group, d_seqs, tmpfix)
			# name
			gene_name, name_count = self.get_gene_name(group, d_names,)
			gene_name = NameInfo.get_trn_name(gene_name, name_count)
			logger.info('get gene name `{}` from {}'.format(gene_name, name_count))
			max_count = max(name_count.values()) #name_count[gene_name]
			#max_count = sum(name_count.values())
			if max_count < self.min_count:
				continue
			logger.info('using gene name `{}`'.format(gene_name))
			gene_id = gene_name
			if gene_id in d_genes:   # for duplicated gene names
				gene_id = '{}-{}'.format(gene_name, i)
			d_genes[gene_id] = (gene_name, group)

			# product
			product, product_count = self.get_gene_name(group, d_products)
#			logger.info('get product `{}` from {}'.format(product, product_count))
			top_products = product_count.most_common(25)
			top_products = ','.join(['{}:{}'.format(name, count) \
							for name, count in top_products])
			top_names = name_count.most_common(25)
			top_names = ','.join(['{}:{}'.format(name, count) \
							for name, count in top_names if len(name)<18])
			# count
			gene_count = len(group)
			species = [gene.split('|')[0] for gene in group]
			species_count = Counter(species)
			species = ','.join(['{}:{}'.format(name, count) \
								for name, count in sorted(species_count.items())])
			species_count = len(species_count)
			# exon number and trans_splicing
			nexon, nexon_count = self.get_gene_name(group, d_nexon)
			trans_splicing, ts_count = self.get_gene_name(group, d_trans_splicing)
#			logger.info('get exon number `{}` from {}'.format(nexon, nexon_count))
#			logger.info('get trans_splicing `{}` from {}'.format(trans_splicing, ts_count))

			# type
			seqtype, _ = self.get_gene_name(group, d_types)
			# write info
			line = [gene_id, gene_name, product, seqtype, 
					nexon, trans_splicing,
					gene_count, species_count, top_names, 
					top_products, species]
			lines += [line]
			# write sequences
			outseq = self.get_seqfile(gene_id)
		#	print >>sys.stderr, gene_id, gene_name, group, outseq
			self.write_orthologs(group, d_seqs, outseq)	# by group
			outseqs += [outseq]
			gene_names += [gene_id]
			groups += [group]

		GeneInfoLine().write(fout)
		for line in sorted(lines):
			GeneInfoLine(line).write(fout)
		return gene_names, outseqs, groups

	def train_augustus(self, train_set, species, gene_number, outdir, ncpu, transl_table=1):
		gene_number = old_div(gene_number,3)
		kfold = 8 if ncpu > 8 else ncpu
		spdir = '{}/species/'.format(AUGUSTUS_CONFIG_PATH,)
		cmd = '''# train
(cd {spdir} && rm {species} -r)
randomSplit.pl {train_set} {gene_number}
new_species.pl --species={species}
etraining --species={species} {train_set}.train --translation_table={transl_table}
augustus --species={species} {train_set}.test --translation_table={transl_table} &> {train_set}.train_test1.out
# optimize
#optimize_augustus.pl --species={species} --cpus={ncpu} --kfold={kfold} --rounds=5 {train_set}.train && rm -r tmp_opt_{species}
#etraining --species={species} {train_set}.train
#augustus --species={species} {train_set}.test &> {train_set}.train_test2.out
(cd {spdir} && tar cvzf {outdir}/{species}.tgz {species})
'''.format(train_set=train_set, gene_number=gene_number, 
			species=species, ncpu=ncpu, kfold=kfold,
			outdir=outdir, spdir=spdir,
			transl_table=transl_table)
		run_cmd(cmd, log=False)
		
	def get_augustus_training_set(self, records, training_features, fout, flank=1000):
		'''get small DNA'''
		gene_number = 0
		for rc in records:
			record = rc.record
			indices = []
			trans_splicing_records = []
			for feat in rc.get_features('protein'):
				if not feat.pep[0] == 'M':  # start codon
					continue
				if feat in training_features:
					if feat.trans_splicing:     # checked all trans_splicing
						trans_splicing_records += [feat.rearrange_trans_splicing(record)]
						continue
					indices += [feat.index]
			features = [record.features[index] for index in indices]
			for i, feature in enumerate(features):
				start = feature.location.start - flank
				end = feature.location.end + flank
				start = 0 if start < 0 else start
				feat_record = record[start:end]		# sub-record
				feat_record.features = [feat for feat in feat_record.features \
								if (feat.location.start == feature.location.start - start and \
								    feat.location.end == feature.location.end - start)]
				if feature.strand == -1:
					feat_record = feat_record.reverse_complement()
				feat_record.name = feat_record.id = '{}_{}-{}'.format(record.id, start, end)
				if not feat_record.features[0].type == 'source':
					source = SeqFeature(
								FeatureLocation(0, len(feat_record.seq)), 
								type='source', strand=1)
					feat_record.features = [source] + feat_record.features
				SeqIO.write(feat_record, fout, 'genbank')
			for feat_record in trans_splicing_records:
				SeqIO.write(feat_record, fout, 'genbank')
			gene_number += len(indices) + len(trans_splicing_records)
		return gene_number
	def get_seqfile(self, gene):
		return '{}/{}.fa'.format(self.seqdir, gene)
	def get_alnfile(self, gene):
		return '{}/{}.aln'.format(self.msadir, gene)
	def get_hmmfile(self, gene):
		return '{}/{}.hmm'.format(self.hmmdir, gene)
	def get_pflfile(self, gene):
		return '{}/{}.prfl'.format(self.pfldir, gene)

	def get_augustus_species(self, gene):
		# 兼容GeneInfoLine对象 / 字符串两种输入
		if hasattr(gene, "name"):
			gene_name = gene.name
		else:
			gene_name = gene
		# 剥离副本后缀 atp1-a1 -> atp1
		if "-" in gene_name:
			base_gene = gene_name.split("-")[0]
			logger.info(f"Map copy gene {gene_name} to base gene {base_gene}")
		else:
			base_gene = gene_name
		# ogtype为 test_mt_11angio，直接拼接，不再额外加mt前缀
		species = '{}-{}'.format(self.ogtype, base_gene)
		spdir = '{}/species/'.format(AUGUSTUS_CONFIG_PATH,)
		if not os.access(spdir, os.W_OK):
			raise ValueError('{} is not writable'.format(spdir))
		logger.info(f"DEBUG get_augustus_species: species={species}")
		return species

	def align_seqs(self, gene_names, seqfiles, outdir):
		alnfiles = []
		for gene_name, seqfile in zip(gene_names, seqfiles):
			outaln = self.get_alnfile(gene_name)
#			#cmd = 'mafft --auto {} > {}'.format(seqfile, outaln)
			cmd = 'mafft --auto {} | prepareAlign | mafft --auto - > {}'.format(seqfile, outaln)
			run_cmd(cmd)
			alnfiles += [outaln]
		return alnfiles
	def hmmbuild_alns(self, gene_names, alnfiles, outdir):
		hmmfiles = []
		for gene_name, alnfile in zip(gene_names, alnfiles):
			outhmm = self.get_hmmfile(gene_name)
			cmd = 'hmmbuild -n {} {} {}'.format(gene_name, outhmm, alnfile)
			run_cmd(cmd)
			hmmfiles += [outhmm]
		return hmmfiles

	def msa2prfl(self, gene_names, alnfiles, outdir):
		prflfiles = []
		for gene_name, alnfile in zip(gene_names, alnfiles):
			outprfl = self.get_pflfile(gene_name)
			cmd = 'msa2prfl.pl {} --setname={} > {}'.format(alnfile, gene_name, outprfl)
			run_cmd(cmd)
			prflfiles += [outprfl]
		return prflfiles

	def write_orthologs(self, group, d_seqs, outseq):
		with open(outseq, 'w') as fout:
			genes = sorted(group)
			seqs = [d_seqs[gene] for gene in genes]
			for gene, seq in self.remove_abnormal_length(genes, seqs):
				print('>{}\n{}'.format(gene, seq), file=fout)
	def remove_abnormal_length(self, genes, seqs):
		lengths = [len(seq) for seq in seqs]
		q1 = np.percentile(lengths, 25)
		q3 = np.percentile(lengths, 75)
		iqr = q3 - q1
		iqr = max(3, iqr)
		lower = q1 - iqr*8 - 10
		upper = q3 + iqr*8 + 10
		for gene, seq in zip(genes, seqs):
			length = len(seq)
			if length < lower or length > upper:
				print('   remove gene {} with length {} out of [{}, {}]'.format(
						gene, length, lower, upper), file=sys.stderr)
				continue
			yield gene, seq
	def filter_paralogs(self, genes, d_seqs, tmpfix):
		'''for paralogs, cluster for respective ones'''
		group = OrthoMCLGroupRecord(genes=genes)
		single_genes = []
		for sp, genes in list(group.spdict.items()):
#			#gene = max(genes, key=lambda x: len(d_seqs[x]))
			if len(genes) > 1:
				tmpseq = '{}.paralogs'.format(tmpfix)
				clstseq = '{}.paralogs.clust'.format(tmpfix)
				with open(tmpseq, 'w') as f:
					for gene in sorted(genes):
						print('>{}\n{}'.format(gene, d_seqs[gene]), file=f)
				genes = self.cluster_genes(self.cdhit, tmpseq, clstseq)
			single_genes += genes
		return single_genes
	def cluster_genes(self, cdhit, input, output):
		cmd = '{} -i {} -o {} -c 0.95 -aS 0.95 -G 0 -T 4 -d 0 -M 2000'.format(cdhit, input, output)
		run_cmd(cmd)
		return [rc.id for rc in SeqIO.parse(output, 'fasta')]
	def get_gene_name(self, genes, d_names, min_count=3):
		'''get gene name by max count'''
		gene_names = [d_names[gene] for gene in genes]
		return self.get_most_common(gene_names)
	def get_most_common(self, names):
		name_count = Counter(names)
		max_name, max_count = max(list(name_count.items()), key=lambda x: x[1])
		return max_name, name_count

	def get_dbdir(self):
		rootdir = os.path.dirname(os.path.realpath(__file__))
		dbdir = '{}/../db'.format(rootdir)
		return os.path.realpath(dbdir)

	def get_tree(self, taxa, outtree, name_to_taxid=None):
		from ete3 import NCBITaxa
		ncbi = NCBITaxa()

		if name_to_taxid is not None:
			taxid_list = []
			for sp_name in taxa:
				key = sp_name.lower()
				tid = name_to_taxid.get(key)
				if tid is None:
					raise KeyError(f"Cannot find taxid for species: {sp_name}")
				taxid_list.append(tid)
			tree = ncbi.get_topology(taxid_list, intermediate_nodes=True)
		else:
			tree = ncbi.get_topology(list(taxa), intermediate_nodes=True)

		tree.write(outfile=outtree)



	def convert_tree(self, ncbi_tree, species_tree, trefmt='newick'):
		tree = Phylo.read(ncbi_tree, trefmt)
		for clade in tree.get_terminals() + tree.get_nonterminals():
			try:
				name = re.compile(r':name=(.*?) \- \d+:rank').search(clade.comment).groups()[0]
				#name = name.replace(' ', '_')
				name = format_taxon(name)
				clade.name = name
				clade.comment = None
				clade.confidence = None
				clade.branch_length = 0.1
			except TypeError:
				#tree.root_with_outgroup(clade)
				pass
	#	depths = [(clade, depth) for clade, depth in tree.depths().items() if depth > 0]	# clade.clades
	#	clade = min(depths, key=lambda x:x[1])[0]
	#	tree.root_with_outgroup(clade)
	#	print >>sys.stderr, 'rooted by {}'.format(clade)
		tree.root_at_midpoint()
		#tree.rooted = True
		Phylo.write(tree, species_tree, trefmt)
	@lazyproperty
	def templete(self):
		return '{}/{}'.format(self.dbrootdir, 'sqn.template')

class Info(object):
	def __init__(self, infofile, ban=None):
		self.infofile = open(infofile)
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in self.infofile:
			if line.startswith('#'):
				continue
			if hasattr(self, 'ban') and self.ban is not None and line.startswith(self.ban):
				continue
			line = re.compile(r'\t+').split((line.strip()))
			yield self._parse_line(line)
	def _parse_line(self, line):
		return line
	@lazyproperty
	def dict(self):
		return OrderedDict([(line.id, line) for line in self])

def safe_int(x):
    if isinstance(x, int):
        return x
    if not x:
        return 0
    if isinstance(x, str) and x.strip() == "":
        return 0
    return int(x)


def safe_float(x):
    if isinstance(x, (int, float)):
        return float(x)
    if not x:
        return 0.0
    if isinstance(x, str) and x.strip() == "":
        return 0.0
    return float(x)

class GeneInfo(Info):
	def __init__(self, infofile, include_orf=False):
		super(GeneInfo, self).__init__(infofile)
		self.include_orf = include_orf
	def _parse_line(self, line):
		return GeneInfoLine(line)
	@lazyproperty
	def genes(self):
		if self.include_orf:
			for gene in list(self.dict.values()):
				if self.is_orf(gene):
					gene.product = 'hypothetical protein'
			return list(self.dict.values())
		else:
			return [gene for gene in list(self.dict.values()) if not self.is_orf(gene)]
	def is_orf(self, gene):
		return gene.name.lower().startswith('orf') \
			or gene.product == 'hypothetical protein' \
			or re.compile(r'\-gene\d+').search(gene.name)

class GeneInfoLine(object):
	def __init__(self, line=None):
		self.title = ['id', 'name', 'product', 'seq_type',
                	  'exon_count', 'trans_splicing',
                      'gene_count', 'species_count', 
					  'names', 'products', 'species']
		self.ctype = [str, str, str, str, safe_int, bool, safe_int, safe_int, str, str, str]
		self.line = line
		self.set_attr()

	def set_attr(self):
		if self.line is None:
			return
		if isinstance(self.line, str):
			# 从文件读取：tab分隔字符串
			parts = self.line.rstrip("\n").split("\t")
		elif isinstance(self.line, list):
			# 内存直接传入list（makedb内部构造）
			parts = self.line
		else:
			parts = []

		expect = len(self.title)
		if len(parts) == expect:
			# 字段数量完全匹配，正常赋值
			for idx, (key, conv) in enumerate(zip(self.title, self.ctype)):
				val = parts[idx]
				setattr(self, key, conv(val))
		else:
			# 不管是读文件还是内存list，只要长度不对就警告，做容错
			logger.warning(f"[GeneInfoLine] field count mismatch! expect {expect}, got {len(parts)}, raw={repr(self.line)}")
			for idx, (key, conv) in enumerate(zip(self.title, self.ctype)):
				if idx < len(parts):
					val = parts[idx]
				else:
					val = ""
				setattr(self, key, conv(val))


	def __str__(self):
		return self.id
	def write(self, fout):
		if self.line is None:
			print('#' + '\t'.join(self.title), file=fout)
		else:
			print('\t'.join(map(str, self.line)), file=fout)
class SpeciesInfo(Info):
	def __init__(self, infofile):
		super(SpeciesInfo, self).__init__(infofile)
	def _parse_line(self, line):
		return SpeciesInfoLine(line)
	@lazyproperty
	def taxa(self):
		return [line.organism for line in list(self.dict.values())]
	@lazyproperty
	def ntaxa(self):
		return len(set(self.taxa))

class SpeciesInfoLine(GeneInfoLine):
    def __init__(self, line=None):
        self.title = ['id', 'organism', 'name_count', 'cds_count', 'trn_count', 'rrn_count', 'length', 'GC', 'taxonomy']
        self.ctype = [str, str, safe_int, safe_int, safe_int, safe_int, safe_float, safe_float, str]
        self.line = line
        self.set_attr()
		
class NameInfo(Info):
	def __init__(self, infofile):
		super(NameInfo, self).__init__(infofile)
	def _parse_line(self, line):
		return NameInfoLine(line)
	@property
	def dict(self):
		d = {}
		for info in self:
			for name in [info.name] + info.alias:
				d[self.__class__.convert_name(name)] = info
				# tRNA U=T
				if re.compile(r'(trn|tRNA)', re.I).match(name):
					anti_codon = re.compile(r'\-([ATUCG]{3})', re.I).search(name)
					if not anti_codon:
						continue
					anti_codon = anti_codon.groups()[0]
					if 'T' in set(anti_codon):
						new_codon = anti_codon.replace('T', 'U')
					elif 'U' in set(anti_codon):
						new_codon = anti_codon.replace('U', 'T')
					elif 't' in set(anti_codon):
						new_codon = anti_codon.replace('t', 'u')
					elif 'u' in set(anti_codon):
						new_codon = anti_codon.replace('u', 't')
					else:
						continue
					name = name.replace(anti_codon, new_codon)
					d[self.__class__.convert_name(name)] = info
		return d

	@classmethod
	def convert_name(cls, name):
		if re.compile(r'(trn|tRNA)', re.I).match(name):
			lower = r'(^[a-z]+)(.*)([atucg]{3}$)'	 # trnT-ugu	fungi, animal
			upper = r'(^[a-z]+)(.*)([ATUCG]{3}$)'	 # trnT-UGU	plant
			if re.compile(lower).match(name):
				prefix, midfix, suffix = re.compile(lower).match(name).groups()
			elif re.compile(upper).match(name):
				prefix, midfix, suffix = re.compile(upper).match(name).groups()
			try:
				return prefix + midfix.lower() + suffix
			except NameError: pass

		lower = r'(^[a-z]+)(.*)$'
		upper = r'(^[A-Z]+)(.*)$'
		if re.compile(lower).match(name):
			prefix, suffix = re.compile(lower).match(name).groups()
		elif re.compile(upper).match(name):
			prefix, suffix = re.compile(upper).match(name).groups()
		else:
			return name 
		return prefix + suffix.lower()
	@classmethod
	def get_trn_name(cls, raw_name, name_count):
		name_count = [name for name,_ in sorted(list(name_count.items()), key=lambda x:-x[1])]
		if re.compile(r'(trn|tRNA)', re.I).match(raw_name):
			trn = tRNAscanRecord()
			for name in name_count:
				if not re.compile(r'(trn|tRNA)', re.I).match(name):
					continue
				if trn.is_cp(name):	# cp firstly
					try:
						return trn.update_name(name)
					except AttributeError:
						# codon not find
						try:
							for name2 in name_count:
								if not re.compile(r'(trn|tRNA)', re.I).match(name):
									continue
								anti_codon = trn.get_codon(name2)
								if anti_codon is None:
									continue
								return trn.update_name(name, anti_codon=anti_codon)
							else: # codon not find
								return '{}{}-cp'.format(trn.get_prefix(name), trn.get_aa(name))
						except AttributeError: # no aa
							return raw_name
			for name in name_count:
				if not re.compile(r'(trn|tRNA)', re.I).match(name):
					continue
				try:
					return trn.update_name(name)
				except AttributeError:
					# codon not find
					try:
						for name2 in name_count:
							if not re.compile(r'(trn|tRNA)', re.I).match(name2):
								continue
							anti_codon = trn.get_codon(name2)
							if anti_codon is None:
									continue
							return trn.update_name(name, anti_codon=anti_codon)
						else: # codon not find
							return '{}{}'.format(trn.get_prefix(name), trn.get_aa(name))
					except AttributeError: # no aa
						return raw_name
			# upper = r'^([a-z]+)(.*)([ATUCG]{3}).*(\-cp)$'
			# for name in name_count:
				# if re.compile(upper, re.I).match(name):
					# return name
			# for name in name_count:
				# if re.compile(r'(\-cp)', re.I).search(name):
					# return name
			# upper = r'^([a-z]+)(.*)([ATUCG]{3})$'
			# for name in name_count:
				# if re.compile(upper, re.I).match(name):
					# return name
		return raw_name

class NameInfoLine(GeneInfoLine):
	def __init__(self, line=None):
		self.title = ['name', 'product', '_alias']
		self.ctype = [str, str, str]
		self.line = line
		self.set_attr()
		try:
			self.alias = self._alias.split(',')
		except AttributeError:
			self.alias = []

class TaxonomyInfo(Info):
	def __init__(self, infofile):
		super(TaxonomyInfo, self).__init__(infofile)
	def _parse_line(self, line):
		return TaxonomyInfoLine(line)
	@lazyproperty
	def transl_table(self):
		return [line for line in self][0].transl_table
class TaxonomyInfoLine(GeneInfoLine):
	def __init__(self, line=None):
		self.title = ['id', 'db', 'transl_table', 'taxonomy']
		self.ctype = [str, str, int, str]
		self.line = line
		self.set_attr()
		
class Ete3TaxonomyInfo(Info):
	def __init__(self, infolines):
		self.infofile = infolines.split('\n')
		self.ban = 'WARN'
	def _parse_line(self, line):
		return Ete3TaxonomyInfoLine(line)
		
class Ete3TaxonomyInfoLine(GeneInfoLine):
	def __init__(self, line=None):
		self.title = ['taxid', 'sci_name', 'rank', 'named_lineage_', 'taxid_lineage_']
		self.ctype = [int, str, str, str, str]
		self.line = line
		self.set_attr()
	@lazyproperty
	def named_lineage(self):
		return self.named_lineage_.split(',')
	@lazyproperty
	def taxid_lineage(self):
		return list(map(int, self.taxid_lineage_.split(',')))

def main():
	'''example:
-gbfiles *.gb[.gz] -organ mt -taxon rosids				  # use records of rosids and name db as mt-rosids
-gbfiles *.gb[.gz] -organ mt -taxon rosids -custom myself # use records of rosids but name db as mt-myself
-gbfiles *.gb[.gz] -organ mt -custom myself               # use all records in gbfiles and name db as mt-myself
'''
	args = makeArgparse()
	print(args.__dict__, file=sys.stderr)
	db = Database(**args.__dict__)
	if args.check:
		db.checkdb(untar=False)
	elif args.list:
		db.listdb()
	elif args.count:
		db.count_taxa()
	else:
		db.makedb()

def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-list', action='store_true', default=False,
						help="list database")
	parser.add_argument('-check', action='store_true', default=False,
						help="check database")
	parser.add_argument('-count', action='store_true', default=False,
						help="count taxa")
	parser.add_argument('-gbfiles', type=str, nargs='+', default=None, 
						help="genbank files")
	parser.add_argument('-organ', type=str, choices=['mt', 'pt'], default='mt',
						help="mt (mitochondrion) or pt (plastid) to build [default=%(default)s]")
	parser.add_argument('-dbdir', type=str, default=None,
						help="directory to output [default=auto]")
	parser.add_argument('-taxon', type=str, default=None,
						help="taxon to select, such as Embryophyta [default=%(default)s]")
	parser.add_argument('-custom', type=str, default=None,
						help="use custom mark instead taxon [default=%(default)s]")
	parser.add_argument('-exclude_taxa', type=str, default=None,
						help="exclude taxa seperated by comma, for example: 'Himantura,Vazella' [default=%(default)s]")
	parser.add_argument('-tmpdir', type=str, default='/dev/shm/tmp',
						help="tmpdir")
	parser.add_argument('-cleanup', action='store_true', default=False,
						help="clean up")
	
	args = parser.parse_args()
	if args.exclude_taxa is not None:
		args.exclude_taxa = set(args.exclude_taxa.split(','))
	return args

if __name__ == '__main__':
	main()
