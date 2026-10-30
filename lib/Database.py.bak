import sys, os
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
from lazy_property import LazyWritableProperty as lazyproperty

from Genbank import GenbankParser, format_taxon
from OrthoFinder import OrthoFinder, OrthoMCLGroupRecord
from tRNA import tRNAscanRecord
from RunCmdsMP import logger, run_cmd, run_job
from small_tools import mkdirs, rmdirs

AUGUSTUS_CONFIG_PATH=os.environ['AUGUSTUS_CONFIG_PATH']

d_lineage = {}
class Database():
	def __init__(self, organ=None, taxon=None, dbrootdir=None,
				gbfiles=None, custom=None, version=None,
				include_orf=False,
				tmpdir='/dev/shm/tmp', 
				ncpu=24,
				upper_limit=500,
				lower_limit=2,
				exclude_taxa=None,
				cleanup=True,
				**kargs):
		if dbrootdir is None:
			dbrootdir = self.get_dbdir()
		self.taxon = taxon
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
		self.vesion_file = '{}/version'.format(self.dbdir)
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
			print >> sys.stdout, db.organ, db.taxon
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
			print >> sys.stderr, 'Using existed', organism
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
		self.version = open(self.vesion_file).read().strip()
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
		for taxon, bin_records in d_records.items():
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
			n_per_genus = int(round(1.0*self.upper_limit / ngenera, 0))
			logger.info('pruning by {} with > {} organisms'.format(rank, n_per_genus))
			return self.prune_by_rank(records, rank, n_per_genus)
		else:
			records = self.prune_by_rank(records, rank, 1)
		
		# prune by family evenly
		rank = 'family'
		families = {record.family for record in records if record.family is not None}
		nfamily = len(families)
		if nfamily <= self.upper_limit:
			n_per_family = int(round(1.0*self.upper_limit / nfamily, 0))
			logger.info('pruning by {} with > {} organisms'.format(rank, n_per_family))
			d_records = {}  # by family
			return self.prune_by_rank(records, rank, n_per_family)
		else:
			records = self.prune_by_rank(records, rank, 1)
		
		# prune by order evenly, regardless of upper limit
		rank = 'order'
		orders = {record.order for record in records if record.order is not None}
		norder = len(orders)
		n_per_order = int(round(1.0*self.upper_limit / norder, 0))
		if n_per_order < 1:
			n_per_order = 1
		logger.info('pruning by {} with > {} organisms'.format(rank, n_per_order))
		return self.prune_by_rank(records, rank, n_per_order)
	
	def get_species_tree(self, taxa):
		tmptree = '{}/species.tree.tmp'.format(self.tmpdir)
		self.get_tree(taxa, tmptree)
		species_tree = '{}/species.treefile'.format(self.dbdir)
		self.convert_tree(tmptree, species_tree)
		return species_tree
		
	def count_taxa(self):
		gb = GenbankParser(self.gbfiles)
		records = [record.prune() for record in gb.filter_by_taxon(taxon=self.taxon)]
		self.count_rank(records)
	
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
			print >>sys.stderr, 'check point: ' + _ckp
			if not os.path.exists(_ckp):
				_cmd = cmd.format(dirname=dirname, basename=basename)
				run_cmd(_cmd, log=True)
	def makedb(self):
		rmdirs(self.pepdir, self.rnadir)  # for OrthoFinder
		mkdirs(self.pepdir, self.rnadir, self.alndir, self.msadir, self.trndir) 
		mkdirs(self.seqdir, self.hmmdir, self.pfldir, self.tnddir)
		logger.info('parsing {}'.format(self.gbfiles))
		gb = GenbankParser(self.gbfiles)
		if self.version is None:
			version = time.strftime("%Y-%m-%d", time.gmtime(os.path.getmtime(self.gbfiles[0])))
			today = time.strftime("%Y-%m-%d", time.gmtime(time.time()))
			if re.compile(r'\S+.genomic.gbff').search(self.gbfiles[0]):
				source = 'RefSeq'
			else:
				source = 'Unknown'
			self.version = '{} {}: build {}'.format(source, version, today)
			logger.info('version: {}'.format(self.version))
			with open(self.vesion_file, 'w') as fout:
				print >>fout, self.version
		ofbin = 'orthofinder.py'
		# read into MEM
		records = list(gb.filter_by_taxon(taxon=self.taxon, 
							exclude=self.exclude_taxa,
							taxonomy_dbfile=self.taxonomy_dbfile))
		# count
		self.count_rank(records)
		# prune
		logger.info('pruning taxa')
		records = self.prune(records)
		self.count_rank(records)

		taxa = [record.organism for record in records]
		logger.info('building taxonomy tree')
		sp_tree = self.get_species_tree(set(taxa))
		logger.info('using speceis tree `{}`'.format(sp_tree))
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
			
			features = gb.write_seqs_by_taxon(seqdir, feat_type=feat_type, records=records, taxon=self.taxon)
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
				taxonomy = [';'.join(tax) for tax in gb.tax]
				taxon, taxon_count = self.get_most_common(taxonomy)
				logger.info('get taxon `{}` from {}'.format(taxon, taxon_count))
				tables = [feat.transl_table for feat in features]
				transl_table, table_count = self.get_most_common(tables)
				logger.info('get transl_table `{}` from {}'.format(transl_table, table_count))
				line = [self.organ, self.ogtype, transl_table, taxon]
				with open(self.taxonomy_file, 'w') as fout:
					TaxonomyInfoLine().write(fout)
					TaxonomyInfoLine(line).write(fout)

			logger.info('running OrthoFinder to cluster genes')
			if sp_tree is not None:
				of_opts += ' -s {}'.format(sp_tree)
		#	of_opts += ' -og'
			cmd = '{} -f {} -t {} {}'.format(ofbin, seqdir, self.ncpu, of_opts)
			run_cmd(cmd, log=True)

			logger.info('parsing OrthoFinder orthologs')
			orfdir = glob.glob('{}/OrthoFinder/*/'.format(seqdir))[0]
			orf = OrthoFinder(orfdir)
			logger.info('parsing orthologs group into {}'.format(self.alndir))
			min_count = ntaxa * min_ratio
			self.min_count = max(self.lower_limit, min_count)
			info_file = prefix + '.raw'

			# parse orthologs, bin sequences by orthologs
			with open(info_file, 'w') as fout:
				gene_names, seqfiles, groups = \
						self.parse_orthologs(orf, features, self.alndir, fout) # names
			self.check_unique(seqfiles)
			#if not os.path.exists(prefix):
			os.rename(info_file, prefix)
			logger.info('extrat {} genes: {}'.format(
					len(gene_names), sorted(gene_names)))	# gene id

			logger.info('aligning sequences')
			self.alnfiles = self.align_seqs(gene_names, seqfiles, outdir=self.msadir)
			self.check_unique(self.alnfiles)
			logger.info('building hmm profiles')
			self.hmmfiles = self.hmmbuild_alns(gene_names, self.alnfiles, outdir=self.hmmdir)
			self.check_unique(self.hmmfiles)
			if build_pfl:
				logger.info('building augustus profiles')
				self.pflfiles = self.msa2prfl(gene_names, self.alnfiles, outdir=self.pfldir)
				logger.info('training augustus models')
				# train augustus
				d_features = {feat.id: feat for feat in features}
				for gene_name, group in zip(gene_names, groups):
					training_features = {d_features[gene] for gene in group}
					train_db = '{}/train-{}.gb'.format(self.trndir, gene_name)
					with open(train_db, 'w') as fout:
						gene_number = self.get_augustus_training_set(records, training_features, fout)
					sepcies = '{}-{}'.format(self.ogtype, gene_name)
					self.train_augustus(train_set=train_db, species=sepcies, outdir=self.tnddir,
						gene_number=gene_number, ncpu=self.ncpu, transl_table=transl_table)
		# tar
		self.tgz_dirs(self.seqdir, self.msadir, self.hmmdir, self.pfldir, self.tnddir)
		if self.cleanup:
			logger.info('cleaning')
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
		gene_number = gene_number/3
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
		species = '{}-{}'.format(self.ogtype, gene)
		spdir = '{}/species/'.format(AUGUSTUS_CONFIG_PATH,)
		if not os.access(spdir, os.W_OK):
			raise ValueError('{} is not writable'.format(spdir))
		species_dir = '{}/{}'.format(spdir, species)
		species_tgz = '{}/{}.tgz'.format(self.tnddir, species)
		if os.path.exists(species_dir):
			pass
		elif os.path.exists(species_tgz):
			cmd = 'cd {spdir} && tar xzf {tgz}'.format(spdir=spdir, tgz=species_tgz)
			run_cmd(cmd, log=False)
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
				print >>fout, '>{}\n{}'.format(gene, seq)
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
				print >> sys.stderr, '   remove gene {} with length {} out of [{}, {}]'.format(
						gene, length, lower, upper)
				continue
			yield gene, seq
	def filter_paralogs(self, genes, d_seqs, tmpfix):
		'''for paralogs, cluster for respective ones'''
		group = OrthoMCLGroupRecord(genes=genes)
		single_genes = []
		for sp, genes in group.spdict.items():
#			#gene = max(genes, key=lambda x: len(d_seqs[x]))
			if len(genes) > 1:
				tmpseq = '{}.paralogs'.format(tmpfix)
				clstseq = '{}.paralogs.clust'.format(tmpfix)
				with open(tmpseq, 'w') as f:
					for gene in sorted(genes):
						print >>f, '>{}\n{}'.format(gene, d_seqs[gene])
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
		max_name, max_count = max(name_count.items(), key=lambda x: x[1])
		return max_name, name_count

	def get_dbdir(self):
		rootdir = os.path.dirname(os.path.realpath(__file__))
		dbdir = '{}/../db'.format(rootdir)
		return os.path.realpath(dbdir)

	def get_tree(self, taxa, outre):
		taxa = ' '.join(map(repr, taxa))
		cmd = 'ete3 ncbiquery --search {} --tree >  {}'.format(taxa, outre)
		run_cmd(cmd, log=True)

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

class GeneInfo(Info):
	def __init__(self, infofile, include_orf=False):
		super(GeneInfo, self).__init__(infofile)
		self.include_orf = include_orf
	def _parse_line(self, line):
		return GeneInfoLine(line)
	@lazyproperty
	def genes(self):
		if self.include_orf:
			for gene in self.dict.values():
				if self.is_orf(gene):
					gene.product = 'hypothetical protein'
			return self.dict.values()
		else:
			return [gene for gene in self.dict.values() if not self.is_orf(gene)]
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
		self.ctype = [str, str, str, str, int, bool, int, int, str, str, str]
		self.line = line
		self.set_attr()
	def set_attr(self):
		if self.line is not None:
			for key, value, type in zip(self.title, self.line, self.ctype):
				setattr(self, key, type(value))
	def __str__(self):
		return self.id
	def write(self, fout):
		if self.line is None:
			print >>fout, '#' + '\t'.join(self.title)
		else:
			print >>fout, '\t'.join(map(str, self.line))
class SpeciesInfo(Info):
	def __init__(self, infofile):
		super(SpeciesInfo, self).__init__(infofile)
	def _parse_line(self, line):
		return SpeciesInfoLine(line)
	@lazyproperty
	def taxa(self):
		return [line.organism for line in self.dict.values()]
	@lazyproperty
	def ntaxa(self):
		return len(set(self.taxa))

class SpeciesInfoLine(GeneInfoLine):
	def __init__(self, line=None):
		self.title = ['id', 'organism', 'name_count', 'cds_count', 'trn_count', 'rrn_count', 'length', 'GC', 'taxonomy']
		self.ctype = [str, str, str]
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
		name_count = [name for name,_ in sorted(name_count.items(), key=lambda x:-x[1])]
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
		return map(int, self.taxid_lineage_.split(','))

def main():
	'''example:
-gbfiles *.gb[.gz] -organ mt -taxon rosids				  # use records of rosids and name db as mt-rosids
-gbfiles *.gb[.gz] -organ mt -taxon rosids -custom myself # use records of rosids but name db as mt-myself
-gbfiles *.gb[.gz] -organ mt -custom myself               # use all records in gbfiles and name db as mt-myself
'''
	args = makeArgparse()
	print >>sys.stderr, args.__dict__
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
