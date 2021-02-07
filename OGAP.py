import sys, os
import copy
import argparse
import uuid
from collections import OrderedDict
from itertools import izip, combinations
from Bio import SeqIO

from lib.Database import Database
from lib.Hmmer import HmmSearch
from lib.tRNA import tRNAscan, tRNAscanStructs
from lib.Gff import ExonerateGffGenes, AugustusGtfGenes
from lib.Psl import PslParser
from lib.Repeat import RepeatPipeline, addRepeatArgs
from lib.RunCmdsMP import run_cmd, run_job, logger
from lib.translate_seq import six_frame_translate
from lib.small_tools import mkdirs, rmdirs
from lib.small_tools import open_file as open

bindir = os.path.dirname(os.path.realpath(__file__))
os.environ['PATH'] = bindir + ':' + os.environ['PATH']

LOCATION = {'pt': 'chloroplast', 'mt': 'mitochondrion'}

def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("genome", action="store",type=str,
					help="input genome sequence in fasta or genbank format [required]")
	parser.add_argument('-organ', type=str, choices=['mt', 'pt'], default=None,
					help="organelle type: mt (mitochondrion) or pt (plastid) [default=%(default)s]")
	parser.add_argument('-mt', action="store_true", default=False,
					help="equal to '-organ mt' [default=%(default)s]")
	parser.add_argument('-pt', action="store_true", default=False,
					help="equal to '-organ pt' [default=%(default)s]")
	parser.add_argument('-taxon', type=str, default=None, nargs='+',
					help="taxon to use as reference, such as rosids [default=auto by -organism]")
	parser.add_argument('-trans', action="store_true", default=False,
					help='transcipt input mode [default=%(default)s]')
	parser.add_argument("-tmpdir", action="store",
					default='/dev/shm/tmp', type=str,
					help="temporary directory [default=%(default)s]")
	parser.add_argument('-cleanup', action="store_true", default=False,
					help='clean up temporary directory [default=%(default)s]')
					
	parser.add_argument('-seqfmt', type=str, choices=['fasta', 'genbank'], default=None,
					help="genome seqence format [default=auto]")
	parser.add_argument("-est", action="store",type=str,
					help="EST sequences as evidence in fasta format for coding genes annotation")

	parser.add_argument('-sp', '-organism', type=str, default=None, dest='organism',
					help="organism to be included in sqn, required for fasta input [default=%(default)s]")
	parser.add_argument('-linear', action="store_true", default=False,
					help="topology to be included in sqn [default=circular]")
	parser.add_argument('-partial', action="store_true", default=False,
					help="completeness to be included in sqn [default=complete]")
	parser.add_argument('-no_sqn', action="store_true", default=False,
                    help="do not generate sqn file [default=%(default)s]")
	parser.add_argument('-wgs', action="store_true", default=False,
                    help="for WGS submission [default=%(default)s]")
	parser.add_argument('-sqn_annot', action="store_true", default=False,
                    help="only include sequences with annotaion into sqn file [default=%(default)s]")
	parser.add_argument('-no_cds', action="store_true", default=False,
                    help="do not annotate coding gene [default=%(default)s]")
	parser.add_argument('-no_rrn', action="store_true", default=False,
                    help="do not annotate rRNA [default=%(default)s]")
	parser.add_argument('-no_trn', action="store_true", default=False,
                    help="do not annotate tRNA [default=%(default)s]")

	group_out = parser.add_argument_group('output',)
	group_out.add_argument('-o', "-outdir", action="store", dest='outdir',
					default='.', type=str,
					help="output directory [default=%(default)s]")
	group_out.add_argument('-pre', "-prefix", action="store", dest='prefix',
					default=None, type=str,
					help="output prefix [default=genome file basename]")
	parser.add_argument('-min_cds_cov', type=float, default=40,
					help="min coverage to filter candidate coding genes [default=%(default)s]")
	parser.add_argument('-min_rrn_cov', type=float, default=40,
					help="min coverage to filter candidate rRNA genes [default=%(default)s]")
	parser.add_argument('-min_trn_cov', type=float, default=40,
					help="min coverage to filter candidate tRNA genes [default=%(default)s]")
	group_out.add_argument('-min_cov', type=float, default=25,
					 help="min coverage to filter out final gene set [default=%(default)s]")
	group_out.add_argument('-min_score', type=float, default=15,
					 help="min score to filter out final gene set [default=%(default)s]")
#	parser.add_argument('-score_cutoff', type=float, default=0.85,
#					help="min score ratio of the highest of multi-copy gene to output [default=%(default)s]")
	group_out.add_argument('-min_cds_score', type=float, default=0.85,
					help="min score ratio to filter duplicated coding genes [default=%(default)s]")
	group_out.add_argument('-min_rrn_score', type=float, default=0.95,
					help="min score ratio to filter duplicated rRNA genes [default=%(default)s]")
	group_out.add_argument('-min_trn_score', type=float, default=0.8,
					help="min score ratio to filter duplicated tRNA genes [default=%(default)s]")
	group_out.add_argument('-orf', '-include_orf', action="store_true", default=False,
					dest='include_orf',
					help="include ORF genes in coding gene annotation [default=%(default)s]")
	group_out.add_argument('-trn_struct', action="store_true", default=False,
					help="output tRNA structure [default=%(default)s]")
	group_out.add_argument('-draw_map', action="store_true", default=False,
					help="draw gene map [default=%(default)s]")
	group_out.add_argument('-compare_map', action="store_true", default=True,
                    help="compare gene map with genbank input [default=%(default)s]")

	group_rep = parser.add_argument_group('repeat', )
	group_rep.add_argument('-repeat', action="store_true", default=False,
					help="output repeats [default=%(default)s]")
	addRepeatArgs(group_rep)

	args = parser.parse_args()
	if args.organism is not None:
		args.organism = args.organism.replace('_', ' ')
	if args.organ is None:
		if args.mt:
			args.organ = 'mt'
		elif args.pt:
			args.organ = 'pt'
		else:
			raise ValueError('no organ type specified')
	return args

class Pipeline():
	def __init__(self, genome, 
				organ, taxon, 
				trans=False,
				outdir='.',
				prefix=None,
				transl_table=None,
				est=None, # EST evidence
				tmpdir='/dev/shm/tmp', 
				organism=None,
				linear=False, 
				partial=False,
				nosqn=False, wgs=False,
				no_cds=False,		# do not annotate CDS
				no_rrn=False,
				no_trn=False,
				sqn_annot=False,	# only sequences with annotation to sqn
				include_orf=False,	# annotate ORF
				exon_diff_penalty = 100,
				min_cds_hmmcov=0,	# min coverage of each hmm hit
				min_rrn_hmmcov=5,
				min_trn_hmmcov=5,
				min_cds_cov=40,		# min coverage of sum of hmm hits
				min_rrn_cov=40,		# 
				min_trn_cov=40,		# 
				cov_cutoff=0.75,	# filter candidate (> highest * cov_cutoff)
				#score_cutoff=0.85,	# filter final set (> highest * score_cutoff) to filter out multi-copy
				min_cds_score=0.85,		# score_cutoff for cds
				min_rrn_score=0.95,		# score_cutoff for rrn
				min_trn_score=0.80,		# score_cutoff for trn
				min_score=20, # filter final set (> hard_score_cutoff) to filter out single copy
				min_cov=50,			# filter final set (> min_cov%) to filter out single copy
				trn_opts=' -O',
				trn_struct=False,
				draw_map=True,
				compare_map=True,
				repeat=False,
				cleanup=True,
				seqfmt='fasta', **kargs):
		
		self.organ = organ
		self.genome = os.path.realpath(genome)
		self.organ = organ
		self.taxon = taxon
		self.trans = trans
		self.outdir = os.path.realpath(outdir)
		self.est = est
		self.organism = organism
		self.linear = linear
		self.partial = partial
		self.nosqn = nosqn
		self.wgs = wgs
		self.sqn_annot = sqn_annot
		self.no_cds = no_cds
		self.no_rrn = no_rrn
		self.no_trn = no_trn
		self.include_orf = include_orf
		self.trn_opts = trn_opts
		self.trn_struct = trn_struct
		self.draw_map = draw_map
		self.compare_map = compare_map
		self.repeat =repeat
		self.cleanup = cleanup
		self.exon_diff_penalty = exon_diff_penalty
		self.kargs = kargs
		if seqfmt is None:
			self.seqfmt = self.guess_seqfmt(self.genome)
			logger.info('genome in {} format'.format(self.seqfmt))
		else:
			self.seqfmt = seqfmt
		if self.seqfmt == 'fasta' and self.organism is None:
			raise ValueError('-sp is required for fasta-format input')
		if self.seqfmt == 'genbank':
			rc = self.get_record(self.genome)
			if self.organism is None:
				try:
					self.organism = rc.annotations['organism'].replace(' Unclassified.', '')
					logger.info('using organism: `{}`'.format(self.organism))
				except KeyError: pass
			if not self.linear:
				try:
					if rc.annotations['topology'] == 'linear':
						self.linear = True
				except KeyError: pass
			if not self.partial:
				if rc.description.find('complete') < 0:
					self.partial = True
		# min coverage of one hmm hit
		self.min_cds_hmmcov = min_cds_hmmcov
		self.min_rrn_hmmcov = min_rrn_hmmcov
		self.min_trn_hmmcov = min_trn_hmmcov
		if self.trans:	# expect complete gene region
			self.min_cds_hmmcov = 90
			self.min_rrn_hmmcov = 90
			self.min_trn_hmmcov = 90
		# min coverage
		self.min_cds_cov = min_cds_cov
		self.min_rrn_cov = min_rrn_cov
		self.min_trn_cov = min_trn_cov
		
		self.cov_cutoff = cov_cutoff
		#self.score_cutoff = score_cutoff
		self.score_cutoff = {'mRNA':min_cds_score, 'rRNA':min_rrn_score, 'tRNA':min_trn_score}
		self.min_score = min_score
		self.min_cov = min_cov

		if taxon is None and self.organism is not None:
			taxon = Database(organ=organ).select_db(self.organism).taxon
			taxon = [taxon]
			logger.info('automatically select db: {}'.format(taxon))
		elif taxon is None:
			logger.info('neither -taxon or -organism must be specified')
		self.taxon = taxon	# taxa
		#self.db = Database(organ=organ, taxon=taxon, include_orf=include_orf)
		
		#self.ogtype = self.db.ogtype
		self.transl_table = transl_table
		
		if prefix is None:
			self.prefix = os.path.basename(genome)
		else:
			self.prefix = prefix
		# folder
		self.tmpdir0 = tmpdir
		#self.tmpdir = '{}/{}/{}'.format(tmpdir, self.ogtype, 
		#						os.path.basename(self.prefix))
		#self.hmmoutdir = '{}/hmmout'.format(self.tmpdir)
		##self.estoutdir = '{}/estout'.format(self.tmpdir)
		##self.exnoutdir = '{}/exnout'.format(self.tmpdir)
		#self.agtoutdir = '{}/augustus'.format(self.tmpdir)
		#self.trndir = '{}/{}.trna'.format(self.outdir, self.prefix)

	def run(self):
		# draw map
#		self.seqs = self.get_seqs(self.genome, self.seqfmt)
#		self.seqlen = sum([len(seq) for seq in self.seqs.values()])
#		if self.drawgenemap:
#			self.draw_map()
#			if self.seqfmt == 'genbank':
#				self.compare_map()
#		return
		mkdirs(self.outdir)
#		rmdirs(self.agtoutdir)
#		mkdirs(self.outdir, self.tmpdir)
#		mkdirs(self.hmmoutdir, self.agtoutdir)
		# check
#		logger.info('checking database: {}'.format(self.db.ogtype))
#		self.db.checkdb()
#		if self.transl_table is None:
#			self.transl_table = self.db.transl_table
	
		# read genome seqs
		self.seqs = self.get_seqs(open(self.genome), self.seqfmt)
		self.seqlen = sum([len(seq) for seq in self.seqs.values()])
		seqids = self.seqs.keys()
		self.nseqs = len(seqids)
		if self.nseqs > 1:
			logger.info('changing partial to True due to nseqs>1')
			self.linear = self.partial = True
			self.sqn_annot = True
		if self.contains_gap(self.seqs.values()):
			logger.info('changing partial to True due to non-ATCG gap(s)')
			self.partial = True

		# to fasta
		self.fsa = self.to_fsa()

		records = []
		for taxon in self.taxon:
			db_records = []
			self.db = Database(organ=self.organ, taxon=taxon, include_orf=self.include_orf)
			self.ogtype = self.db.ogtype
			# folder
			uid = uuid.uuid1()
			self.tmpdir = '{}/{}/{}-{}'.format(self.tmpdir0, self.ogtype,
                               os.path.basename(self.prefix), uid)
			self.hmmoutdir = '{}/hmmout'.format(self.tmpdir)
			self.agtoutdir = '{}/augustus'.format(self.tmpdir)
			self.trndir = '{}/{}.trna'.format(self.outdir, self.prefix)
			rmdirs(self.agtoutdir)
			mkdirs(self.tmpdir)
			mkdirs(self.hmmoutdir, self.agtoutdir)
			# check
			logger.info('checking database: {}'.format(self.db.ogtype))
			self.db.checkdb()
			if self.transl_table is None:
				self.transl_table = self.db.transl_table
			# cds-protein finder
			if self.est is not None:
				self.hmmsearch_est()
			logger.info('finding protein-coding genes by HMM+enonerate+augustus')
			db_records += self.hmmsearch_protein()
		
			# rna finder
			logger.info('finding non-coding genes by HMM+exonerate (rRNA) or HMM+tRNAscan-SE(tRNA)')
			db_records += self.hmmsearch_rna()
			for record in db_records:
				record[0].attributes['db'] = taxon
			records += db_records
		# score
		#logger.info('scoring genes by HMM')
		#records = [self.score_record(record) for record in records]
		# remove duplicates
		logger.info('removing duplicated genes with the identical coordinate')
		records = self.remove_duplicates(records)
		# remove low quality of the same gene
		logger.info('removing duplicated genes with lower score: {} * highest score'.format(self.score_cutoff))
		records = self.remove_lowqual(records, cutoff=self.score_cutoff, 
					hard_cutoff=self.min_score, min_cov=self.min_cov)
		
		# out fasta of gene, sorted by score
		self.to_fasta(records)
		# tRNA structure
		if self.trn_struct:
			self.plot_struct(self.trndir, records)
		
		# repeat	
		if self.repeat:
			rep = RepeatPipeline(genome=self.fsa, tmpdir=self.tmpdir, prefix=self.prefix, **self.kargs)
			records += rep.run()
		# sort by coordinate
		records = sorted(records, key=lambda x: (seqids.index(x.chrom), x.start))
		
		# to gff
		self.to_gff3(records)

		# out sqn
		if not self.nosqn:
			if self.sqn_annot:
				self.re_fsa(records)
			logger.info('outputing sqn for submitting Genbank')
			self.to_sqn(records)

			# draw map
			if self.draw_map:
				self.draw_gene_map()
			if self.compare_map and self.seqfmt == 'genbank':
				self.compare_gene_map()
	
		# summary
		self.summary_source(records)
		self.summary_records(records)
		
		# clean up
		if self.cleanup:
			logger.info('cleaning')
			rmdirs(self.tmpdir)

	def remove_lowqual(self, records, cutoff={}, hard_cutoff=20, min_cov=50):
		'''remove duplicates with same name and low qual'''
		d_group = OrderedDict()
		for record in records:
			key = record.name		# id or name?
			try: d_group[key] += [record]	# duplicates from homologous gene
			except KeyError: d_group[key] = [record]
		better_records = []
		for key, records in d_group.items():
			count = len(records)
			if count == 1:	# unique
				better_records += records
				continue
			
			# remove that with low score
			highest_record = max(records, key=lambda x:x.score)
			highest_score = highest_record.score
			good_cutoff = highest_score * cutoff[highest_record.rna_type]
			better_record = [record for record in records if record.score >= good_cutoff]
			#print >>sys.stderr, key, count, '->', len(better_record)
			# remove that high qual but with more parts
			top_record = [record for record in records if record.score >= highest_score*0.96]
			min_npart = min([record.npart for record in top_record])
			better_record = [record for record in better_record if record.npart <= min_npart]
			
			# remove overlap
			if len(better_record) >1:
				better_record = self.remove_overlaps(better_record)
			
			print >>sys.stderr, key, count, '->', len(better_record)
			better_records += better_record
		
		
		filtered_records = []
		for record in better_records:
			if record.score < hard_cutoff:
				print >>sys.stderr, record, record.name, 'removed with too low score: {}'.format(record.score)
				continue
			if record.cov < min_cov:
				print >>sys.stderr, record, record.name, 'removed with too low coverage: {}'.format(record.cov)
				continue
			filtered_records += [record]
		return filtered_records
	
	def remove_duplicates(self, records):
		'''remove duplicates with same coordinate'''
		keys = set([])
		d_group = OrderedDict()
		for record in records:
			key1 = (str(record),)	# coordinate
			key2 = key1 + (record.id, )
			if key2 in keys:	# duplicates from the same gene
				print >>sys.stderr, key2, 'removed'
				continue
			try: d_group[key1] += [record]	# duplicates from homologous gene
			except KeyError: d_group[key1] = [record]
			keys.add(key2)
		unique_records = []
		for key, records in d_group.items():	# duplicates from different genes with same coordinate
			if len(records) == 1:	# unique
				unique_records += records
				continue
			best_record = max(records, key=lambda x:x.score)
			print >>sys.stderr, key, best_record.name, len(records), '->', 1
			unique_records += [best_record]
		return unique_records
		
	def remove_overlaps(self, records):
		while True:
			records = sorted(records, key=lambda x: (x.start, -x.end))
			overalped_records = []
			for rc1, rc2 in combinations(records, 2):
				#print >>sys.stderr, 'remove_overlaps'
				if rc1.overlaps(rc2):
					overalped_record = min([rc1, rc2], key=lambda x:x.score)
					print >>sys.stderr, overalped_record, overalped_record.name, 'removed with overlap'
					overalped_records += [overalped_record]
			if len(overalped_records) == 0:
				break
			records = set(records) - set(overalped_records)
		return records
	def score_record(self, record):
		rnafa = self.get_filename(self.hmmoutdir, record.gene_id, 'fasta')
		with open(rnafa, 'w') as fout:
			try:
				print >> fout, '>{} {}\n{}'.format(record.gene_id, record, record.pep_seq)
			except:
				print >> fout, '>{} {}\n{}'.format(record.gene_id, record, record.rna_seq)
		hmmfile = self.db.get_hmmfile(record.id)
		domtblout = rnafa + '.domtbl'
		self.hmmsearch(hmmfile, rnafa, domtblout)
		hmm_best = HmmSearch(domtblout).get_best_hit()
		try:
			record.score = round(hmm_best.edit_score, 1)
			record.cov = round(hmm_best.cov, 1)
		except AttributeError:
			record.score = 0
			record.cov = 0
		record[0].score = record.score
		record[0].attributes.update(cov=record.cov)
		return record
		
	def summary_records(self, records):
		d_smy = OrderedDict()
		for record in records:
			try: d_smy[record.rna_type] += [record]
			except KeyError: d_smy[record.rna_type] = [record]
		logger.info('summary by gene type:')
		line = ['type', 'copy number', 'gene number', 'gene names']
		print >>sys.stdout, '\t'.join(line)
		self.print_summary(d_smy)
	def print_summary(self, d_smy):
		for rna_type, records in d_smy.items():
			genes = [record.id for record in records]
			names = list({record.name for record in records})
			names = sorted(names)
			line = [rna_type, len(genes), len(names), ','.join(names)]
			line = map(str, line)
			print >>sys.stdout, '\t'.join(line)
	def summary_source(self, records):
		d_smy = OrderedDict()
		for record in records:
			try: d_smy[record.source] += [record]
			except KeyError: d_smy[record.source] = [record]
		logger.info('summary by source:')
		line = ['source', 'copy number', 'gene number', 'gene names']
		print >>sys.stdout, '\t'.join(line)
		self.print_summary(d_smy)

	def hmmsearch_rna(self):
		records = []
		if self.no_rrn and self.no_trn:
			return records
		na_seq = '{}/{}.genome.na'.format(self.tmpdir, self.prefix)
		na_seq = '{}/{}.genome.na'.format(self.tmpdir, self.prefix)
		with open(na_seq, 'w') as fout:
			d_length = self.double_seqs(self.fsa, fout, seqfmt='fasta') #self.seqfmt)
		#trn_records = []
		for gene in self.db.rna_genes:
			if gene.seq_type == 'rRNA' and self.no_rrn:
				continue
			if gene.seq_type == 'tRNA' and self.no_trn:
				continue
			print >>sys.stderr, '\n   >> {}: {}'.format(gene, gene.name)
			hmmfile = self.db.get_hmmfile(gene)
			domtblout = self.get_domtblout(gene, src='g')
			self.hmmsearch(hmmfile, na_seq, domtblout)
			if gene.seq_type == 'rRNA':
				structs = HmmSearch(domtblout).get_gene_structure(d_length, 
								min_ratio=self.cov_cutoff, 
								min_hmmcov=self.min_rrn_hmmcov, min_part=2,
								min_cov=self.min_rrn_cov, seq_type='nucl', flank=2000)
				for i, parts in enumerate(structs):
					parts.id = '{}-{}'.format(gene, i+1)
					parts = parts.link_part()
					print >> sys.stderr, 'old', parts.to_str()
					genefa = self.get_filename(self.agtoutdir, parts, 'fa')
					with open(genefa, 'w') as fout:
						parts.write_seq(self.seqs, fout)
					#self.exonerate_gene_est(genefa, gene, parts)
					best_exons = self.exonerate_est2genome(genefa, gene, parts, minintron=500)
					if best_exons is None:
						continue
					rrn, new_parts = parts.map_coord(best_exons)	# GffExons
					new_parts = new_parts.link_part()
					print >> sys.stderr, 'new', new_parts.to_str()
					#rrn = rrn.link_exons(minintron=200)
					new_parts.source = 'blat'
					record = rrn.extend_gene(gene, new_parts, rna_type=gene.seq_type)
					# new_parts = parts
					# parts.source = 'hmmsearch'
					# exons = parts.to_exons()
					# #exons.write(sys.stderr)
					# record = exons.extend_gene(gene, parts, rna_type=gene.seq_type)
					print >> sys.stderr, ''
					rna_seq = record.extract_seq(self.seqs)
					record.rna_seq = rna_seq
					record = self.score_record(record)
					record.npart = len(new_parts)
					record.write(sys.stderr)
					records += [record]
			elif gene.seq_type == 'tRNA':
				#continue
				structs = HmmSearch(domtblout).get_gene_structure(d_length,
								min_ratio=self.cov_cutoff, 
								min_hmmcov=self.min_trn_hmmcov, min_part=2, # diff
							    min_cov=self.min_trn_cov, seq_type='nucl', flank=200)
				c = 0
				for i, parts in enumerate(structs):
					parts.id = '{}-{}'.format(gene, i+1)
					parts = parts.link_part()
					print >> sys.stderr, 'old', parts.to_str()
					genefa = self.get_filename(self.hmmoutdir, parts, 'fa')
					with open(genefa, 'w') as fout:
						parts.write_seq(self.seqs, fout)
					output = self.get_filename(self.hmmoutdir, parts, 'trn')
					struct_file = output + '.struct'
					self.trnascan(genefa, output, opts='{} -Q -f {}'.format(self.trn_opts, struct_file))
					for trn, struct in izip(tRNAscan(output), tRNAscanStructs(struct_file)):
						if not trn.is_trn(gene.name):
							continue
						c += 1
						trna, new_parts = parts.map_coord(trn.to_exons())
						new_parts = new_parts.link_part()
						print >> sys.stderr, 'new', new_parts.to_str()
						new_parts.source = 'tRNAscan'
						rename = trn.update_name(gene.name)
						new_gene = copy.deepcopy(gene)
						if rename != gene.name:
							logger.info('renaming `{}` to `{}`'.format(gene.name, rename))
							new_gene.name = rename
						new_parts.id = '{}-{}'.format(gene, c)
						record = trna.extend_gene(new_gene, new_parts, rna_type=gene.seq_type)	 # GffExons
						#record.write(sys.stderr)
						rna_seq = record.extract_seq(self.seqs)
						record.rna_seq = rna_seq
						record.struct = struct
						record = self.score_record(record)
						record.npart = len(new_parts)
						record.write(sys.stderr)
					#	trn_records += [record]
						if record.trans_splicing:
							logger.warn('This tRNA is annotated as trans_splicing. Discarded')
							continue
						records += [record]
		#records += self.remove_duplicates(trn_records)
			#break
		return records


	def hmmsearch_protein(self):
		records = []
		if self.no_cds:
			return records
		# translate
		aa_seq = '{}/{}.genome.aa'.format(self.tmpdir, self.prefix)
		with open(aa_seq, 'w') as fout:
			d_length = six_frame_translate(self.fsa, fout, seqfmt='fasta', #self.seqfmt, 
											transl_table=self.transl_table)
			#print >>sys.stderr,d_length
		# hmmsearch
		for gene in self.db.cds_genes:
			print >>sys.stderr, '\n   >> {}: {}'.format(gene, gene.name)
			hmmfile = self.db.get_hmmfile(gene)
			domtblout = self.get_domtblout(gene, src='g')
			self.hmmsearch(hmmfile, aa_seq, domtblout)
			structs = HmmSearch(domtblout).get_gene_structure(d_length,
							min_ratio=self.cov_cutoff,
							min_hmmcov=self.min_cds_hmmcov,
							min_cov=self.min_cds_cov, seq_type='prot', flank=5000)
			for i, parts in enumerate(structs):
				parts.id = '{}-{}'.format(gene, i+1)	# parts is a copy
				parts = parts.link_part()
				genefa = self.get_filename(self.agtoutdir, parts, 'fa')
				with open(genefa, 'w') as fout:
					parts.write_seq(self.seqs, fout)
				if self.est is not None:
					self.exonerate_gene_est(genefa, gene, parts)
				# hints
				self.prepare_gene_hints(genefa, gene, parts)
				# gene by exonerate
				ex_gtf, ex_domtblout = self.exonerate_gene_predict(genefa, gene, parts)
				# gene by augustus
				ag_gtf, ag_domtblout = self.augustus_gene_predict(genefa, gene, parts)
				best_gtf, source = self.get_best_gene(ag_gtf, ag_domtblout, ex_gtf, ex_domtblout, id=gene.id)
				pseudo = False
				if best_gtf is None:
					pseudo = True
					ex_gtf, ex_domtblout = self.exonerate_gene_predict(genefa, gene, parts, completed=False)
					best_gtf, source = self.get_best_gene(ag_gtf, ag_domtblout, 
												ex_gtf, ex_domtblout, id=gene.id)
					if best_gtf is None:
						continue
			#	print >> sys.stderr, source
			#	best_gtf.write(sys.stderr)	# GffRecord -> AugustusGtfLine -> Gtf
			#	print >> sys.stderr, ''
				print >> sys.stderr, 'old', parts.to_str()
				cds, new_parts = parts.map_coord(best_gtf.to_exons().filter('CDS'))	# GffExons
				new_parts = new_parts.link_part()
				print >> sys.stderr, 'new', new_parts.to_str()
			#	print >> sys.stderr, parts.to_str()
				new_parts.source = source
			#	cds.write(sys.stderr)
				record = cds.extend_gene(gene, new_parts, rna_type=gene.seq_type, pseudo=pseudo)
				print >> sys.stderr, ''
				#record.write(sys.stderr)
				cds_seq = record.extract_seq(self.seqs)
				pep_seq = record.translate_cds(cds_seq, transl_table=self.transl_table)
				record.cds_seq, record.pep_seq = cds_seq, pep_seq
				record.pseudo = pseudo
			#	record.rna_seq = record.cds_seq
			#	print >> sys.stderr, '# coding sequence = [{}]'.format(cds_seq)
			#	print >> sys.stderr, '# protein sequence = [{}]'.format(pep_seq)
				#record = cds.to_gff_record()
				record = self.score_record(record)
				record.npart = len(new_parts)
				record.write(sys.stderr)
				records += [record]
			#break
			#print >>sys.stderr, '\n'
		return records
		
	def prepare_hints(self):
		for gene in self.db.cds_genes:
			self.prepare_gene_hints(self.fsa, gene)
	def prepare_gene_hints(self, reference, gene, copy=None):
		if copy is None:
			copy = gene
		seqfile = self.db.get_seqfile(gene)
		exn_gff = self.get_exnfile(copy, 'p')
		self.exonerate(seqfile, reference, exn_gff,		# protein
				model='protein2genome', percent=20,
				maxintron=500000,
				showtargetgff='T')
		hintfile = self.get_hintfile(copy)
		with open(hintfile, 'w') as fout:
			ExonerateGffGenes(exn_gff).to_hints(fout, src='P', pri=4)
			if self.est is not None:
				est_exn_gff = self.get_exnfile(copy, 'e')
				ExonerateGffGenes(est_exn_gff).to_hints(fout, src='E', pri=4)
		return exn_gff
	def exonerate_est2genome(self, reference, gene, copy=None, minintron=500):
		if copy is None:
			copy = gene
		seqfile = self.db.get_seqfile(gene)
		exn_gff = self.get_exnfile(copy, 'e')
		# self.exonerate(seqfile, reference, exn_gff,		# est
				# model='est2genome', bestn=5, percent=70, 
				# maxintron=10000,
				# minintron=200,
				# showtargetgff='T')
		self.blat(seqfile, reference, exn_gff, maxIntron=10000)
		
		outfa = exn_gff + '.fa'
		gene_seqs = self.get_seqs(reference)
		with open(outfa, 'w') as fout:
			# ex_exons = ExonerateGffGenes(exn_gff).to_exons(gene_seqs, fout)
			ex_exons = PslParser(exn_gff).to_exons(minintron=minintron, d_seqs=gene_seqs, fout=fout)
		hmmfile = self.db.get_hmmfile(gene)
		domtblout = outfa + '.domtbl'
		self.hmmsearch(hmmfile, outfa, domtblout)
		ex_hmm_best = HmmSearch(domtblout).get_best_hit(score=True) if os.path.exists(domtblout) else None
		if ex_hmm_best is None:
			ex_best = None
		else:
			ex_best = [record for record in ex_exons \
							if record.id == ex_hmm_best.tname][0]
		return ex_best
	
	def penalize_exon_diff(self, record, id):
		#id = record.id
		exon_count = record.count_type('cds', 'CDS') #count_exon()
		try: db_exon_count = self.db.gene_info[id].exon_count
		except KeyError as e:
			print >>sys.stderr, self.db.gene_info
			raise KeyError(e)
		diff = abs(exon_count - db_exon_count)
		return diff * self.exon_diff_penalty
		
	def get_best_gene(self, ag_gtf, ag_domtblout, ex_gtf, ex_domtblout, id=None, ex_weight=0.95):
		ag_hmm_best = HmmSearch(ag_domtblout).get_best_hit() if os.path.exists(ag_domtblout) else None
		ex_hmm_best = HmmSearch(ex_domtblout).get_best_hit() if os.path.exists(ex_domtblout) else None
		if ag_hmm_best is None:
			ag_best = None
		else:
			try:
				ag_best = [record for record in AugustusGtfGenes(ag_gtf) \
							if record.id == ag_hmm_best.tname][0]
				both_support = ag_best.annotations.supported == ag_best.annotations.total_exons \
						and ag_best.annotations.fully_obeyed > 0
			except IndexError:	# should not to here
				ex_hmm_best = None
				ag_best = None
		if ex_hmm_best is None:
			ex_best = None
		else:
			try:
				ex_best = [record for record in AugustusGtfGenes(ex_gtf) \
							if record.id == ex_hmm_best.tname][0]
			except IndexError:
				ex_hmm_best = None
				ex_best = None
		ag_best = (ag_best, 'augustus')
		ex_best = (ex_best, 'exonerate')
		none = (None, None)
		if ag_hmm_best is None and ex_hmm_best is None:	# both no hit
			return none
		elif ag_hmm_best is None:	# augustus no hit
			return ex_best
		elif ex_hmm_best is None:	# exonerate no hit
			if both_support:
				return ag_best
			else:
				return ag_best
		ex_hmm_best.edit_score -= self.penalize_exon_diff(ex_best[0], id)
		ag_hmm_best.edit_score -= self.penalize_exon_diff(ag_best[0], id)
		print >>sys.stderr, ex_hmm_best.edit_score, ag_hmm_best.edit_score
		if ex_hmm_best.edit_score*ex_weight > ag_hmm_best.edit_score:
			return ex_best
		else:
			if both_support:	# augustus must be support by both all exons and fully hints
				return ag_best	# strict?
			else:
				return ag_best	# change
	def exonerate_gene_predict(self, reference, gene, copy=None, completed=True):
		exn_gff = self.get_exnfile(copy, 'p')
		outgff = exn_gff + '.gff'
		gene_seqs = self.get_seqs(reference)
		with open(outgff, 'w') as fout:
			exons = ExonerateGffGenes(exn_gff).get_gene_gtf(gene_seqs, fout,
											  transl_table=self.transl_table)
		pepfaa = outgff + '.faa'
		with open(pepfaa, 'w') as fout:
			self.check_augustus_gff(outgff, fout, completed=completed)
		hmmfile = self.db.get_hmmfile(gene)
		domtblout = pepfaa + '.domtbl'
		self.hmmsearch(hmmfile, pepfaa, domtblout)
		return outgff, domtblout


	def augustus_predict(self):
		for gene in self.db.cds_genes:
			self.augustus_gene_predict(self.fsa, gene)
	def augustus_gene_predict(self, reference, gene, copy=None):
		if copy is None:
			copy = gene
		augusuts_species = self.db.get_augustus_species(gene)
		pflfile = self.db.get_pflfile(gene)
		hintfile = self.get_hintfile(copy)
		outgff = self.get_augustus_gff(copy)
		kargs = {'translation_table': self.transl_table,
				'hintsfile': hintfile, 'extrinsicCfgFile': 'extrinsic.MPE.cfg',
				'proteinprofile': pflfile, '/ExonModel/minexonlength': 20,
				'codingseq': 1, 'noInFrameStop': 1,
				}
		self.augustus(reference, augusuts_species, outgff, kargs=kargs)
		pepfaa = outgff + '.faa'
		with open(pepfaa, 'w') as fout:
			self.check_augustus_gff(outgff, fout)
		hmmfile = self.db.get_hmmfile(gene)
		domtblout = pepfaa + '.domtbl'
		self.hmmsearch(hmmfile, pepfaa, domtblout)
		return outgff, domtblout
	def check_augustus_gff(self, gff, fout, completed=True):
		genes, has_block, has_support, has_obey = 0,0,0,0
		full_support = 0
		both_support = 0
		for record in AugustusGtfGenes(gff):
			if completed and not record.is_complete:	# only use compelte gene
				continue
			genes += 1
			if record.annotations.blocks:
				has_block += 1
			if record.annotations.supported > 0:
				has_support += 1
			if record.annotations.supported == record.annotations.total_exons:
				full_support += 1
			if record.annotations.fully_obeyed > 0:
				has_obey += 1
			if record.annotations.supported == record.annotations.total_exons \
				and record.annotations.fully_obeyed > 0:
				both_support += 1
			seq = record.annotations.protein_sequence
			desc = 'block:{} CDS_exons:{}/{} P:{} E:{} fully_obeyed:{}'.format(
				len(record.annotations.blocks), record.annotations.supported,
				record.annotations.total_exons,
				record.annotations.supported_P, record.annotations.supported_E,
				record.annotations.fully_obeyed)
			if completed and {'X', '*'} & set(seq):	# stop codon in CDS
				continue
			print >>fout, '>{} {}\n{}'.format(record.id, desc, seq)
		return genes, has_block, has_support, full_support, has_obey, both_support
	def hmmsearch_est(self):
		aa_seq = '{}/{}.est.aa'.format(self.tmpdir, self.prefix)
		with open(aa_seq, 'w') as fout:
			d_length = six_frame_translate(self.est, fout, transl_table=self.transl_table)	
		for gene in self.db.cds_genes:
			hmmfile = self.db.get_hmmfile(gene)
			domtblout = self.get_domtblout(gene, src='e')
			self.hmmsearch(hmmfile, aa_seq, domtblout)
			est_seq = self.get_domtblfa(gene, src='e')
			HmmSearch(domtblout).get_hit_nuclseqs(self.est, est_seq)
			#self.exonerate_gene_est(self.genome, gene)
		os.remove(aa_seq)
	def exonerate_gene_est(self, reference, gene, copy=None):
		if copy is None:
			copy = gene
		est_seq = self.get_domtblfa(gene, src='e')
		exn_gff = self.get_exnfile(copy, 'e')
		self.exonerate(est_seq, reference, exn_gff, 		# est
				model='est2genome', bestn=5, percent=70, 
				maxintron=500000,
				geneticcode=self.transl_table,
				showtargetgff='T')
		
	def hmmsearch(self, hmmfile, seqdb, domtblout):
		cmd = 'hmmsearch --nobias --domtblout {domtblout} {hmmfile} {seqdb} > /dev/null'.format(
				hmmfile=hmmfile, seqdb=seqdb, domtblout=domtblout)
		run_cmd(cmd, log=True)
		return cmd
	def exonerate(self, queryfile, targetfile, outhit, **kargs):
		cmd = ['exonerate {query} {target}'.format(
						query=queryfile, target=targetfile)]
		for key, value in kargs.items():
			if value is not None:
				cmd += ['--{key} {value}'.format(key=key, value=value)]
		cmd += ['> {}'.format(outhit)]
		cmd = ' '.join(cmd)
		run_cmd(cmd, log=True)
		return cmd
	def augustus(self, queryfile, species, outgff, kargs={}):
		cmd = ['augustus --species={}'.format(species)]
		for key, value in kargs.items():
			if value is not None:
				cmd += ['--{key}={value}'.format(key=key, value=value)]
		cmd += [queryfile]
		cmd += ['> {}'.format(outgff)]
		cmd = ' '.join(cmd)
		run_cmd(cmd, log=True)
		return cmd
	def blat(self, database, query, output, **kargs):
		cmd = ['blat {} {} {}'.format(database, query, output)]
		for key, value in kargs.items():
			if value is not None:
				cmd += ['-{key}={value}'.format(key=key, value=value)]
		cmd = ' '.join(cmd)
		run_cmd(cmd, log=True)
		return cmd
		
	def trnascan(self, queryfile, output, opts='-O'):
		cmd = ['tRNAscan-SE {}'.format(queryfile)]
		cmd += [opts]
		cmd += ['> {}'.format(output)]
		cmd = ' '.join(cmd)
		run_cmd(cmd, log=True)
		return cmd

	def get_filename(self, _dir, gene, *field):
		return '{}/{}.{}'.format(_dir, gene, '.'.join(field))
	def get_domtblfa(self, gene, **kargs):
		return self.get_domtblout(gene, **kargs) + '.fa'
	def get_augustus_gff(self, gene):
		return '{}/{}.gff'.format(self.agtoutdir, gene)
	def get_domtblout(self, gene, src='g'):
		outdir = self.hmmoutdir
		return '{}/{}.{}.domtblout'.format(outdir, gene, src)
	def get_exnfile(self, gene, src='p'):
		return '{}/{}.{}.exgff'.format(self.agtoutdir, gene, src)
	def get_hintfile(self, gene):
		return '{}/{}.hints'.format(self.agtoutdir, gene)
		
	def get_seqs(self, seqfile, seqfmt='fasta'):
		return OrderedDict([(rc.id, rc.seq) for rc in SeqIO.parse(seqfile, seqfmt)])
	def double_seqs(self, seqfile, fout, seqfmt='fasta'):
		d_length = {}
		for rc in SeqIO.parse(seqfile, seqfmt):
			for seq, suffix0 in zip([rc.seq, rc.seq.reverse_complement()], ['fwd', 'rev']):
				suffix = '|{}'.format(suffix0)
				print >> fout, '>{}{}\n{}'.format(rc.id, suffix, seq)
			d_length[rc.id] = len(rc.seq)
		return d_length
		
	def to_fasta(self, records):
		
		cds_fa = '{}/{}.cds.fasta'.format(self.outdir, self.prefix)
		pep_fa = '{}/{}.pep.fasta'.format(self.outdir, self.prefix)
		rna_fa = '{}/{}.rna.fasta'.format(self.outdir, self.prefix)
		f_cds = open(cds_fa, 'w')
		f_pep = open(pep_fa, 'w')
		f_rna = open(rna_fa, 'w')
		ids = set([])
		d_count = {}
		for i, record in enumerate(sorted(records, key=lambda x:(x.name, -x.score))):
			id = record.name
			if id in d_count:
				xid = '{}-{}'.format(id, d_count[id] + 1)
			else:
				xid = id
			try: d_count[id] += 1
			except KeyError: d_count[id] = 1
			
			desc = 'gene={};id={};product={};exons={};score={};cov={}'.format(record.name, 
						record.rna_id, record.product, record, record.score, record.cov)
			if getattr(record, 'pseudo', None):
				desc += ';pseudo=true'
			try:
				print >> f_cds, '>{} {}\n{}'.format(xid, desc, record.cds_seq)
			except AttributeError: pass
			try:
				print >> f_pep, '>{} {}\n{}'.format(xid, desc, record.pep_seq)
			except AttributeError: pass
			try:
				print >> f_rna, '>{} {}\n{}'.format(xid, desc, record.rna_seq)
			except AttributeError: pass
		f_cds.close()
		f_pep.close()
		f_rna.close()

	def to_gff3(self, records):
		gff = '{}/{}.gff3'.format(self.outdir, self.prefix)
		with open(gff, 'w') as fout:
			print >> fout, '##gff-version 3'
			for record in records:
				record.write(fout)
				try:
					print >> fout, '# coding sequence = [{}]'.format(record.cds_seq)
					print >> fout, '# protein sequence = [{}]'.format(record.pep_seq)
				except AttributeError: pass
	def to_fsa(self):
		desc = []
		if self.organism is not None:
			organism = self.organism
			desc += ['[organism={}]'.format(organism)]
			
		desc += ['[location={}]'.format(LOCATION[self.organ])]
		if not self.linear:
			desc += ['[topology=circular]']
#		else:
#			desc += ['[topology=linear]']
		if not self.partial:
			desc += ['[completeness=complete]']
#		else:
#			desc += ['[completeness=partial]']
		if self.trans:
			desc += ['[moltype=transcribed_RNA]']
		desc += ['gcode={}'.format(self.transl_table)]	# transl_table in tbl do not work well

		desc = ' '.join(desc)
		fsa = '{}/{}.fsa'.format(self.outdir, self.prefix)
		fout = open(fsa, 'w')
		for id, seq in self.seqs.items():	
			print >> fout, '>{} {}\n{}'.format(id, desc, seq)
		fout.close()
		return fsa

	def re_fsa(self, records):
		d_seqs = OrderedDict([(rc.id, rc) for rc in SeqIO.parse(self.fsa, 'fasta')])
		chroms = {exon.chrom for record in records for exon in record}
		fout = open(self.fsa, 'w')
		for chrom, rc in d_seqs.items():
			if chrom in chroms:
				SeqIO.write(rc, fout, 'fasta')
		fout.close()

	def to_sqn(self, records):
		g,s = self.organism.split()[:2]
		prefix = g[:1] + s[:2] + self.organ #[0]
		prefix = prefix.upper()
		# to tbl
		tbl = '{}/{}.tbl'.format(self.outdir, self.prefix)
		fout = open(tbl, 'w')
		d_tags = {}
		i = 0
		for seqid in self.seqs.keys():
			my_records = [record for record in records if seqid in record.chroms]
			if not my_records:
				continue
			print >>fout, '>Feature {}'.format(seqid)
			for record in my_records:
				if record in d_tags:
					locus_tag = d_tags[record]
				else:
					i += 1
					locus_tag = '{}_{:0>4d}'.format(prefix, i)
					d_tags[record] = locus_tag
				record.to_tbl(fout, seqid, transl_table=self.transl_table, 
						locus_tag=locus_tag, wgs=self.wgs)
		fout.close()
		# fsa + tbl -> sqn -> genbank
		sqn = '{}/{}.sqn'.format(self.outdir, self.prefix)
		templete = self.db.templete
		cmd = 'tbl2asn -i {} -o {} -t {}'.format(self.fsa, sqn, templete)
		if self.nseqs > 1:
			cmd += ' -M b'
		run_cmd(cmd, log=True)
		gb = '{}/{}.gb'.format(self.outdir, self.prefix)
		cmd = 'asn2gb -i {} -o {}'.format(sqn, gb)
		run_cmd(cmd, log=True)
	def draw_gene_map(self):
		gb = '{}/{}.gb'.format(self.outdir, self.prefix)
		outfig = '{}/{}.map.pdf'.format(self.outdir, self.prefix)
		self._draw_map(gb, outfig)
		
	def _draw_map(self, ingb, outfig, opts='--density 300'):
		tmpfig = '{}/{}.ps'.format(self.tmpdir, os.path.basename(outfig))
		cmd = 'drawgenemap --infile {} --format ps --outfile {} {}'.format(ingb, tmpfig, opts)
		run_cmd(cmd, log=True)
		cmd = 'ps2pdf {} {}'.format(tmpfig, outfig)
		run_cmd(cmd, log=True)
		
	def compare_gene_map(self):
		opts = '--force_linear'
		prefix = self.prefix.replace('.', '_')
		gb = self.genome
		outfig_ref = '{}/{}_ref.pdf'.format(self.tmpdir, prefix)
		self._draw_map(gb, outfig_ref, opts=opts)
		gb = '{}/{}.gb'.format(self.outdir, self.prefix)
		outfig_cmp = '{}/{}_cmp.pdf'.format(self.tmpdir, prefix)
		self._draw_map(gb, outfig_cmp, opts=opts)
		acc = '/'.join(self.seqs.keys()).replace('_', '\_')
		tex = r'''%% !Mode:: "TeX:UTF-8:Hard"
\documentclass[a4paper]{article}
\usepackage[margin=3mm]{geometry}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{subfig}
\usepackage{float}
%%\usepackage{pdflscape}
\begin{document}
%%\begin{landscape}
\begin{figure}[!htb]
  \centering
  \subfloat[Reference]{\includegraphics[angle=-90,totalheight=270mm]{%s}}
  \subfloat[OGAP]{\includegraphics[angle=-90,totalheight=270mm]{%s}}
  \caption{%s %s / %s (%s bp)}
\end{figure}
%%\end{landscape}
\end{document}
''' % (outfig_ref, outfig_cmp, self.organism, LOCATION[self.organ], acc, self.seqlen)
		texfile = '{}/{}.merge.tex'.format(self.tmpdir, self.prefix)
		with open(texfile, 'w') as f:
			print >>f, tex
		basename = os.path.splitext(os.path.basename(texfile))[0]
		cmd = 'cd {} && pdflatex {} && rm {basename}.aux {basename}.log'.format(
				self.outdir, texfile, basename=basename)
		run_cmd(cmd, log=True)
	
	def plot_struct(self, trndir, records):
		mkdirs(trndir)
		ids = set([])
		for record in sorted(records, key=lambda x:(x.name, -x.score)):
			if not record.rna_type == 'tRNA':
				continue
			id = '{}-{}'.format(record.product, record.name)
			if id in ids:
				print >>sys.stderr, '{} duplicates, ignored'.format(id)
				continue
			ids.add(id)
			struct_file = '{}/{}.struct'.format(trndir, id)
			with open(struct_file, 'w') as fout:
				print >> fout, '{seq}\n{struct}'.format(
					seq=record.struct.rna_seq, struct=record.struct.rna_struct)
			struct_file = os.path.realpath(struct_file)
			outfig = '{}/{}.pdf'.format(trndir, id)
			# cmd = 'cd {dir} && cat {struct_file} | \
					# RNAplot --pre "{start} {end} 8 GREEN omark" && mv rna.ps {id}.ps'.format(
					# dir=trndir, struct_file=struct_file, 
					# start=record.struct.codon_start, end=record.struct.codon_end,
					# id=id)
			cmd = 'cd {dir} && cat {struct_file} | \
					RNAplot --pre "{start} {end} 8 GREEN omark" && ps2pdf rna.ps {outfig}'.format(
					dir=self.tmpdir, struct_file=struct_file, 
					start=record.struct.codon_start, end=record.struct.codon_end,
					outfig=outfig)
			run_cmd(cmd, log=True)

	def _check_database(self):
		pass
	def _check_dependencies(self):
		pass
	def guess_seqfmt(self, seqfile):
		head = open(seqfile).read(10)
		if head.startswith('>'):
			return 'fasta'
		elif head.startswith('LOCUS'):
			return 'genbank'
		else:
			raise ValueError('unrecognized sequence format')
	def get_record(self, genome):
		for rc in SeqIO.parse(genome, self.seqfmt):
			return rc
	def contains_gap(self, seqs):
		for seq in seqs:
			if set(seq.upper()) - set('ATCG'):
				return True
		return False

def main():
	args = makeArgparse()
	print >>sys.stderr, args.__dict__
	pipeline = Pipeline(**args.__dict__)
	pipeline.run()


if __name__ == '__main__':
	main()

