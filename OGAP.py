import sys, os
import copy
import argparse
from collections import OrderedDict
from itertools import izip
from Bio import SeqIO

from lib.Database import Database
from lib.Hmmer import HmmSearch
from lib.tRNA import tRNAscan, tRNAscanStructs
from lib.Gff import ExonerateGffGenes, AugustusGtfGenes
from lib.RunCmdsMP import run_cmd, run_job, logger
from lib.translate_seq import six_frame_translate
from lib.small_tools import mkdirs, rmdirs

LOCATION = {'pt': 'chloroplast', 'mt': 'mitochondrion'}
class Pipeline():
	def __init__(self, genome, 
				organ, taxon, 
				transl_table=None,
				est=None, # EST evidence
				prefix=None, tmpdir='/dev/shm/tmp', 
				organism=None,
				linear=False, 
				partial=False,
				include_orf=False,
				min_cds_cov=80,
				min_rrn_cov=90,
				min_trn_cov=50,
				score_cutoff=0.9,
				cov_cutoff=0.9,
				trn_struct=False,
				cleanup=True,
				seqfmt='fasta', **kargs):
		
		self.organ = organ
		self.genome = os.path.realpath(genome)
		self.organ = organ
		self.taxon = taxon
		self.est = est
		self.organism = organism
		self.linear = linear
		self.partial = partial
		self.include_orf = include_orf
		self.trn_struct = trn_struct
		self.cleanup = cleanup
		if seqfmt is None:
			self.seqfmt = self.guess_seqfmt(self.genome)
			logger.info('genome in {} format'.format(self.seqfmt))
		else:
			self.seqfmt = seqfmt
		if self.seqfmt == 'genbank':
			rc = self.get_record()
			if self.organism is None:
				try:
					self.organism = rc.annotations['organism']
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

		# min coverage
		self.min_cds_cov = min_cds_cov
		self.min_rrn_cov = min_rrn_cov
		self.min_trn_cov = min_trn_cov
		self.score_cutoff = score_cutoff
		self.cov_cutoff = cov_cutoff

		if taxon is None and self.organism is not None:
			taxon = Database(organ=organ).select_db(self.organism).taxon
			logger.info('automatically select db: {}'.format(taxon))
		elif taxon is None:
			logger.info('neither -taxon or -organism must be specified')
		self.db = Database(organ=organ, taxon=taxon, include_orf=include_orf)
		
		self.ogtype = self.db.ogtype
		self.transl_table = transl_table
		
		if prefix is None:
			self.prefix = os.path.basename(genome)
		else:
			self.prefix = prefix
		# folder
		self.tmpdir = '{}/{}/{}'.format(tmpdir, self.ogtype, 
								os.path.basename(self.prefix))
		self.hmmoutdir = '{}/hmmout'.format(self.tmpdir)
		#self.estoutdir = '{}/estout'.format(self.tmpdir)
		#self.exnoutdir = '{}/exnout'.format(self.tmpdir)
		self.agtoutdir = '{}/augustus'.format(self.tmpdir)
		self.trndir = '{}/{}.trna'.format('.', self.prefix)

	def run(self):
		rmdirs(self.agtoutdir)
		mkdirs(self.tmpdir)
		mkdirs(self.hmmoutdir, self.agtoutdir)
		# check
		logger.info('checking database: {}'.format(self.db.ogtype))
		self.db.checkdb()
		if self.transl_table is None:
			self.transl_table = self.db.transl_table
		# read genome seqs
		self.seqs = self.get_seqs(self.genome, self.seqfmt)
		seqids = self.seqs.keys()
		# to fasta
		self.fsa = self.to_fsa()

		records = []
		# cds-protein finder
		if self.est is not None:
			self.hmmsearch_est()
		logger.info('finding protein-coding genes by HMM + enonerate + augustus')
		records += self.hmmsearch_protein()
		
		# rna finder
		logger.info('finding non-coding genes by HMM (rRNA) or HMM + tRNAscan-SE(tRNA)')
		records += self.hmmsearch_rna()

		# score
		logger.info('scoring genes by HMM')
		records = [self.score_record(record) for record in records]
		# remove duplicates
		logger.info('removing duplicated genes with the identical coordinate')
		records = self.remove_duplicates(records)
		# remove low quality of the same gene
		logger.info('removing duplicated genes with lower score: {} * highest score'.format(self.score_cutoff))
		records = self.remove_lowqual(records, cutoff=self.score_cutoff)
		
		# sort by coordinate
		records = sorted(records, key=lambda x: (seqids.index(x.chrom), x.start))
		# to gff
		self.to_gff3(records)

		# out sqn
		logger.info('outputing sqn for submitting Genbank')
		self.to_sqn(records)

		# out fasta of gene, sorted by score
		self.to_fasta(records)
	
		# tRNA structure
		if self.trn_struct:
			self.plot_struct(self.trndir, records)
	
		# summary
		self.summary_records(records)
		
		# clean up
		if self.cleanup:
			logger.info('cleaning')
			rmdirs(self.tmpdir)
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
			cmd = 'cd {dir} && cat {struct_file} | \
					RNAplot --pre "{start} {end} 8 GREEN omark" && mv rna.ps {id}.ps'.format(
					dir=trndir, struct_file=struct_file, 
					start=record.struct.codon_start, end=record.struct.codon_end,
					id=id)
			run_cmd(cmd, log=True)

	def remove_lowqual(self, records, cutoff=0.9):
		d_group = OrderedDict()
		for record in records:
			key = record.id
			try: d_group[key] += [record]	# duplicates from homologous gene
			except KeyError: d_group[key] = [record]
		better_records = []
		for key, records in d_group.items():
			if len(records) == 1:	# unique
				better_records += records
				continue
			highest_score = max([record.score for record in records])
			good_cutoff = highest_score * cutoff
			better_record = [record for record in records if record.score >= good_cutoff]
			print >>sys.stderr, key, len(records), '->', len(better_record)
			better_records += better_record
		return better_records
	def remove_duplicates(self, records):
		keys = set([])
		d_group = OrderedDict()
		for record in records:
			key1 = (str(record),)
			key2 = key1 + (record.id, )
			if key2 in keys:	# duplicates from the same gene
				print >>sys.stderr, key2, 'removed'
				continue
			try: d_group[key1] += [record]	# duplicates from homologous gene
			except KeyError: d_group[key1] = [record]
			keys.add(key2)
		unique_records = []
		for key, records in d_group.items():
			if len(records) == 1:	# unique
				unique_records += records
				continue
			best_record = max(records, key=lambda x:x.score)
			print >>sys.stderr, key, best_record.name, len(records), '->', 1
			unique_records += [best_record]
		return unique_records
		
	def score_record(self, record):
		rnafa = self.get_filename(self.hmmoutdir, record.gene_id, 'fasta')
		with open(rnafa, 'w') as fout:
			try:
				print >> fout, '>{}\n{}'.format(record.gene_id, record.pep_seq)
			except:
				print >> fout, '>{}\n{}'.format(record.gene_id, record.rna_seq)
		hmmfile = self.db.get_hmmfile(record.id)
		domtblout = rnafa + '.domtbl'
		self.hmmsearch(hmmfile, rnafa, domtblout)
		hmm_best = HmmSearch(domtblout).get_best_hit()
		record.score = round(hmm_best.edit_score, 1)
		record[0].score = record.score
		return record
		
	def summary_records(self, records):
		d_smy = OrderedDict()
		for record in records:
			try: d_smy[record.rna_type] += [record]
			except KeyError: d_smy[record.rna_type] = [record]
			
		line = ['type', 'copy number', 'gene number', 'gene names']
		print >>sys.stderr, '\t'.join(line)
		for rna_type, records in d_smy.items():
			genes = [record.id for record in records]
			names = list({record.name for record in records})
			names = sorted(names)
			line = [rna_type, len(genes), len(names), ','.join(names)]
			line = map(str, line)
			print >>sys.stderr, '\t'.join(line)
	
	def hmmsearch_rna(self):
		na_seq = '{}/{}.genome.na'.format(self.tmpdir, self.prefix)
		na_seq = '{}/{}.genome.na'.format(self.tmpdir, self.prefix)
		with open(na_seq, 'w') as fout:
			d_length = self.double_seqs(self.genome, fout, seqfmt=self.seqfmt)
		records = []
		#trn_records = []
		for gene in self.db.rna_genes:
			print >>sys.stderr, '	>>', gene
#			if not gene.seq_type == 'rRNA':
#				continue
			hmmfile = self.db.get_hmmfile(gene)
			domtblout = self.get_domtblout(gene, src='g')
			self.hmmsearch(hmmfile, na_seq, domtblout)
			if gene.seq_type == 'rRNA':
				structs = HmmSearch(domtblout).get_gene_structure(d_length, 
								min_ratio=self.cov_cutoff,
								min_cov=self.min_rrn_cov, seq_type='nucl', flank=0)
				for i, parts in enumerate(structs):
					parts.id = '{}-{}'.format(gene, i+1)
					parts.source = 'hmmsearch'
					exons = parts.to_exons()
					#exons.write(sys.stderr)
					record = exons.extend_gene(gene, parts, rna_type=gene.seq_type)
					print >> sys.stderr, ''
					record.write(sys.stderr)
					rna_seq = record.extract_seq(self.seqs[record.chrom])
					record.rna_seq = rna_seq
					record = self.score_record(record)
					records += [record]
			elif gene.seq_type == 'tRNA':
				structs = HmmSearch(domtblout).get_gene_structure(d_length,
								min_ratio=self.cov_cutoff,
							    min_cov=self.min_trn_cov, seq_type='nucl', flank=200)
				for i, parts in enumerate(structs):
					parts.id = '{}-{}'.format(gene, i+1)
					genefa = self.get_filename(self.hmmoutdir, parts, 'fa')
					with open(genefa, 'w') as fout:
						parts.write_seq(self.seqs, fout)
					output = self.get_filename(self.hmmoutdir, parts, 'trn')
					struct_file = output + '.struct'
					self.trnascan(genefa, output, opts='-O -Q -f {}'.format(struct_file))
					for trn, struct in izip(tRNAscan(output), tRNAscanStructs(struct_file)):
						if not trn.is_trn(gene.name):
							continue
						trna, new_parts = parts.map_coord(trn.to_exons())
						new_parts.source = 'tRNAscan'
						rename = trn.update_name(gene.name)
						new_gene = copy.deepcopy(gene)
						if rename != gene.name:
							logger.info('renaming `{}` to `{}`'.format(gene.name, rename))
							new_gene.name = rename
						record = trna.extend_gene(new_gene, new_parts, rna_type='tRNA')	 # GffExons
						record.write(sys.stderr)
						rna_seq = record.extract_seq(self.seqs[record.chrom])
						record.rna_seq = rna_seq
						record.struct = struct
					#	record = self.score_record(record)
					#	trn_records += [record]
						records += [record]
		#records += self.remove_duplicates(trn_records)
		return records


	def hmmsearch_protein(self):
		# translate
		aa_seq = '{}/{}.genome.aa'.format(self.tmpdir, self.prefix)
		with open(aa_seq, 'w') as fout:
			d_length = six_frame_translate(self.genome, fout, seqfmt=self.seqfmt, 
											transl_table=self.transl_table)
		# hmmsearch
		records = []
		for gene in self.db.cds_genes:
			print >>sys.stderr, '   >>', gene
			hmmfile = self.db.get_hmmfile(gene)
			domtblout = self.get_domtblout(gene, src='g')
			self.hmmsearch(hmmfile, aa_seq, domtblout)
			structs = HmmSearch(domtblout).get_gene_structure(d_length,
							min_ratio=self.cov_cutoff,
							min_cov=self.min_cds_cov, seq_type='prot', flank=1000)
			for i, parts in enumerate(structs):
				parts.id = '{}-{}'.format(gene, i+1)	# parts is a copy
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
				best_gtf, source = self.get_best_gene(ag_gtf, ag_domtblout, ex_gtf, ex_domtblout)
				if best_gtf is None:
					continue
			#	print >> sys.stderr, source
			#	best_gtf.write(sys.stderr)	# GffRecord -> AugustusGtfLine -> Gtf
			#	print >> sys.stderr, ''
				cds, new_parts = parts.map_coord(best_gtf.to_exons().filter('CDS'))	# GffExons
			#	print >> sys.stderr, parts.to_str()
				new_parts.source = source
			#	cds.write(sys.stderr)
				record = cds.extend_gene(gene, new_parts, rna_type='mRNA')
				print >> sys.stderr, ''
				record.write(sys.stderr)
				cds_seq = record.extract_seq(self.seqs[record.chrom])
				pep_seq = record.translate_cds(cds_seq, transl_table=self.transl_table)
				record.cds_seq, record.pep_seq = cds_seq, pep_seq
			#	record.rna_seq = record.cds_seq
			#	print >> sys.stderr, '# coding sequence = [{}]'.format(cds_seq)
			#	print >> sys.stderr, '# protein sequence = [{}]'.format(pep_seq)
				#record = cds.to_gff_record()
				#record = self.score_record(record)
				records += [record]
			print >>sys.stderr, '\n'
		return records
		
	def prepare_hints(self):
		for gene in self.db.cds_genes:
			self.prepare_gene_hints(self.genome, gene)
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

	def get_best_gene(self, ag_gtf, ag_domtblout, ex_gtf, ex_domtblout, ex_weight=0.99):
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
	#	print >>sys.stderr, ex_hmm_best.edit_score, ag_hmm_best.edit_score
		if ex_hmm_best.edit_score*ex_weight > ag_hmm_best.edit_score:
			return ex_best
		else:
			if both_support:	# augustus must be support by both all exons and fully hints
				return ag_best	# strict?
			else:
				return ag_best	# change
	def exonerate_gene_predict(self, reference, gene, copy=None):
		exn_gff = self.get_exnfile(copy, 'p')
		outgff = exn_gff + '.gff'
		gene_seqs = self.get_seqs(reference)
		with open(outgff, 'w') as fout:
			exons = ExonerateGffGenes(exn_gff).get_gene_gtf(gene_seqs, fout,
											  transl_table=self.transl_table)
		pepfaa = outgff + '.faa'
		with open(pepfaa, 'w') as fout:
			self.check_augustus_gff(outgff, fout)
		hmmfile = self.db.get_hmmfile(gene)
		domtblout = pepfaa + '.domtbl'
		self.hmmsearch(hmmfile, pepfaa, domtblout)
		return outgff, domtblout


	def augustus_predict(self):
		for gene in self.db.cds_genes:
			self.augustus_gene_predict(self.genome, gene)
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
			print >>sys.stderr, copy, self.check_augustus_gff(outgff, fout)
		hmmfile = self.db.get_hmmfile(gene)
		domtblout = pepfaa + '.domtbl'
		self.hmmsearch(hmmfile, pepfaa, domtblout)
		return outgff, domtblout
	def check_augustus_gff(self, gff, fout):
		genes, has_block, has_support, has_obey = 0,0,0,0
		full_support = 0
		both_support = 0
		for record in AugustusGtfGenes(gff):
			if not record.is_complete:	# only use compelte gene
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
			if {'X', '*'} & set(seq):	# stop codon in CDS
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
		cmd = 'hmmsearch --domtblout {domtblout} {hmmfile} {seqdb} > /dev/null'.format(
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
		
		cds_fa = '{}/{}.cds.fasta'.format('.', self.prefix)
		pep_fa = '{}/{}.pep.fasta'.format('.', self.prefix)
		rna_fa = '{}/{}.rna.fasta'.format('.', self.prefix)
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
			
			desc = 'gene={};id={};product={};exons={};score={}'.format(record.name, 
						record.rna_id, record.product, record, record.score)
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
		gff = '{}/{}.gff3'.format('.', self.prefix)
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
		if not self.partial:
			desc += ['[completeness=complete]']
		desc = ' '.join(desc)
		fsa = '{}/{}.fsa'.format('.', self.prefix)
		fout = open(fsa, 'w')
		for id, seq in self.seqs.items():
			
			print >> fout, '>{} {}\n{}'.format(id, desc, seq)
		fout.close()
		return fsa

	def to_sqn(self, records):
		# to tbl
		tbl = '{}/{}.tbl'.format('.', self.prefix)
		fout = open(tbl, 'w')
		for seqid in self.seqs.keys():
			print >>fout, '>Feature {}'.format(seqid)
			for record in records:
				record.to_tbl(fout, seqid, transl_table=self.transl_table)
		fout.close()
		# fsa + tbl -> sqn -> genbank
		sqn = '{}/{}.sqn'.format('.', self.prefix)
		templete = self.db.templete
		cmd = 'tbl2asn -i {} -o {} -t {}'.format(self.fsa, sqn, templete)
		run_cmd(cmd, log=True)
		gb = '{}/{}.gb'.format('.', self.prefix)
		cmd = 'asn2gb -i {} -o {}'.format(sqn, gb)
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
	def get_record(self):
		for rc in SeqIO.parse(self.genome, self.seqfmt):
			return rc

def main():
	args = makeArgparse()
	print >>sys.stderr, args.__dict__
	pipeline = Pipeline(**args.__dict__)
	pipeline.run()

def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("genome", action="store",type=str,
					help="input genome sequence in fasta or genbank format [required]")
	parser.add_argument('-organ', type=str, choices=['mt', 'pt'], default='mt',
					help="organelle type: mt (mitochondrion) or pt (plastid) [default=%(default)s]")
	parser.add_argument('-taxon', type=str, default=None, 
					help="taxon to filter out, such as Embryophyta [default=auto by -organism]")
	parser.add_argument("-tmpdir", action="store",
					default='/dev/shm/tmp', type=str,
					help="temporary directory [default=%(default)s]")
	parser.add_argument('-cleanup', action="store_true", default=False,
					help='clean up temporary directory')
	parser.add_argument("-prefix", action="store",
					default=None, type=str,
					help="output prefix [default='genome']")

	parser.add_argument('-seqfmt', type=str, choices=['fasta', 'genbank'], default=None,
					help="genome seqence format [default=auto]")
	parser.add_argument('-include_orf', action="store_true", default=False,
					help="include orf in gene predictione [default=%(default)s]")
	parser.add_argument("-est", action="store",type=str,
					help="EST sequences for evidence in fasta format")

	parser.add_argument('-organism', type=str, default=None,
					help="organism to be included in sqn [default=%(default)s]")
	parser.add_argument('-linear', action="store_true", default=False,
					help="topology to be included in sqn [default=circular]")
	parser.add_argument('-partial', action="store_true", default=False,
					help="completeness to be included in sqn [default=complete]")
	
	parser.add_argument('-min_cds_cov', type=float, default=80,
					help="min coverage to filter candidate coding genes [default=%(default)s]")
	parser.add_argument('-min_rrn_cov', type=float, default=90,
					help="min coverage to filter candidate rRNA genes [default=%(default)s]")
	parser.add_argument('-min_trn_cov', type=float, default=50,
					help="min coverage to filter candidate tRNA genes [default=%(default)s]")
	parser.add_argument('-score_cutoff', type=float, default=0.9,
					help="min score to output [default=%(default)s]")
	parser.add_argument('-trn_struct', action="store_true", default=False,
					help="output tRNA structure [default=%(default)s]")
	args = parser.parse_args()
	if args.organism is not None:
		args.organism = args.organism.replace('_', ' ')
	return args

if __name__ == '__main__':
	main()

