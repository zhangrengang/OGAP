import sys, os
import argparse
from collections import OrderedDict
from Bio import SeqIO
from lib.Database import Database
from lib.Hmmer import HmmSearch
from Gff import ExonerateGffGenes, AugustusGtfGenes
from RunCmdsMP import run_cmd, run_job, logger
from translate_seq import six_frame_translate
from small_tools import mkdirs, rmdirs

LOCATION = {'pt': 'chloroplast', 'mt': 'mitochondrion'}
class Pipeline():
	def __init__(self, genome, 
				organ, taxon, 
				est=None, # EST evidence
				prefix=None, tmpdir='/dev/shm/tmp', 
				organism=None,
				linear=False, 
				partial=False,
				include_orf=False,
				max_evalue=1e-5, 
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
		self.seqfmt = seqfmt

		self.db = Database(organ=organ, taxon=taxon, include_orf=include_orf)
		self.ogtype = self.db.ogtype
		if prefix is None:
			self.prefix = os.path.basename(genome)
		else:
			self.prefix = prefix
		self.tmpdir = '{}/{}'.format(tmpdir, self.ogtype)
		self.hmmoutdir = '{}/hmmout'.format(self.tmpdir)
		self.estoutdir = '{}/estout'.format(self.tmpdir)
		self.exnoutdir = '{}/exnout'.format(self.tmpdir)
		self.agtoutdir = '{}/augustus'.format(self.tmpdir)

	def run(self):
		rmdirs(self.agtoutdir)
		mkdirs(self.tmpdir)
		mkdirs(self.hmmoutdir, self.estoutdir, self.agtoutdir)
		# check
		logger.info('checking database: {}'.format(self.db.ogtype))
		self.db.checkdb()
		self.seqs = self.get_seqs(self.genome, self.seqfmt)
		seqids = self.seqs.keys()
		# est hmmsearch for evidence
	#	if self.est is not None:
	#		logger.info('HMMsearch EST')
	#		self.hmmsearch_est()
		# hints
	#	logger.info('prepare hints from protein [and EST] sequences')
		#self.prepare_hints()
		records = self.hmmsearch_rrn()
		#return
		# cds-protein hmmsearch
		logger.info('HMMsearch protein-encoding genes')
#		records = self.hmmsearch_protein()
		# cds augustus
	#	logger.info('using AUGUSTUS to predict protein-encoding genes')
		#self.augustus_predict()
		# cds integrate
		
		# rna hmmsearch
		
		# tRNAScan

		# rna integrate
		# sort
		records = sorted(records, key=lambda x: (seqids.index(x.chrom), x.start))
		# to gff
		self.to_gff3(records)

		# out sqn
		self.to_sqn(records)
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

	def to_sqn(self, records):
		# to tbl
		tbl = '{}/{}.tbl'.format('.', self.prefix)
		fout = open(tbl, 'w')
		for seqid in self.seqs.keys():
			print >>fout, '>Feature {}'.format(seqid)
			for record in records:
				record.to_tbl(fout, seqid, transl_table=self.db.transl_table)
		fout.close()
		# to fsa
		desc = []
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
		sqn = '{}/{}.sqn'.format('.', self.prefix)
		templete = self.db.templete
		cmd = 'tbl2asn -i {} -o {} -t {}'.format(fsa, sqn, templete)
		run_cmd(cmd, log=True)
		gb = '{}/{}.gb'.format('.', self.prefix)
		cmd = 'asn2gb -i {} -o {}'.format(sqn, gb)
		run_cmd(cmd, log=True)

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
	def get_exnfile(self, gene, src='p'):
		return '{}/{}.{}.exgff'.format(self.agtoutdir, gene, src)
	def get_hintfile(self, gene):
		return '{}/{}.hints'.format(self.agtoutdir, gene)
	def hmmsearch_est(self):
		aa_seq = '{}/{}.est.aa'.format(self.tmpdir, self.prefix)
		with open(aa_seq, 'w') as fout:
			d_length = six_frame_translate(self.est, fout, transl_table=self.db.transl_table)	
		for gene in self.db.cds_genes:
			hmmfile = self.db.get_hmmfile(gene)
			domtblout = self.get_domtblout(gene, src='e')
			self.hmmsearch(hmmfile, aa_seq, domtblout)
			est_seq = self.get_domtblfa(gene, src='e')
			HmmSearch(domtblout).get_hit_nuclseqs(self.est, est_seq)
			self.exonerate_gene_est(self.genome, gene)
		os.remove(aa_seq)
	def exonerate_gene_est(self, reference, gene, copy=None):
		if copy is None:
			copy = gene
		est_seq = self.get_domtblfa(gene, src='e')
		exn_gff = self.get_exnfile(copy, 'e')
		self.exonerate(est_seq, reference, exn_gff, 		# est
				model='est2genome', bestn=5, percent=70, 
				maxintron=500000,
				geneticcode=self.db.transl_table,
				showtargetgff='T')
	def hmmsearch_rrn(self):
		na_seq = '{}/{}.genome.na'.format(self.tmpdir, self.prefix)
		with open(na_seq, 'w') as fout:
			d_length = self.double_seqs(self.genome, fout, seqfmt=self.seqfmt)
		records = []
		for gene in self.db.rna_genes:
#			if not gene.seq_type == 'rRNA':
#				continue
			hmmfile = self.db.get_hmmfile(gene)
			domtblout = self.get_domtblout(gene, src='g')
			self.hmmsearch(hmmfile, na_seq, domtblout)
			structs = HmmSearch(domtblout).get_gene_structure(d_length, 
									min_cov=90, seq_type='nucl', flank=0)
			for i, parts in enumerate(structs):
				parts.id = '{}-{}'.format(gene, i+1)
				parts.source = 'hmmsearch'
				exons = parts.to_exons()
				#exons.write(sys.stderr)
				record = exons.extend_gene(gene, parts, rna_type=gene.seq_type)
				print >> sys.stderr, ''
				record.write(sys.stderr)
				records += [record]
		return records
	def hmmsearch_protein(self):
		# translate
		aa_seq = '{}/{}.genome.aa'.format(self.tmpdir, self.prefix)
		with open(aa_seq, 'w') as fout:
			d_length = six_frame_translate(self.genome, fout, seqfmt=self.seqfmt, 
											transl_table=self.db.transl_table)
		d_seqs = self.get_seqs(self.genome, self.seqfmt)
		# hmmsearch
		records = []
		for gene in self.db.cds_genes:
			print >>sys.stderr, '    >>', gene
			hmmfile = self.db.get_hmmfile(gene)
			domtblout = self.get_domtblout(gene, src='g')
			self.hmmsearch(hmmfile, aa_seq, domtblout)
			for i, parts in enumerate(HmmSearch(domtblout).get_gene_structure(d_length)):
				parts.id = '{}-{}'.format(gene, i+1)	# parts is a copy
				genefa = self.get_filename(self.agtoutdir, parts, 'fa')
				with open(genefa, 'w') as fout:
					parts.write_seq(d_seqs, fout)
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
				record = cds.extend_gene(gene, new_parts)
				print >> sys.stderr, ''
				record.write(sys.stderr)
				cds_seq = record.extract_seq(self.seqs[record.chrom])
				pep_seq = record.translate_cds(cds_seq, transl_table=self.db.transl_table)
				record.cds_seq, record.pep_seq = cds_seq, pep_seq
			#	print >> sys.stderr, '# coding sequence = [{}]'.format(cds_seq)
			#	print >> sys.stderr, '# protein sequence = [{}]'.format(pep_seq)
				#record = cds.to_gff_record()
				records += [record]
			print >>sys.stderr, '\n'
		return records
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
                                              transl_table=self.db.transl_table)
		pepfaa = outgff + '.faa'
		with open(pepfaa, 'w') as fout:
			self.check_augustus_gff(outgff, fout)
		hmmfile = self.db.get_hmmfile(gene)
		domtblout = pepfaa + '.domtbl'
		self.hmmsearch(hmmfile, pepfaa, domtblout)
		return outgff, domtblout

	def get_seqs(self, seqfile, seqfmt='fasta'):
		return OrderedDict([(rc.id, rc.seq) for rc in SeqIO.parse(seqfile, seqfmt)])
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
		kargs = {'translation_table': self.db.transl_table,
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
	def get_augustus_gff(self, gene):
		return '{}/{}.gff'.format(self.agtoutdir, gene)		
	def get_domtblout(self, gene, src='g'):
		outdir = self.hmmoutdir
		return '{}/{}.{}.domtblout'.format(outdir, gene, src)
	def get_filename(self, _dir, gene, *field):
		return '{}/{}.{}'.format(_dir, gene, '.'.join(field))
	def get_domtblfa(self, gene, **kargs):
		return self.get_domtblout(gene, **kargs) + '.fa'

	def hmmsearch(self, hmmfile, seqdb, domtblout):
		cmd = 'hmmsearch --domtblout {domtblout} {hmmfile} {seqdb} > /dev/null'.format(
				hmmfile=hmmfile, seqdb=seqdb, domtblout=domtblout)
		run_cmd(cmd, log=True)
		return cmd
	def double_seqs(self, seqfile, fout, seqfmt='fasta'):
		d_length = {}
		for rc in SeqIO.parse(seqfile, seqfmt):
			for seq, suffix0 in zip([rc.seq, rc.seq.reverse_complement()], ['fwd', 'rev']):
				suffix = '|{}'.format(suffix0)
				print >> fout, '>{}{}\n{}'.format(rc.id, suffix, seq)
			d_length[rc.id] = len(rc.seq)
		return d_length

	def _check_database(self):
		pass
	def _check_dependencies(self):
		pass

def main():
	args = makeArgparse()
	print >>sys.stderr, args.__dict__
	pipeline = Pipeline(**args.__dict__)
	pipeline.run()

def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("genome", action="store",type=str,
					help="input genome sequence in fasta format [required]")
	parser.add_argument("-prefix", action="store",
					default=None, type=str,
					help="output prefix [default='genome']")
	parser.add_argument("-est", action="store",type=str,
                    help="EST sequences for evidence")
	parser.add_argument('-organ', type=str, choices=['mt', 'pt'], default='mt',
						help="mt (mitochondrion) or pt (plastid) [default=%(default)s]")
	parser.add_argument('-taxon', type=str, default=None,
						help="taxon to filter out, such as Embryophyta [default=%(default)s]")
	
	args = parser.parse_args()
	return args

if __name__ == '__main__':
	main()

