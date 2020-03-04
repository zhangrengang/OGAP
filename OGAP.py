import sys, os
import argparse
from Bio import SeqIO
from lib.Database import Database
from lib.Hmmer import HmmSearch
from Gff import ExonerateGffGenes, AugustusGtfGenes
from RunCmdsMP import run_cmd, run_job, logger
from translate_seq import six_frame_translate
from small_tools import mkdirs, rmdirs

class Pipeline():
	def __init__(self, genome, 
				organ, taxon, 
				est=None, # EST evidence
				prefix=None, tmpdir='/dev/shm/tmp', 
				linear=False, 
				include_orf=False,
				max_evalue=1e-5, 
				seqfmt='fasta', **kargs):
		self.genome = os.path.realpath(genome)
		self.db = Database(organ=organ, taxon=taxon, include_orf=include_orf)
		self.ogtype = self.db.ogtype
		self.est = est
		if prefix is None:
			self.prefix = os.path.basename(genome)
		else:
			self.prefix = prefix
		self.tmpdir = '{}/{}'.format(tmpdir, self.ogtype)
		self.seqfmt = seqfmt
		self.hmmoutdir = '{}/hmmout'.format(self.tmpdir)
		self.estoutdir = '{}/estout'.format(self.tmpdir)
		self.exnoutdir = '{}/exnout'.format(self.tmpdir)
		self.agtoutdir = '{}/augustus'.format(self.tmpdir)

	def run(self):
		mkdirs(self.tmpdir)
		mkdirs(self.hmmoutdir, self.estoutdir, self.agtoutdir)
		# check
		logger.info('checking database: {}'.format(self.db.ogtype))
		self.db.checkdb()
		# est hmmsearch for evidence
		if self.est is not None:
			logger.info('HMMsearch EST')
			self.hmmsearch_est()
		# hints
		logger.info('prepare hints from protein [and EST] sequences')
		self.prepare_hints()
		# cds-protein hmmsearch
		logger.info('HMMsearch protein-encoding genes')
#		self.hmmsearch_protein()
		# cds augustus
		logger.info('using AUGUSTUS to predict protein-encoding genes')
		self.augustus_predict()
		# cds integrate
		
		# rna hmmsearch
		
		# tRNAScan

		# rna integrate
	def prepare_hints(self):
		for gene in self.db.cds_genes:
			seqfile = self.db.get_seqfile(gene)
			exn_gff = self.get_exnfile(gene, 'p')
			self.exonerate(seqfile, self.genome, exn_gff,		# protein
					model='protein2genome', percent=20,
					maxintron=500000,
					showtargetgff='T')
			hintfile = self.get_hintfile(gene)
			with open(hintfile, 'w') as fout:
				ExonerateGffGenes(exn_gff).to_hints(fout, src='P', pri=4)
				if self.est is not None:
					exn_gff = self.get_exnfile(gene, 'e')
					ExonerateGffGenes(exn_gff).to_hints(fout, src='E', pri=4)
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
			exn_gff = self.get_exnfile(gene, 'e')
			self.exonerate(est_seq, self.genome, exn_gff, 		# est
					model='est2genome', bestn=5, percent=70, 
					maxintron=500000,
					geneticcode=self.db.transl_table,
					showtargetgff='T')
		os.remove(aa_seq)
	def hmmsearch_protein(self):
		# translate
		aa_seq = '{}/{}.genome.aa'.format(self.tmpdir, self.prefix)
		with open(aa_seq, 'w') as fout:
			d_length = six_frame_translate(self.genome, fout, transl_table=self.db.transl_table)
		# hmmsearch
		for gene in self.db.cds_genes:
			hmmfile = self.db.get_hmmfile(gene)
			domtblout = self.get_domtblout(gene, src='g')
			self.hmmsearch(hmmfile, aa_seq, domtblout)
	def augustus_predict(self):
		for gene in self.db.cds_genes:
			augusuts_species = self.db.get_augustus_species(gene)
			pflfile = self.db.get_pflfile(gene)
			hintfile = self.get_hintfile(gene)
			outgff = self.get_augustus_gff(gene)
			kargs = {'translation_table': self.db.transl_table,
				'hintsfile': hintfile, 'extrinsicCfgFile': 'extrinsic.MPE.cfg',
				'proteinprofile': proteinprofile, '/ExonModel/minexonlength': 20,
				'codingseq': 1,
				}
			self.augustus(self.genome, augusuts_species, outgff, kargs=kargs)
			print >>sys.stderr, gene, self.check_augustus_gff(outgff)
	def check_augustus_gff(self, gff):
		genes, has_block, has_support, has_obey = 0,0,0,0
		full_support = 0
		both_support = 0
		for record in AugustusGtfGenes(gff):
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
				cmd += ['--{key} {value}'.format(key=key, value=value)]
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
	def get_domtblfa(self, gene, **kargs):
		return self.get_domtblout(gene, **kargs) + '.fa'

	def hmmsearch(self, hmmfile, seqdb, domtblout):
		cmd = 'hmmsearch --domtblout {domtblout} {hmmfile} {seqdb} > /dev/null'.format(
				hmmfile=hmmfile, seqdb=seqdb, domtblout=domtblout)
		run_cmd(cmd)
		return cmd
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

