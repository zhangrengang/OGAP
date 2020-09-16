import sys, os
import re
import glob
import argparse
from RunCmdsMP import run_cmd, logger
def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("genome", action="store",type=str,
		help="input genome sequence in fastaformat [required]")
	addReditArgs(parser)
	args = parser.parse_args()
	return args

def addReditArgs(parser):
	parser.add_argument('-rna_bam', type=str, default=None,
			help="sorted RNA-Seq BAM file for REDItools (-i) (required for RNA editing analysis) [default='%(default)s']")
	parser.add_argument('-dna_bam', type=str, default=None,
			help="sorted DNA-Seq BAM file for REDItools (-j) (will switch to REDItoolDnaRna) [default='%(default)s']")
	parser.add_argument('-prefix', type=str, default='',
			help="output prefix of REDItools [default='%(default)s']")
	parser.add_argument('-dnarna_opts', type=str, 
			default='-c 10,10 -q 20,20 -m 30,30 -g 1 -u -a 6-0 -v 3 -n0.0 -N0.0 -V',
			help="options for REDItoolDnaRna when both -rna_bam and -dna_bam are specified [default='%(default)s']")
	parser.add_argument('-denovo_opts', type=str, 
			default='',
			help="options for REDItoolDenovo when only -rna_bam is specified [default='%(default)s']")
	parser.add_argument('-select_opts', type=str,
			default='-d 12 -c 10 -C 10 -v 5 -V 0 -f 0.1 -F 1.0 -e -u',
			help="options for selectPositions [default='%(default)s']")

class ReditPipeline():
	def __init__(self, genome, 
				rna_bam, dna_bam=None,
				tmpdir='/tmp/',
				dnarna_opts='',
				denovo_opts='',
				select_opts='',
				prefix='redit_out',
				gff=None,
				):
		self.genome = genome
		self.rna_bam = rna_bam
		self.dna_bam = dna_bam
		self.tmpdir = tmpdir
		self.dnarna_opts = dnarna_opts
		self.denovo_opts = denovo_opts
		self.select_opts = select_opts
		self.prefix = prefix
		self.out_folder = '{}/{}'.format(tmpdir, prefix)
		self.gff = gff
		if self.dna_bam is None:
			self.tools = 'REDItoolDenovo.py'
			self.options = '{}'.format(self.denovo_opts)
		else:
			self.tools = 'REDItoolDnaRna.py'
			self.options = '-j {} {}'.format(self.dna_bam, self.dnarna_opts)
			self.select_opts += ' -d 12'
		self.options += ' -o {}'.format(self.out_folder)
		if gff is not None:
			
		self.main_cmd = '{} -i {} -f {} {}'.format(self.tools,
				self.rna_bam, self.genome, self.options)
		
	def run(self):
		rmdirs(self.out_folder)
		# REDItools
		out, err, status = run_cmd(self.main_cmd, log=True)
		out_table = self.get_output(err)
		# select
		cmd = 'selectPositions.py -i {} -o {} {}'.format(out_table, )
	def get_output(self, err):
		outfile = re.compile(r'Results saved on (.+)').search(err).groups()[0]
		return outfile
