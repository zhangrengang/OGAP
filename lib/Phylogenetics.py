import sys, os
import re
import glob
from Bio import SeqIO
from RunCmdsMP import run_cmd, logger
from OrthoFinder import catAln

class PhyloPipeline:
	def __init__(self, indir, 
			outprefix=None,
			tmpdir='/tmp',
			types=['cds'], 
			min_shared=50):
		self.indir = indir
		self.tmpdir = tmpdir
		self.types = types
		self.min_shared = min_shared
		self.outprefix = outprefix
		if outprefix is None:
			self.outprefix = '{}-{}'.format(
				os.path.basename(indir.strip('/')), '-'.join(types))
		
	def run(self):
		# bin
		self.genes = self.bin_seqs()
		# align
		alnfiles = self.align_seqs()
		# concat
		concat_alignment = '{}.aln'.format(self.outprefix)
		with open(concat_alignment, 'w') as fout:
			catAln(alnfiles, fout)
		# trim
		trimed_alignment = concat_alignment + '.trimal'
		self.trimal(concat_alignment, trimed_alignment)
		# tree
		self.iqtree(trimed_alignment)

	def iqtree(self, alnfile, opts='-nt AUTO'):
		cmd = 'iqtree -s {} -bb 1000 {} > /dev/null'.format(alnfile, opts)
		run_cmd(cmd, log=True)

	def trimal(self, inaln, outaln):
		cmd = 'trimal -in {} -automated1 > {}'.format(inaln, outaln)
		run_cmd(cmd)

	def align_seqs(self):
		alnfiles = []
		for gene in self.genes:
			seqfile = self.get_seqfile(gene)
			alnfile = self.get_alnfile(gene)
			cmd = 'mafft --auto {seq} > {aln}'.format(seq=seqfile, aln=alnfile)
			run_cmd(cmd)
			alnfiles += [alnfile]
		return alnfiles

	def bin_seqs(self):
		fastas = []
		prefixs = []
		for type in self.types:
			suffix = '.{}.fasta'.format(type)
			pattern = r'(\S+){}'.format(suffix)
			fastafiles = glob.glob('{}/*{}'.format(self.indir, suffix))
			prefixs += [self.get_prefix(fasta, pattern) for fasta in fastafiles]
			fastas += fastafiles
		d_seqs = {}
		for fasta, prefix in zip(fastas, prefixs):
			for rc in SeqIO.parse(fasta, 'fasta'):
				if re.compile('\-\d+$').search(rc.id):
					continue
				gene = rc.id
				rc.id = prefix
				try: d_seqs[gene] += [rc]
				except KeyError: d_seqs[gene] = [rc]

		nspecies = len(set(prefixs))
		# write
		for gene, records in d_seqs.items():
			shared = 1e2*len(records) / nspecies
			if shared < self.min_shared:
				logger.info('gene {} is only shared by {}% species, removed.\
							'.format(gene, shared))
				continue
			seqfile = self.get_seqfile(gene)
			with open(seqfile, 'w') as fout:
				for record in records:
					SeqIO.write(record, fout, 'fasta')
		return d_seqs.keys()

	def get_alnfile(self, gene):
		return '{}/{}.aln'.format(self.tmpdir, gene)
	def get_seqfile(self, gene):
		return '{}/{}.fa'.format(self.tmpdir, gene)
	
	def get_prefix(self, filename, pattern):
		filename = os.path.basename(filename)
		return re.compile(pattern).search(filename).groups()[0]

def main():
	indir = sys.argv[1]
	PhyloPipeline(indir).run()

if __name__ == '__main__':
	main()
