import sys, os
import re
import glob
import itertools
from Bio import SeqIO
from RunCmdsMP import run_cmd, logger
from OrthoFinder import catAln
from small_tools import mkdirs

class PhyloPipeline(object):
	def __init__(self, indir, 
			outprefix=None,
			gnid = True,
			tmpdir='./tmp',
			types=['cds'], 
			min_shared=50):
		self.indir = indir
		self.tmpdir = tmpdir
		self.types = types
		self.min_shared = min_shared
		self.outprefix = outprefix
		self.gnid = gnid
		if outprefix is None:
			self.outprefix = '{}-{}'.format(
				os.path.basename(indir.strip('/')), '-'.join(types))
		
	def run(self):
		mkdirs(self.tmpdir)
		self.d_gnid = self.get_genome_id()
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
			fastafiles, prefies = self.get_fastas(type)
			prefixs += prefies
			fastas += fastafiles
		d_seqs = {}
		
		for rc, prefix in self.get_seqs(fastas, prefixs):
			gene = rc.id
			rc.id = prefix
			if self.gnid:
				rc.id += '-{}'.format(self.d_gnid[prefix])
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
	
	def get_prefix(self, filename, suffix):
		pattern = r'(\S+){}'.format(suffix)
		filename = os.path.basename(filename)
		return re.compile(pattern).search(filename).groups()[0]
	def get_fastas(self, type):
		suffix = '.{}.fasta'.format(type)
		fastafiles = sorted(glob.glob('{}/*{}'.format(self.indir, suffix)))
		prefixs = [self.get_prefix(fasta, suffix) for fasta in fastafiles]
		return fastafiles, prefixs
	def get_genome_id(self):
		suffix = '.gb'
		gbfiles = sorted(glob.glob('{}/*{}'.format(self.indir, suffix)))
		prefixs = [self.get_prefix(fasta, suffix) for fasta in gbfiles]
		d = {}
		for prefix, gbfile in zip(prefixs, gbfiles):
			for rc in SeqIO.parse(gbfile, 'genbank'):
				break
			d[prefix] = rc.id
		return d
	def get_seqs(self, fastas, prefixs):
		for fasta, prefix in zip(fastas, prefixs):
			for rc in SeqIO.parse(fasta, 'fasta'):
				if re.compile('\-\d+$').search(rc.id):
					continue
				yield rc, prefix
class KaKsPipeline(PhyloPipeline):
	def __init__(self, indir, tmpdir='/dev/shm/tmp/', ncpu=20):
		self.indir = indir
		self.tmpdir = tmpdir
		self.ncpu = ncpu
		self.types = ['cds', 'pep']
		self.output = '{}.kaks.xls'.format(self.indir.rstrip('/'))
	def run(self):
		# cat
		files = self.cat_seqs()
		# run
		self.run_ParaAT(files)
	def run_ParaAT(self, files):
		outdir = '{}/paraat'.format(self.tmpdir)
		cmd = 'mkdir -p {outdir} && echo {ncpu} > {outdir}/../proc && \
	ParaAT.pl -h {hom} -n {cds} -a {pep} -o {outdir} -p {outdir}/../proc -f axt -kaks -nogap -m mafft '.format(
			ncpu=self.ncpu, outdir=outdir, output=self.output, **files)
		run_cmd(cmd, log=True)
		kaksfiles = sorted(glob.glob('{}/*.kaks'.format(outdir)))
		with open(self.output, 'w') as fout:
			self.merge_kaks(kaksfiles, fout)

	def merge_kaks(self, kaksfiles, fout):
		for i, kaksfile in enumerate(kaksfiles):
			for j, line in enumerate(open(kaksfile)):
				if i>0 and j==0:
					continue
				fout.write(line)
	def cat_seqs(self):
		d_genes = {}
		d_files = {}
		for type in self.types:
			catfasta = '{}/{}.fasta'.format(self.tmpdir, type)
			d_files[type] = catfasta
			fastas, prefixs = self.get_fastas(type)
			fout = open(catfasta, 'w')
			for rc, prefix in self.get_seqs(fastas, prefixs):
				gene = rc.id
				rc.id = '{}|{}'.format(prefix, gene)
				try: d_genes[gene] += [rc.id]
				except KeyError: d_genes[gene] = [rc.id]
				SeqIO.write(rc, fout, 'fasta')
			fout.close()
		homfile = '{}/homology'.format(self.tmpdir)
		d_files['hom'] = homfile
		fout = open(homfile, 'w')
		for gene, ids in sorted(d_genes.items()):
			ids = sorted(set(ids))
			for g1, g2 in itertools.combinations(ids, 2):
				print >>fout, '{}\t{}'.format(g1, g2)
		fout.close()
		return d_files

def main():
	subcmd = sys.argv[1]
	indir = sys.argv[2]
	if subcmd == 'phylo':
		PhyloPipeline(indir).run()
	elif subcmd == 'kaks':
		KaKsPipeline(indir).run()
	else:
		raise ValueError

if __name__ == '__main__':
	main()
