# coding: utf-8
import sys
import itertools
#from CAI import RSCU
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Data.IUPACData import protein_letters_1to3
from Bio.SeqUtils import CodonUsage as CU

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.table

class CodonUsage(CU.CodonAdaptationIndex):
	def __init__(self, cdsfile, transl_table=1):
		self.cdsfile = cdsfile
		
	def __iter__(self):
		return self.parse()
	def parse(self):
		self._count_codons(self.cdsfile)
		total_codon = sum(self.codon_count.values())
		for _aa in SynonymousCodons():
			aa, codons = _aa.aa, _aa.codons
			rcsu = _aa
			total = sum([self.codon_count[codon] for codon in codons])
			denominator = float(total) / len(codons)
			rcsu.conut = [self.codon_count[codon] for codon in codons]
			rcsu.rcsu = [self.codon_count[codon] / denominator for codon in codons]
			rcsu.percent = [1e2*self.codon_count[codon] / total_codon for codon in codons]
			yield rcsu
def main(cdsfile, outable, transl_table=1):
	line = ['AA', 'Codons', 'Count', 'Percent', 'RSCU', 'sum_Count', 'sum_Percent']
	print >> outable, '\t'.join(line)
	for rcsu in CodonUsage(cdsfile, transl_table=transl_table):
		aa = '{} ({})'.format(rcsu.aa3, rcsu.aa)
		codons = '/'.join(rcsu.codons)
		count = '/'.join(map(str, rcsu.conut))
		percent = '/'.join(map(lambda x:'{:.2f}%'.format(x), rcsu.percent))
		_rcsu = '/'.join(map(lambda x:'{:.2f}'.format(x), rcsu.rcsu))
		sum_count = sum(rcsu.conut)
		sum_pervent = '{:.2f}%'.format(sum(rcsu.percent))
		line = [aa, codons, count, percent, _rcsu, sum_count, sum_pervent]
		line = map(str, line)
		print >> outable, '\t'.join(line)

def main_plot(cdsfiles, transl_table=1):
	taxa = [list(CodonUsage(cdsfile, transl_table=transl_table)) 
				for cdsfile in cdsfiles]
	bar_plot(taxa)
colors_lst = ['#980000','#00ffff','#4a86e8','#ff9900', #'#ffff00', '#0000ff',
    '#9900ff','#ff00ff','#274e13','#000000','#cccccc','#7f6000',
    '#a64d79','#6aa84f','#fff2cc','#47a952','#3ea6b6','#a5b805','#8f9276','#ca8d7c',
    '#ff0000', '#00ff00','#20124d',]

def bar_plot(taxa, ylabel='RSCU'):
	iwidth = 0.7
	clist = colors_lst
	aa3s = [rcsu.aa3 for rcsu in taxa[0]]
	xl = aa3s
	xind = np.arange(len(xl))
	width = iwidth/len(taxa)
	fig=plt.figure(figsize=(10,6))
	ax=plt.subplot(111)
	plt.tick_params(top= 'off', right= 'off',bottom= 'off', left= 'on')
	ax.set_xticks(xind)
	ax.set_xticklabels('')
	ax.set_xlim(-0.5,len(xind)-0.5)
	plt.ylabel(ylabel)

	vx, vy = 0, 0
	# bar plots
	step = 0
	for i, rcsus in enumerate(taxa):	# per taxon
		for j, rcsu in enumerate(rcsus):			# per AA
			lasty = 0
			x = j - iwidth/2 + step
			for k, value in enumerate(rcsu.rcsu):	# per codon
				ax.bar(x, value, width, facecolor=clist[k], 
					edgecolor='w', linewidth=0, align='edge', bottom=lasty)
				lasty += value
			vy = max(vy, k+1)
		vx = max(vx, j+1)
		step += width*1.1

	# table text
	print vx, vy
	xx = vx * vy
	atable = np.array([''] *xx, dtype='|S10').reshape(vy, vx)
	acolor = np.array(['w']*xx, dtype='|S10').reshape(vy, vx)
	for j, rcsu in enumerate(rcsus):
		for k, codon in enumerate(rcsu.codons):
			atable[k, j] = codon
			acolor[k, j] = clist[k]
	print atable
	print acolor
	fig.linewidth=0
	the_table = ax.table(cellText=atable,cellColours=acolor,rowLabels=None,colLabels=xl,loc='bottom',
					cellLoc='center', edges='closed')
	for key, cell in the_table.get_celld().items():
		cell.set_linewidth(0.0)
	plt.savefig("CodonUsage.pdf",bbox_inches='tight',format="pdf")
class SynonymousCodons:
	def __init__(self, transl_table=1):
		try: 
			transl_table = int(transl_table)
			table = CodonTable.unambiguous_dna_by_id[transl_table]
		except ValueError: 
			table = CodonTable.unambiguous_dna_by_name[transl_table]
		self.table = table.forward_table
	def __iter__(self):
		return self.parse_table()
	def parse_table(self):
		aa_dict = {}
		for codon, aa in self.table.items():
			try: aa_dict[aa] += [codon]
			except KeyError: aa_dict[aa] = [codon]
		for aa, codons in sorted(aa_dict.items(), key=lambda x:protein_letters_1to3[x[0]]):
			yield AA(aa, codons)

class AA:
	def __init__(self, aa, codons=None):
		self.aa = aa
		self.codons = codons
	@property
	def aa3(self):
		return protein_letters_1to3[self.aa]

def calculate_RSCU(cdsfile, outable):
	seqs = [str(rc.seq) for rc in SeqIO.parse(cdsfile, 'fasta')]
	print(RSCU(seqs))

if __name__ == '__main__':
	#cdsfile = sys.argv[1]
	cdsfiles = sys.argv[1:]
	outable = sys.stdout
	#calculate_RSCU(cdsfile, outable)
	#main(cdsfile, outable)
	main_plot(cdsfiles)
