#coding: utf-8
import sys, os
import glob
import itertools
from Bio import SeqIO
from Bio import Phylo
try: from xopen import xopen as open
except ImportError: from small_tools import open_file as open
from collections import Counter
from RunCmdsMP import run_cmd, run_job, logger

def catAln(inALNs, outALN, allow_missing=True):
	'''首尾相连alignments'''
	if allow_missing:
		species = set([])
		for inALN in inALNs:
			for rc in SeqIO.parse(inALN, 'fasta'):
				sp, g = gene_format_common(rc.id)
				species.add(sp)
	names = [os.path.basename(aln).split('.')[0] for aln in inALNs]
	d_seqs = {}
	lens = []
	for inALN in inALNs:
		ntax = 0
		sps = set([])
		for rc in SeqIO.parse(inALN, 'fasta'):
			ntax += 1
			sp, g = gene_format_common(rc.id)
			seq = str(rc.seq)
			try: d_seqs[sp] += [seq]
			except KeyError: d_seqs[sp] = [seq]
			sps.add(sp)
		lens += [len(seq)]
		if allow_missing:
			for sp in species-sps:
				ntax += 1
				seq = ''.join(['-']*len(seq))
				try: d_seqs[sp] += [seq]
				except KeyError: d_seqs[sp] = [seq]
	xlens = ','.join(map(str, lens))
	names = ','.join(names)
	description = 'taxa:{} genes:{} sites:{} blocks:{} names:{}'.format(ntax, len(lens), sum(lens), xlens, names)
	for sp, seqs in d_seqs.items():
		seqs = ''.join(seqs)
		print >> outALN, '>{} {}\n{}'.format(sp, description, seqs)

class Group():  # 解析groups.tsv，迭代返回每行
	def __init__(self, inGrp):
		self.inGrp = inGrp
	def __iter__(self):
		return self.parse()
	def parse(self):
		i = 0
		for line in open(self.inGrp):
			i += 1
			if i == 1:
				self.species = line.strip().split('\t')
				continue
			yield GroupRecord(line, self.species)
class GroupRecord(object): # 解析每行
	def __init__(self, line, species):
		line = line.strip('\n\r').split('\t')
		self.ogid = self.id = line[0]
		self.genes = [genes.split(', ') for genes in line[1:]]
		self.species = species
		self.genes = [self._strip(genes) for genes in self.genes]
		self.raw_genes = line[1:]
	
	def get_group(self):
		for genes in self.genes:
			for gene in genes:
				yield gene
	@property
	def counts(self):
		return [len(genes) for genes in self.genes]
	@property
	def spdict(self):
		return dict(zip(self.species, self.genes))
	@property
	def counter(self):
		return {sp: len(genes) for sp, genes in sorted(self.spdict.items())}
#		return dict(zip(self.species, self.counts))
	@property
	def singlecopy_ratio(self):
		singles = [v for v in self.counts if v == 1]
		return 1.0*len(singles)/len(self.counts)
	@property
	def singlecopy_dict(self):
		return {genes[0]: sp for sp, genes, count in zip(self.species, self.genes, self.counts) if count==1}
	def _strip(self, values):
		return [v for v in values if v]

class OrthoMCLGroup():  # 解析groups.txt，迭代返回每行
	def __init__(self, inGrp):
		self.inGrp = inGrp
	def __iter__(self):
		return self.parse()
	def parse(self):
		for line in open(self.inGrp):
			yield OrthoMCLGroupRecord(line)
class OrthoMCLGroupRecord(GroupRecord): # 解析每行
	def __init__(self, line=None, genes=None):
		if line is not None:
			line = line.strip().split()
			self.ogid = line[0].strip(':')
			self.genes = line[1:]
		if genes is not None:
			self.genes = genes
		self.species = [gene.split('|')[0] for gene in self.genes]
	@property
	def spdict(self):
		d = {}
		for sp, gene in zip(self.species, self.genes):
			try: d[sp] += [gene]
			except KeyError:  d[sp] = [gene]
		return d

def to_hybpiper(ResultsDir, cdsSeq=None, outOGSeq=None, species=None, min_singlecopy=0.7, only_stats=False):
	if species:
		if isinstance(species, str):
			species = [line.strip().split()[0] for line in open(species)]
	Orthogroups = OrthoFinder(ResultsDir).Orthogroups
	if not only_stats:
		d_seqs = seq2dict(cdsSeq)
	ratios = []
	i = 0
	for og in Group(Orthogroups):
		if species:
			og.genes = [og.spdict[sp] for sp in species]
			og.species = species
		if only_stats:
			ratios += [og.singlecopy_ratio]
			continue
		if not og.singlecopy_ratio >= min_singlecopy:
			continue
		i += 1
		for gene, sp in og.singlecopy_dict.items():
			rc = d_seqs[gene]
			rc.id = '{}-{}'.format(sp.replace('-', '_'), og.id)
			SeqIO.write(rc, outOGSeq, 'fasta')
	if only_stats:
		print >>sys.stderr, '{}\t{}'.format('total', len(ratios))
		for cutoff in range(50, 105, 5):
			cutoff = cutoff/1e2
			ratios = filter(lambda x: x >= cutoff, ratios)
			print >>sys.stderr, '>={}\t{}'.format(cutoff, len(ratios))
		return
	print >>sys.stderr, '{} OGs'.format(i)

def to_astral(ResultsDir, pepSeq, outTrees, tmpdir='/io/tmp/share', min_singlecopy=0.7):
	from RunCmdsMP import run_job
	tmpdir = '{}/to_astral.{}'.format(tmpdir, os.getpid())
	logger.info('change tmpdir to {}'.format(tmpdir))
	if not os.path.exists(tmpdir):
		os.mkdir(tmpdir)
	Orthogroups = OrthoFinder(ResultsDir).Orthogroups
	d_seqs = seq2dict(pepSeq)
	cmd_list = []
	treefiles = []
	for og in Group(Orthogroups):
		if not og.singlecopy_ratio >= min_singlecopy:
			continue
		d_singlecopy = {genes[0]: sp for sp, genes, count in zip(og.species, og.genes, og.counts) if count==1}
		outSeq = '{}/{}.pep'.format(tmpdir, og.ogid)
		f = open(outSeq, 'w')
		for gene, sp in d_singlecopy.items():
			rc = d_seqs[gene]
			rc.id = sp
			SeqIO.write(rc, f, 'fasta')
		f.close()
		alnSeq = outSeq + '.aln'
		alnTrim = alnSeq + '.trimal'
		treefile = alnTrim + '.treefile'
		cmd = '[ ! -s {} ]'.format(treefile)
		cmds = [cmd]
		cmd = 'mafft --auto {} > {} 2> /dev/null'.format(outSeq, alnSeq)
		cmds += [cmd]
		cmd = 'trimal -automated1 -in {} -out {} &> /dev/null'.format(alnSeq, alnTrim)	# -gt 0.8 (before 2020-2-21)
		cmds += [cmd]
		cmd = 'iqtree -s {} -bb 1000 -mset JTT -nt 1 &> /dev/null'.format(alnTrim)
		cmds += [cmd]
		cmds = ' && '.join(cmds)
		cmd_list += [cmds]
		treefiles += [treefile]
	cmd_file = '{}/cmds.list'.format(tmpdir)
	run_job(cmd_file, cmd_list=cmd_list, tc_tasks=20)
	for treefile in treefiles:
		if not os.path.exists(treefile):
			logger.warn('{} do not exist, check the log file `{}`'.format(treefile, treefile.replace('.treefile', '.log')))
			continue
		for line in open(treefile):
			outTrees.write(line)
class OrthoFinder:
	def __init__(self, ResultsDir):
		self.ResultsDir = ResultsDir
#		self.SpeciesTreeAlignment = '{}/WorkingDirectory/Alignments_ids/SpeciesTreeAlignment.fa'.format(ResultsDir)
		self.SpeciesTreeAlignment = '{}/MultipleSequenceAlignments/SpeciesTreeAlignment.fa'.format(ResultsDir)
		self.SpeciesIDs = '{}/WorkingDirectory/SpeciesIDs.txt'.format(ResultsDir)
		self.SpeciesTree_rooted = '{}/Species_Tree/SpeciesTree_rooted_node_labels.txt'.format(ResultsDir)
		self.Orthologues = glob.glob('{}/Orthologues/*/*__v__*.tsv'.format(ResultsDir))
		self.Duplications = '{}/Gene_Duplication_Events/Duplications.tsv'.format(ResultsDir)
		self.Orthogroups = '{}/Orthogroups/Orthogroups.tsv'.format(ResultsDir)
		self.SequenceIDs = '{}/WorkingDirectory/SequenceIDs.txt'.format(ResultsDir)

	def orthogroups(self):
		return Group(self.Orthogroups)
	def get_orthogroups(self):
		for group in Group(self.Orthogroups):
			yield list(group.get_group())
	@property
	def OGDict(self):
		'''OG id -> gene id的字典'''
		d = {}
		for rc in Group(self.Orthogroups):
			for genes in rc.genes:
				for gene in genes:
					d[gene] = rc.ogid
		return d
	@property
	def Single_Copy_Orthologue(self):
		'''单拷贝OG'''
	#	Single_Copy_OGs = [fa.split('/')[-1].split('.')[0] for fa in \
	#			glob.glob('{}/Single_Copy_Orthologue_Sequences/*.fa'.format(self.ResultsDir))]
		Single_Copy_OGs = [line.strip() for line in \
			open('{}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt'.format(self.ResultsDir))]
		return Single_Copy_OGs
	def Single_Copy_Codon_Align(self, cdsSeqs, tmpdir='/tmp'):
		'''生成单拷贝OG的密码子比对'''
		Single_Copy_OGs = self.Single_Copy_Orthologue
		d_cds = seq2dict(cdsSeqs)
		ALNs = []
		for OG in Single_Copy_OGs:
			pepSeq = '{}/Single_Copy_Orthologue_Sequences/{}.fa'.format(self.ResultsDir, OG)
			cdsSeq = '{}/{}.cds'.format(tmpdir, OG)
			f = open(cdsSeq, 'w')
			for rc in SeqIO.parse(pepSeq, 'fasta'):
				SeqIO.write(d_cds[rc.id], f, 'fasta')
			f.close()
			pepAln = '{}/{}.pep.aln'.format(tmpdir, OG)
			cmd = 'mafft --auto {} > {} 2> /dev/null'.format(pepSeq, pepAln)
			os.system(cmd)
			pepTrim = pepAln + '.trimal'
			
			cdsAln = cdsSeq + '.aln'
			cmd = 'pal2nal.pl -output fasta {} {} > {} 2> /dev/null'.format(pepAln, cdsSeq, cdsAln)
			os.system(cmd)
			if os.path.getsize(cdsSeq) >0 and os.path.getsize(cdsAln) == 0:
				print >> sys.stderr, 'Error in CMDS `{}`'.format(cmd) 
				continue
			cdsTrim = cdsAln + '.trimal'
			cmd = 'trimal -gt 0.8 -in {} -out {} &> /dev/null'.format(cdsAln, cdsTrim)
			os.system(cmd)
			ALNs += [cdsTrim]
		return ALNs
	def Single_Copy_Pep_Align(self):
		return ['{}/MultipleSequenceAlignments/{}.fa'.format(self.ResultsDir, OG) \
					for OG in self.Single_Copy_Orthologue]
	@property
	def SequenceIDdict(self):
		'''蛋白序列新编id与原id的映射'''
		d = {}
		for line in open(self.SequenceIDs):
			if line.startswith('#'):
				continue
			id, geneName = line.strip().split(': ')
			d[id] = geneName
		return d
	def get_SequenceDict(self):
		'''蛋白序列id和record的映射关系'''
		self.Sequences = glob.glob('{}/WorkingDirectory/Species*.fa'.format(self.ResultsDir))
		d_seq_ids = self.SequenceIDdict
		d_seqs = {}
		for Sequence in self.Sequences:
			for rc in SeqIO.parse(Sequence, 'fasta'):
				rc.id = d_seq_ids[rc.id]
				d_seqs[rc.id] = rc
		return d_seqs
		
	@property
	def SpeciesIDdict(self):
		'''物种新编id和名称的映射关系'''
		d = {}
		for line in open(self.SpeciesIDs):
			if line.startswith('#'):
				continue
			temp = line.strip().split()
			id, spName = temp[:2]
			spName = '.'.join(spName.split('.')[:-1])
			id = id.strip(':')
			d[id] = spName
		return d
	@property
	def Species(self):
		'''物种名称的列表'''
		species = []
		for line in open(self.SpeciesIDs):
			if line.startswith('#'):
				continue
			temp = line.strip().split()
			id, spName = temp[:2]
			spName = '.'.join(spName.split('.')[:-1])
			species += [spName]
		return species
	@property
	def reverse_SpeciesIDdict(self):
		'''物种名称和新编id的字典'''
		return dict([(spName, id) for id, spName in self.SpeciesIDdict.items()])
	def spName2Id(self, *sps):
		'''获取指定物种的新编id'''
		d_sp = self.reverse_SpeciesIDdict
		return [d_sp[sp] for sp in sps]
	def spId2Name(self, *sps):
		'''获取指定物种id的名称'''
		d_sp = self.SpeciesIDdict
		return [d_sp[sp] for sp in sps]
	def get_species(self, genes, sep='|'):
		'''由基因名获取物种名，继承自OrthoMCL的编号格式'''
		return [gene.split(sep)[0] for gene in genes]
	def get_species_specific(self, sp):
		'''获取指定物种特有的所有基因'''
		i = 0
		d_genes = {}
		all_sp_genes = []
		for line in open(self.Orthogroups):
			i += 1
			temp = line.rstrip('\r\n').split('\t')
			if i == 1:
				spidx = temp.index(sp)
				continue
			od_id = temp[0]
			sp_genes = temp[spidx]
			temp.pop(spidx)
			genes = temp[1:]
#			print >> sys.stderr, genes
			sp_genes = sp_genes.split(', ')
			all_sp_genes += sp_genes
#			if sp_genes:
#				print >> sys.stderr,genes
			if not (sp_genes and set(genes) == set([''])):
				continue
			for gene in sp_genes:
				d_genes[gene] = od_id
		all_sp_genes = set(all_sp_genes)
		sp_id0 = self.reverse_SpeciesIDdict[sp]
		for seq_id, gene in self.SequenceIDdict.items():
			sp_id, sid = seq_id.split('_')
			if sp_id == sp_id0 and gene not in all_sp_genes:
				d_genes[gene] = None
		return d_genes
	def get_blast_files(self, sp1, sp2, byName=True):
		'''获取两物种间或物种内的blast文件'''
		if byName:
			sp1, sp2 = self.spName2Id(sp1, sp2)
		blast_files = []
		if sp1 == sp2:
			SPs = [(sp1, sp2)]
		else:
			SPs = [(sp1, sp2), (sp2, sp1)]
		for sp1, sp2 in SPs:
			blast_files += ['{}/WorkingDirectory/Blast{}_{}.txt.gz'.format(self.ResultsDir, sp1, sp2)]
		return blast_files
	def get_paralogs(self, sp=None, byName=True, min_support=0.5):
		'''由Gene_Duplication_Events/获取paralogs'''
		if sp is not None and not byName:
			sp, = self.spId2Name(sp)
		Duplications = self.Duplications
		para_pairs = set([])
		i = 0
		for line in open(Duplications):
			i += 1
			if i == 1:
				continue
			temp = line.strip().split('\t')
			Orthogroup, Species_Trer_Node, Gene_Tree_Node, Support, Type, Genes_1, Genes_2 = temp
			Support = float(Support)
			if Support < min_support:
				continue
			Genes_1 = Genes_1.split(', ')
			Genes_2 = Genes_2.split(', ')
			Genes_1 = map(gene_format_p, Genes_1)
			Genes_2 = map(gene_format_p, Genes_2)
			for (sp1, g1), (sp2, g2) in itertools.product(Genes_1, Genes_2):
				if sp1 == sp2:
					if sp is not None and sp != sp1:
						continue
					if (g2, g1) in para_pairs:
						continue
					para_pairs.add( (g1, g2) )
		return para_pairs
	def get_paralogs2(self, sp=None, byName=True):
		'''由Orthologues/获取paralogs'''
		if sp is not None:
			if not byName:
				sp, = self.spId2Name(sp)
			orthoFiles = glob.glob('{}/Orthologues/*{}*.tsv'.format(self.ResultsDir, sp))
		else:
			orthoFiles = self.Orthologues
		para_pairs = set([])
		idx = 0
		for orthoFile in orthoFiles:
			idx += 1
			i = 0
			for line in open(orthoFile):
				i += 1
				if i == 1:
					continue
				temp = line.strip().split('\t')
				Orthogroup, Genes_1, Genes_2 = temp
				Genes_1 = Genes_1.split(', ')
				Genes_2 = Genes_2.split(', ')
				Genes_1 = map(gene_format_o, Genes_1)
				Genes_2 = map(gene_format_o, Genes_2)
				for Genes in [Genes_1, Genes_2]:
					if len(Genes) < 2:
						continue
					for (sp1, g1), (sp2, g2) in itertools.combinations(Genes, 2):
						assert sp1 == sp2
						if sp is not None and sp != sp1:
							continue
						if (g2, g1) in para_pairs:
							continue
						para_pairs.add( (g1, g2) )
		return para_pairs
	def get_paralogs3(self):
		'''由Orthogroups/Orthogroups.tsv获取paralogs'''
		para_pairs = set([])
		for rc in Group(self.Orthogroups):
		#	print vars(rc)
			for genes in rc.genes:
				for g1, g2 in itertools.combinations(genes, 2):
					para_pairs.add( (g1, g2) )
		return para_pairs
	def get_orthologs(self, sps=None, sp1=None, sp2=None, byName=True):
		'''获取成对的orthologs'''
		if sps is not None:
			orthoFiles = []
			for sp1, sp2 in itertools.permutations(sps, 2):
				if not byName:
					sp1, sp2 = self.spId2Name(sp1, sp2)
				orthoFiles += ['{}/Orthologues/Orthologues_{}/{}__v__{}.tsv'.format(self.ResultsDir, sp1, sp1, sp2)]
		elif sp1 is not None or sp2 is not None:
			if not byName:
				sp1, sp2 = self.spId2Name(sp1, sp2)
			orthoFiles = ['{}/Orthologues/Orthologues_{}/{}__v__{}.tsv'.format(self.ResultsDir, sp1, sp1, sp2),
						  '{}/Orthologues/Orthologues_{}/{}__v__{}.tsv'.format(self.ResultsDir, sp2, sp2, sp1),]
		else:
			orthoFiles = self.Orthologues
		ortho_pairs = set([])
		idx = 0
		for orthoFile in orthoFiles:
			idx += 1
			i = 0
			for line in open(orthoFile):
				i += 1
				if i == 1:
					continue
				temp = line.strip().split('\t')
				Orthogroup, Genes_1, Genes_2 = temp
				Genes_1 = Genes_1.split(', ')
				Genes_2 = Genes_2.split(', ')
				Genes_1 = map(gene_format_o, Genes_1)
				Genes_2 = map(gene_format_o, Genes_2)
				for (sp1, g1), (sp2, g2) in itertools.product(Genes_1, Genes_2):
					assert sp1 != sp2
					if (g2, g1) in ortho_pairs:
						continue
					ortho_pairs.add( (g1, g2) )
		return ortho_pairs
	def get_orthologs_cluster(self, **kargs):
		'''获取orthologs的cluster，属OG的子集'''
		import networkx as nx
		G = nx.Graph()
		for g1, g2 in self.get_orthologs(**kargs):
			G.add_edge(g1, g2)
#		for g1, g2 in self.get_paralogs2():
#			G.add_edge(g1, g2)
		for cmpt in nx.connected_components(G):
			yield cmpt
		
	def get_singlecopy_orthologs(self, **kargs):
		for cmpt in self.get_orthologs_cluster(**kargs):
			if self.is_singlecopy(cmpt):
				yield cmpt
	def is_singlecopy(self, cmpt):	
		sps = [gene_format_o(g)[0] for g in cmpt]
		if len(sps) == len(set(sps)):
			return True
		else:
			return False
	def convert_seq_id(self, seqfile, out_seqfile,  fmt='fasta'):
		'''convert species ID to species NAME in sequence file'''
		d_species = self.SpeciesIDdict
		if out_seqfile == self.SpeciesTreeAlignment:
			raise ValueError('stop to write {}'.format(out_seqfile))
		f = open(out_seqfile, 'w')
		for rc in SeqIO.parse(seqfile, fmt):
			rc.id = d_species[rc.id]
			rc.description += ' {} sites'.format(len(rc.seq))
			SeqIO.write(rc, f, fmt)
	def get_root(self, treefile, fmt='newick'):
		'''从树文件获取root'''
		tree = Phylo.read(treefile, fmt)
		root  = tree.root
		if root.name == 'N0':
			for clade in root.clades:
				if not clade.name == 'N1':
					return clade.name
		return root.name
	def re_root(self, treefile, root, out_treefile, fmt='newick'):
		'''由Phylo进行reroot'''
		tree = Phylo.read(treefile, fmt)
		for clade in tree.find_clades(root):
			if root == clade.name:
				root = clade
				break
		else:
			print >>sys.stderr, 'root {} is not found'.format(root)
		tree.root_with_outgroup(root)
		Phylo.write(tree, out_treefile, fmt)
	def get_aln_len(self, alnfile, fmt='fasta'):
		'''获取alignment的长度'''
		for rc in SeqIO.parse(alnfile, fmt):
			return len(rc.seq)
		
def to_paml(OFdir, outDir, cdsSeq):
	'''生成PAML正选择分析所需文件；取共有基因；多拷贝基因取树枝最短的那个'''
	if not os.path.exists(outDir):
		os.mkdir(outDir)

	result = OrthoFinder(OFdir)
	species = result.Species
	d_seqs = result.get_SequenceDict()
	i, j, s = 0,0,0
	groups = []
	for genes in result.get_orthologs_cluster():
		i += 1
		g_species = result.get_species(genes)
		if not (set(species) & set(g_species) == set(species)):
			continue
		j += 1
		if len(species) == len(g_species):	# single copy
			s += 1
			groups += [genes]
			continue
#		if j >3:
#			break
		groups += [select_genes_bytree(genes, d_seqs, j)]
#	return
	print 'total {} groups, shared {}, single copy {}'.format(i, j, s)
	i = 0
	d_cds = seq2dict(cdsSeq)
	outGroup = '{}/groups.txt'.format(outDir)
	f = open(outGroup, 'w')
	for group in groups:
		i += 1
		og = 'OG_{}'.format(i)
		print >>f, '{}: {}'.format(og, ' '.join(group))
		outSeq = '{}/{}.pep'.format(outDir, og)
		outCds = '{}/{}.cds'.format(outDir, og)
		f1 = open(outSeq, 'w')
		f2 = open(outCds, 'w')
		for gene in group:
			rc = d_seqs[gene]
			SeqIO.write(rc, f1, 'fasta')
			rc = d_cds[gene]
			SeqIO.write(rc, f2, 'fasta')
		f1.close()
		f2.close()

		alnSeq = outSeq + '.aln'
		cmd = 'mafft --auto {} > {} 2> /dev/null'.format(outSeq, alnSeq)
		os.system(cmd)
		cdsSeq = outCds + '.aln'
		cmd = 'pal2nal.pl -output fasta {} {} > {} 2> /dev/null'.format(alnSeq, outCds, cdsSeq)
		os.system(cmd)
	f.close()
def seq2dict(seq):
	return dict([(rc.id, rc)for rc in SeqIO.parse(seq, 'fasta')])
def select_genes_bytree(genes, d_seqs, i):
	'''从树上选择多拷贝基因中的一个：树枝最短的那个'''
	outSeq = '/io/tmp/share/xxx{}.fa'.format(i)
	with open(outSeq, 'w') as f:
		for gene in genes:
			rc = d_seqs[gene]
			rc.id = rc.id.replace('|', '-')
			SeqIO.write(rc, f, 'fasta')
	alnSeq = outSeq + '.aln'
	print "if [ $SGE_TASK_ID -eq {} ]; then".format(i)
	cmd = 'mafft --auto {} > {} 2> /dev/null'.format(outSeq, alnSeq)
	print cmd
#	os.system(cmd)
	cmd = 'iqtree -s {} -pre {} -nt AUTO &> /dev/null'.format(alnSeq,alnSeq)
#	os.system(cmd)
	print cmd
	print 'fi'
#	return
	treefile = alnSeq + '.treefile'
	d_dist = {}
	tree = Phylo.read(treefile, 'newick')
	for gene1 in genes:
		sp1 = gene1.split('|')[0]
		gene0 = gene1
		gene1 = gene1.replace('|', '-')
		dist = 0
		for gene2 in genes:
			sp2 = gene2.split('|')[0]
			gene2 = gene2.replace('|', '-')
			if sp1 == sp2:
				continue
			dist += tree.distance(gene1,gene2)
		d_dist[gene0] = dist
	node0 = min(genes, key=lambda x:d_dist[x])
	sp0 = node0.split('|')[0]
	nodes = [node0]
	d_dist = {}
	node0 = node0.replace('|', '-')
	for gene1 in genes:
		sp1 = gene1.split('|')[0]
		gene0 = gene1
		if sp1 == sp0:
			continue
		gene1 = gene1.replace('|', '-')
	#	print gene1,node0
		dist = tree.distance(gene1,node0)
		try: d_dist[sp1][gene0] = dist
		except KeyError: d_dist[sp1] = {gene0: dist}
	for sp, d_sp_dist in d_dist.items():
		node = min(d_sp_dist.keys(), key=lambda x:d_sp_dist[x])
		nodes += [node]
	print sorted(genes),sorted(nodes)
	return nodes
def single_copy_cds_align(OFdir, cdsSeqs, outALN, tmpdir='./tmp'):
	'''首尾连接单拷贝的CDS alignment'''
	if not os.path.exists(tmpdir):
		os.mkdir(tmpdir)

	result = OrthoFinder(OFdir)
	inALNs = result.Single_Copy_Codon_Align(cdsSeqs, tmpdir=tmpdir)
	catAln(inALNs, outALN)
def single_copy_pep_align(OFdir, outALN, ):
	result = OrthoFinder(OFdir)
	inALNs = result.Single_Copy_Pep_Align()
	catAln(inALNs, outALN)
def cafe_count(OFdir, outCount):
	'''生成CAFE所需的基因家族计数文件'''
	def _count(genes):
		return len([v for v in genes.split(', ') if v.strip()])
	result = OrthoFinder(OFdir)
	ogfile = result.Orthogroups
	i = 0
	for line in open(ogfile):
		i += 1
		temp = line.rstrip('\r\n').split('\t')
		if i == 1:
			species = temp[1:]
			line = ['Desc', 'Family ID'] + species
			print >> outCount, '\t'.join(line)
			continue
		ogId = temp[0]
		group = temp[1:]
		count = [_count(genes) for genes in group]
		#if ogId == 'OG0000008':
		#	print ogId, group, count
		line = ['(null)', ogId] + count
		line = map(str, line)
		print >> outCount, '\t'.join(line)
def species_specific_genes(OFdir, sp, outTsv):
	'''物种特有基因'''
	result = OrthoFinder(OFdir)
	print >> outTsv, '\t'.join(['gene', 'group'])
	for gene, group in result.get_species_specific(sp).items():
		print >> outTsv, '\t'.join([gene.split('|')[1], str(group)])
	
def bootstrap_species_tree(OFdir, outdir, bootstrap=1000, iqtree_options='-mset JTT'):
	'''重新用iqtree建树'''
	result = OrthoFinder(OFdir)
	msa = result.SpeciesTreeAlignment
	treefile = result.SpeciesTree_rooted
	new_msa = '{}/singlecopy_aligned.faa'.format(outdir)
	#os.link(msa, new_msa)
	cmd = 'trimal -in {} -gt 0.8 > {}'.format(msa, new_msa)
	run_cmd(cmd)
	prefix = new_msa
#	result.convert_seq_id(msa, new_msa)
#	print >>sys.stderr, 'alignments have {} sites'.format(result.get_aln_len(new_msa))
	root = result.get_root(result.SpeciesTree_rooted)
	cmd = "iqtree -s {} -pre {} -bb {} -bnni -nt AUTO  {} -o {} > /dev/null".format(new_msa, prefix, bootstrap, iqtree_options, root)
	print >>sys.stderr, 'running cmd: {}'.format(cmd)
	os.system(cmd)
	new_treefile = '{}.treefile'.format(prefix)
	new_treefile_rooted = '{}.rooted.tre'.format(prefix)
	print >>sys.stderr, 're-root with {}'.format(root)
	#result.re_root(new_treefile, root, new_treefile_rooted)
	cmd = 'nw_reroot {} {} > {}'.format(new_treefile, root, new_treefile_rooted)
	os.system(cmd)
	return new_treefile_rooted
def singlecopy_tree(OFdir, outdir, bootstrap=1000, iqtree_options='-mset JTT'):
	'''取单拷贝基因建树'''
	result = OrthoFinder(OFdir)
	new_msa = '{}/singlecopy_aligned.faa'.format(outdir)
	prefix = new_msa
	root = result.get_root(result.SpeciesTree_rooted)
	with open(new_msa, 'w') as outALN:
		single_copy_pep_align(OFdir, outALN)
	cmd = "iqtree -s {} -pre {} -bb {} -bnni -nt AUTO  {} -o {} > /dev/null".format(new_msa, prefix, bootstrap, iqtree_options, root)
	os.system(cmd)
	new_treefile = '{}.treefile'.format(prefix)
	new_treefile_rooted = '{}.rooted.tre'.format(prefix)
	cmd = 'nw_reroot {} {} > {}'.format(new_treefile, root, new_treefile_rooted)
	os.system(cmd)
	return new_treefile_rooted
def get_singlecopy_orthologs(OFdir, outHomo, **kargs):
	result = OrthoFinder(OFdir)
	for genes in result.get_singlecopy_orthologs(**kargs):
		print >> outHomo, '\t'.join(genes)
def get_orthologs(OFdir, outHomo, **kargs):
	result = OrthoFinder(OFdir)
	for g1, g2 in result.get_orthologs(**kargs):
		line = [g1, g2]
		print >> outHomo, '\t'.join(line)
	return 
	orthoFiles = result.Orthologues
	idx = 0
	for orthoFile in orthoFiles:
		idx += 1
		i = 0
		for line in open(orthoFile):
			i += 1
			if i == 1:
				continue
			temp = line.strip().split('\t')
			Orthogroup, Genes_1, Genes_2 = temp
			Genes_1 = Genes_1.split(', ')
			Genes_2 = Genes_2.split(', ')
			Genes_1 = map(gene_format_o, Genes_1)
			Genes_2 = map(gene_format_o, Genes_2)
			info = '%s.%s.%s' % (Orthogroup, idx, i-1)
			for (sp1, g1), (sp2, g2) in itertools.product(Genes_1, Genes_2):
				if sp1 != sp2:
					line = [g1, g2, info]
					print >> outHomo, '\t'.join(line)
def gene_format_o(gene):
	sp, g = gene.split('|')
	return (sp, gene)
def get_paralogs(OFdir, outHomo, min_support=0.5):
	result = OrthoFinder(OFdir)
	Duplications = result.Duplications
	i = 0
	for line in open(Duplications):
		i += 1
		if i == 1:
			continue
		temp = line.strip().split('\t')
		Orthogroup, Species_Trer_Node, Gene_Tree_Node, Support, Type, Genes_1, Genes_2 = temp
		Support = float(Support)
		if Support < min_support:
			continue
		Genes_1 = Genes_1.split(', ')
		Genes_2 = Genes_2.split(', ')
		Genes_1 = map(gene_format_p, Genes_1)
		Genes_2 = map(gene_format_p, Genes_2)
		info = [Orthogroup, Species_Trer_Node, Gene_Tree_Node, Support, Type]
		info = map(str, info)
		info = '_'.join(info)
		for (sp1, g1), (sp2, g2) in itertools.product(Genes_1, Genes_2):
			if sp1 == sp2:
				line = [g1, g2, info]
				print >> outHomo, '\t'.join(line)
def gene_format_p(gene):
	sp, g = gene.split('|')
	sp = sp[:len(sp)/2]
	g = sp + '|' + g
	return (sp, g)
def gene_format_common(gene):
	if not '|' in gene:
		sp, g = gene, gene
		return (sp, g)
	sp, g = gene.split('|')
	sp1 = sp[:len(sp)/2]
	sp2 = sp[len(sp)/2+1:]
	if sp1 == sp2:
		sp = sp1
		g = sp + '|' + g
	return (sp, g)
def MCScanX_transposed(OFdir, tsp, cspp, spmap, gff, datadir='data', outdir='result',):
	'''MCScanX_transposed流程'''
	if not os.path.exists(datadir):
		os.mkdir(datadir)
	d_spmap = spmap2dict(spmap)
	result = OrthoFinder(OFdir)
	#prepare t.gff, t.blast
	t_abr = d_spmap[tsp]
	t_gff = '%s/%s.gff' % (datadir, t_abr)
	checkpoint = t_gff + '.ok'
	if not os.path.exists(checkpoint):
		prepare_gff(gff, t_gff, [tsp], d_spmap)
		os.mknod(checkpoint)
	else:
		print 'checkpoint: %s exists, skipping prepare %s' % (checkpoint, t_gff)

	t_blast = '%s/%s.blast' % (datadir, t_abr)
	d_geneIds = {}
	checkpoint = t_blast + '.ok'
	if not os.path.exists(checkpoint):
		d_geneIds = result.SequenceIDdict
		print d_geneIds.items()[:10]
		blast_files = result.get_blast_files(tsp, tsp)
		pairs = result.get_paralogs(tsp)
		outHomo = '%s/%s.homology' % (datadir, t_abr)
		write_homo(pairs, outHomo)
		prepare_blast(tsp, tsp, pairs, d_geneIds, blast_files, t_blast)
		os.mknod(checkpoint)
	else:
		print 'checkpoint: %s exists, skipping prepare %s' % (checkpoint, t_blast)

	# prepare c-t.gff, c-t.blast
	for csp in cspp:
		c_abr = d_spmap[csp]
		c_gff = '%s/%s_%s.gff' % (datadir, t_abr, c_abr)
		checkpoint = c_gff + '.ok'
		if not os.path.exists(checkpoint):
			prepare_gff(gff, c_gff, [tsp,csp], d_spmap)
			os.mknod(checkpoint)
		else:
			print 'checkpoint: %s exists, skipping prepare %s' % (checkpoint, c_gff)

		c_blast = '%s/%s_%s.blast' % (datadir, t_abr, c_abr)
		checkpoint = c_blast + '.ok'
		if not os.path.exists(checkpoint):
			if not d_geneIds:
				d_geneIds = result.SequenceIDdict
			blast_files = result.get_blast_files(tsp, csp)
			pairs = result.get_orthologs(tsp, csp)
			outHomo = '%s/%s_%s.homology' % (datadir, t_abr, c_abr)
			write_homo(pairs, outHomo)
			prepare_blast(tsp, csp, pairs, d_geneIds, blast_files, c_blast)
			os.mknod(checkpoint)
		else:
			print 'checkpoint: %s exists, skipping prepare %s' % (checkpoint, c_blast)
	# run
	c_abrs = ','.join([d_spmap[csp] for csp in cspp])
	suffix = '_'.join([d_spmap[sp] for sp in [tsp]+cspp])
	outdir += '.' + suffix
	log = 'run_%s.log' % (suffix,)
	cmd = 'MCScanX_h-transposed.pl -i %s -t %s -c %s -o %s -x %s &> %s' % (datadir, t_abr, c_abrs, outdir, len(cspp), log)
	print 'CMD: %s' % cmd
	os.system(cmd)
	outcount = 'run_%s.count.xls' % (suffix,)
	count_mcscan(t_abr, outdir, outcount)

def count_mcscan(t_abr, outdir, outcount):
	def _count(inFile):
		return sum([1 for line in open(inFile)]) - 1
	def _out_pairs(inFile, Type, outf):
		i = 0
		for line in open(inFile):
			i += 1
			if i == 1:
				continue
			temp = line.split()
			g1, g2 = temp[0], temp[2]
			line = [g1, g2, Type]
			print >> outf, '\t'.join(line)
	pairType = outdir + '.pair.class'
	f1 = open(pairType, 'w')
	f = open(outcount, 'w')
	line = ['mode', 'genes', 'pairs']
	print >>f, '\t'.join(line)
	for type in ['proximal', 'segmental', 'tandem', 'transposed']:
		genes = '%s/%s.%s.genes' % (outdir, t_abr, type)
		pairs = '%s/%s.%s.pairs' % (outdir, t_abr, type)
		_out_pairs(pairs, type, f1)
		gene_num = _count(genes)
		pair_num = _count(pairs)
		line = [type, gene_num, pair_num]
		line = map(str, line)
		print >>f, '\t'.join(line)
	for enope in ['after', 'between']:
		genesx, pairsx = glob.glob('{}/{}.transposed_{}_*.genes'.format(outdir, t_abr, enope)), \
						 glob.glob('{}/{}.transposed_{}_*.pairs'.format(outdir, t_abr, enope))
		for genes, pairs in zip(genesx, pairsx):
			type1 = genes.split('/')[-1].split('.')[-2]
			type2 = pairs.split('/')[-1].split('.')[-2]
			assert type1 == type2
			_out_pairs(pairs, type1, f1)
			gene_num = _count(genes)
			pair_num = _count(pairs)
			line = [type1, gene_num, pair_num]
			line = map(str, line)
			print >>f, '\t'.join(line)
	f1.close()
	f.close()
def write_homo(pairs, outHomo):
	f = open(outHomo, 'w')
	for line in pairs:
		print >>f, '\t'.join(line)
	f.close()
def prepare_blast(sp01, sp02, pairs, d_geneIds, blast_files, outblast):
	print 'extract blast of {} from {}'.format([sp01, sp02], blast_files)
	sp0x = sorted([sp01, sp02])
	gene_pair_set = pairs

	d_blast = {}
	for blasts in blast_files:
		for line in open(blasts):
			temp= line.strip().split('\t')
			g1, g2 = temp[:2]
			g1, g2 = d_geneIds[g1], d_geneIds[g2]
			if (g1, g2) in gene_pair_set or (g2, g1) in gene_pair_set:
				temp[:2] = g1, g2
				bscore = float(temp[11])
				if (g1, g2) in d_blast and d_blast[(g1, g2)][1] < bscore:
					d_blast[(g1, g2)] = [temp, bscore]
				elif (g1, g2) not in d_blast:
					d_blast[(g1, g2)] = [temp, bscore]
	f = open(outblast, 'w')
	for key, (temp, bscore) in d_blast.items():
		print >>f, '\t'.join(temp)
	f.close()

def prepare_gff(inGff, outGff, spp, d_spmap):
	print 'extract gff of {} from {}'.format(spp, inGff)
	spp = set(spp)
	f = open(outGff, 'w')
	for line in open(inGff):
		temp= line.strip().split()
		chr = temp[0]
		gene = temp[1]
		if d_spmap[chr] in spp:
			f.write(line)
	f.close()
def spmap2dict(spmap):
	d = {}
	for line in open(spmap):
		temp= line.strip().split()
		d[temp[0]] = temp[2]
	return d
def get_unique_logs(OFdir, outPrefix=''):
	result = OrthoFinder(OFdir)
	out_orth = '{}orthologs.txt'.format(outPrefix)
	out_para = '{}inparalogs.txt'.format(outPrefix)
	out_para2 = '{}inparalogs2.txt'.format(outPrefix)
	out_para3 = '{}inparalogs3.txt'.format(outPrefix)
	with open(out_orth, 'w') as f:
		for g1, g2 in result.get_orthologs():
			print >>f, '\t'.join([g1, g2])
	with open(out_para, 'w') as f:
		for g1, g2 in result.get_paralogs():
			print >>f, '\t'.join([g1, g2])
	with open(out_para2, 'w') as f:
		for g1, g2 in result.get_paralogs2():
			print >>f, '\t'.join([g1, g2])
	with open(out_para3, 'w') as f:
		for g1, g2 in result.get_paralogs3():
			print >>f, '\t'.join([g1, g2])
def aln2beast(inALN, outNEX, partify=True):
	import re
	i = 0
	ntax = sum([1 for rc in SeqIO.parse(inALN, 'fasta')])
	for rc in SeqIO.parse(inALN, 'fasta'):
		i += 1
		if i == 1:
			datatype = guess_seqtype(rc.seq)
			desc = rc.description
			try:
				genes = re.compile(r'genes:(\d+)').search(desc).groups()[0]
			except:
				genes = 0
			try:
				sites = re.compile(r'sites:(\d+)').search(desc).groups()[0]
			except:
				sites = len(rc.seq)
			try:
				blocks = re.compile(r'blocks:(.*)').search(desc).groups()[0]
				partitions = re.compile(r'(\d+)').findall(blocks)
				partitions = map(int, partitions)
			except:
				partitions = []
			assert len(partitions) == int(genes)
			print >> outNEX, '''#NEXUS
begin data;
dimensions ntax={ntax} nchar={nchar};
format datatype={datatype} interleave=no gap=-;
matrix'''.format(ntax=ntax, nchar=sites, datatype=datatype)
		print >> outNEX, '{id}\t{seq}'.format(id=rc.id, seq=rc.seq)
	print >> outNEX, ''';
end;'''
	if partify and partitions:
		print >> outNEX, 'begin assumptions;'
		last_start = 1
		i = 0
		for partition in partitions:
			i +=1
			end = last_start + partition-1
			print >> outNEX, '	charset part{part} = {start}-{end};'.format(part=i, start=last_start, end=end)
			last_start = end + 1
		assert end == int(sites)
		print >> outNEX, 'end;'

def guess_seqtype(seq, gap='-'):
	char_count = Counter(seq.upper())
	print >>sys.stderr, char_count
	nACTG = sum([char_count.get(char, 0) for char in 'ACTG'])
	nACUG = sum([char_count.get(char, 0) for char in 'ACUG'])
	gap = char_count.get('-', 0)
	if 1e2*nACTG / (len(seq)-gap) > 0.8:
		return 'dna'
	elif 1e2*nACUG / (len(seq)-gap) > 0.8:
		return 'rna'
	else:
		return 'protein'
def main():
	subcommand = sys.argv[1]
	if subcommand == 'reTree': # 重新建树，加上bootstrap
		OFdir=sys.argv[2]
		outdir=sys.argv[3]
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		bootstrap_species_tree(OFdir=OFdir, outdir=outdir)
	elif subcommand == 'reTree2': # 只用单拷贝基因建树
		OFdir=sys.argv[2]
		outdir=sys.argv[3]
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		singlecopy_tree(OFdir=OFdir, outdir=outdir)
	elif subcommand == 'orthologs':
		OFdir=sys.argv[2]
		try: sp1, sp2 =sys.argv[3:5]
		except: sp1, sp2 = None, None
		outHomo = sys.stdout
		get_orthologs(OFdir, outHomo, sp1=sp1, sp2=sp2)
	elif subcommand == 'paralogs':
		OFdir=sys.argv[2]
		outHomo = sys.stdout
		get_paralogs(OFdir, outHomo)
	elif subcommand == 'uniqlogs':
		OFdir=sys.argv[2]
		get_unique_logs(OFdir)
	elif subcommand == 'to_cafe':
		OFdir=sys.argv[2]
		outTsv = sys.stdout
		cafe_count(OFdir, outTsv)
	elif subcommand == 'MCScanX_transposed': # 封装MCScanX_transposed【当归】
		OFdir=sys.argv[2]
		tsp, cspp, spmap, gff = sys.argv[3:7]
		cspp = cspp.split(',')
		MCScanX_transposed(OFdir, tsp, cspp, spmap, gff)
	elif subcommand == 'single_copy_cds':   # single_copy_cds alignments【水青树】
		OFdir=sys.argv[2]
		cdsSeqs =sys.argv[3]
		outALN = sys.stdout
		single_copy_cds_align(OFdir, cdsSeqs, outALN)
	elif subcommand == 'single_copy_pep':  # single_copy_pep alignments
		OFdir=sys.argv[2]
		outALN = sys.stdout
		single_copy_pep_align(OFdir, outALN)
	elif subcommand == 'catAln':		  # 合并alignments
		inALNs = sys.argv[2:]
		outALN = sys.stdout
		catAln(inALNs, outALN)
	elif subcommand == 'to_paml':		# 为PAML准备OG【垫柳】
		OFdir =sys.argv[2]
		outDir, cdsSeq = sys.argv[3:5]
		to_paml(OFdir, outDir, cdsSeq)
	elif subcommand == 'species_specific':	# 物种特有基因【垫柳】
		OFdir =sys.argv[2]
		sp =sys.argv[3]
		outTsv = sys.stdout
		species_specific_genes(OFdir, sp, outTsv)
	elif subcommand == 'aln2beast':	# alignment文件转换为BEAST的输入nex文件【醉鱼草】
		inALN = sys.argv[2]
		outNEX = sys.stdout
		aln2beast(inALN, outNEX)
	elif subcommand == 'singlecopy_orthologs':	# 准备两物种间单拷贝OG【醉鱼草】
		OFdir=sys.argv[2]
		try: sp1, sp2 =sys.argv[3:5]
		except: sp1, sp2 = None, None
		outHomo = sys.stdout
		print >> sys.stderr, OFdir, sp1, sp2
		get_singlecopy_orthologs(OFdir, outHomo, sp1=sp1, sp2=sp2)
	elif subcommand == 'to_astral': # 生成单拷贝基因树【润楠】
		OFdir=sys.argv[2]
		pepSeq =sys.argv[3]
		outTrees = sys.stdout
		to_astral(OFdir, pepSeq, outTrees)
	elif subcommand == 'to_hybpiper':
		OFdir=sys.argv[2]
		cdsSeq =sys.argv[3]
		outOGSeq = sys.stdout
		try: species = sys.argv[4]
		except IndexError: species = None
		try: cutoff = float(sys.argv[5])
		except IndexError: cutoff = 0.7
		to_hybpiper(OFdir, cdsSeq, outOGSeq, species=species, min_singlecopy=cutoff)
	elif subcommand == 'singlecopy_stats':
		OFdir=sys.argv[2]
		try: species = sys.argv[3]
		except IndexError: species = None
		to_hybpiper(OFdir, species=species, only_stats=True)
	else:
		raise ValueError('Unknown command: {}'.format(subcommand))

if __name__ == '__main__':
	main()

