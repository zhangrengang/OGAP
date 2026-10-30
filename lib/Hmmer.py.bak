import sys
import itertools
import copy
import networkx as nx
try:
	from math import inf
except ImportError:
	inf = float("inf")
from lazy_property import LazyWritableProperty as lazyproperty
from Gff import GffLine, GtfExons
from get_record import get_records

class HmmSearch(object):
	def __init__(self, hmmout, hmmfmt='domtbl'):
		self.hmmout = hmmout
		self.hmmfmt = hmmfmt
	def __iter__(self):
		return self.parse()
	def parse(self):
		for line in open(self.hmmout):
			if line.startswith('#'):
				continue
			if self.hmmfmt == 'domtbl':
				yield HmmSearchDomHit(line)
	def get_best_hit(self, score=False):
		best_hit = None
		d_cov = {}
		for record in self:
			if best_hit is None:
				best_hit = record
			elif not score and record.edit_score > best_hit.edit_score:
				best_hit = record
			elif score and record.score > best_hit.score:
				best_hit = record
			try: d_cov[record.tname] += record.hmmcov
			except KeyError: d_cov[record.tname] = record.hmmcov
		if best_hit is not None:
			best_hit.cov = d_cov[best_hit.tname]
		return best_hit
	def get_hit_nuclseqs(self, inseq, outseq):
		nucl_names = []
		for hit in self:
			hit.split_name()
			nucl_names += [hit.nucl_name]
		get_records(inseq, outseq, nucl_names, type='fasta')
	def get_gene_structure(self, d_length, max_copies=4, min_part=1,
			min_hmmcov=0, min_contained_cov=0.95,
			min_cov=80, min_ratio=0.9, seq_type='prot', flank=1000):
		graph = HmmStructueGraph()
		for hit in self:
			if hit.hmmcov < min_hmmcov:
				continue
			hit.convert_to_nucl_coord(d_length, seq_type=seq_type)
	#		print >>sys.stderr, hit.nucl_hit
			graph.add_node(hit)
		print >>sys.stderr, 'Contained nodes:'
		graph.prune_contained(min_cov=min_contained_cov)
		graph.link_nodes()
		print >>sys.stderr, 'Nodes:'
		for node in sorted(graph.nodes(), key=lambda x: x.hmmcoord):
			print >>sys.stderr, node.short
		print >>sys.stderr, 'Edges:'
		for n1, n2 in sorted(graph.edges(), key=lambda x: (x[0].hmmcoord, x[1].hmmcoord)):
			print >>sys.stderr, n1.short, '-',  n2.short, graph[n1][n2]['dist']
		graph.prune_graph()
		graph.break_circle()
	#	print >>sys.stderr, graph.nodes()

		print >>sys.stderr, 'after pruned:'
		for n1, n2 in sorted(graph.edges(), key=lambda x: (x[0].hmmcoord, x[1].hmmcoord)):
			print >>sys.stderr, n1.short, '-',  n2.short

		hmmcovs = []
		copies = []
		for path in graph.linearize_path():
	#		print >>sys.stderr, 'path:', path
			if path.hmmcov < min_cov:
				print >>sys.stderr, 'discarded path for low coverage:', path
				continue
			print >>sys.stderr, path, len(path.group_nodes()), path.hmmcov
			parts = path.get_parts(d_length, flank=flank)
		#	print >>sys.stderr, parts
#			seq = parts.combine_seq()
			copies += [parts]
			hmmcovs += [path.hmmcov]
		# filter out too many parts
		if len(copies) > 0:
	#		min_part_number = min(map(len, copies))
			max_hmmcov = max(hmmcovs)
			cov_cutoff = max_hmmcov * min_ratio
			best_copies = [parts for parts, hmmcov in zip(copies, hmmcovs) \
								if not hmmcov < max_hmmcov*0.99]
			min_part_number = min(map(len, best_copies))

			i = 0
			for parts, hmmcov in sorted(zip(copies, hmmcovs), key=lambda x:-x[1]):
				if len(parts) > min_part_number+min_part:
					continue
				if hmmcov < cov_cutoff:
					continue
				i += 1
				if i > max_copies and hmmcov < max_hmmcov*0.99:
					continue
		#		print >>sys.stderr, min_part_number, parts
				yield parts		# maybe multi-copy

class HmmPath():
	def __init__(self, path):
		self.path = path
	def __iter__(self):
		return iter(self.path)
	def __len__(self):
		return len(self.path)
	def __getitem__(self, index):
		if isinstance(index, int):
			return self.path[index]
		else:
			return self.__class__(path=self.path[index])
	def __str__(self):
		return ' -> '.join(['{}:{}-{}({})'.format(p.nucl_name, p.nucl_alnstart, 
				p.nucl_alnend, p.strand) for p in self.path])
	def group_nodes(self, max_dist=20000, min_dist=-200):
		if len(self) == 1:
			return [self]
		groups = [[self[0]]]
		for node in self[1:]:
			last_node = groups[-1][-1]
			if node.strand != last_node.strand or node.nucl_name != last_node.nucl_name:
				groups += [[node]]
			else:
				if node.strand == '-':
					dist = last_node.nucl_start - node.nucl_end #start
				else:	# +
					dist = node.nucl_start - last_node.nucl_end #start
				if min_dist < dist < max_dist:
					groups[-1] += [node]
				else:
					groups += [[node]]
		return [HmmPath(path=group) for group in groups]
	def get_parts(self, d_length, flank=1000):
		parts = []	# parts = gene/copy
		for path in self.group_nodes():		# node = exon
			first, last = path[0], path[-1]
			assert first.strand == last.strand and first.nucl_name == last.nucl_name
			seqid = first.nucl_name
			strand = first.strand
			if strand == '-':
				first, last = last, first
			start, end = first.nucl_start-1, last.nucl_end
				
			start = max(0, start-flank)
			end = min(d_length[seqid], end+flank)
			part = SeqPart(seqid, start, end, strand)	# part = continous exons, 0-based
			part.path = path
			parts += [part]
		parts = SeqParts(parts)
#		seq = parts.combine_seq()
		return parts
	@property
	def hmmcov(self):
		first, last = self.path[0], self.path[-1]
		return 1e2* (last.hmmend - first.hmmstart +1) / first.hmm_length
class SeqParts():
	def __init__(self, parts):
		self.parts = parts
	def __iter__(self):
		return iter(self.parts)
	def __len__(self):
		return len(self.parts)
	def __str__(self):
		try:
			return str(self.id)
		except AttributeError:
			return self.to_str()
	def to_str(self):
		return ' -> '.join([str(part) for part in self])
	def combine_seq(self, d_seqs):
		seqs = []
		for part in self:
			seq = part.get_seq(d_seqs)
			if part.strand == '-':
				seq = seq.reverse_complement()
			seqs += [str(seq)]
		return ''.join(seqs)
	def write_seq(self, d_seqs, fout):
		seq = self.combine_seq(d_seqs)
		try:
			blocks = [len(part) for part in self]
		except ValueError:
			print >> sys.stderr, [str(part) for part in self]

		desc = blocks
		print >> fout, '>{} {}\n{}'.format(self, desc, seq)
	def merge(self):
		starts  = [part.start  for part in self]
		ends    = [part.end    for part in self]
		strands = [part.strand for part in self]
		chroms  = [part.seqid  for part in self]
		assert len(set(chroms)) == 1
		assert len(set(strands)) == 1
		chrom  = chroms[0]
		strand = strands[0]
		start  = min(starts + ends)
		end    = max(starts + ends)
		return SeqPart(chrom, start, end, strand)

	def link_part(self, max_dist=10000):
		'''link adj parts'''
		# "--> -->" -> "----->"
		parts = list(self)
		if len(parts) < 2:
			return self
		groups = [[parts[0]]]
		for part in parts[1:]:
			last_part = groups[-1][-1]
			dist = max_dist
			if last_part.seqid == part.seqid and last_part.strand == part.strand:
				if last_part.strand == '-':
					dist = last_part.start - part.end
				else:
					dist = part.start - last_part.end
			if 0 < dist < max_dist:
				groups[-1] += [part]
			else:
				groups += [[part]]
		new_parts = []
		for group in groups:
			if len(group) == 1:
				new_parts += group
				continue
			part = SeqParts(group).merge()	# SeqPart
			new_parts += [part]
		seqparts = copy.deepcopy(self)
		seqparts.parts = new_parts
		seqparts = seqparts.link_overlaped_part()
		return seqparts
	def link_overlaped_part(self):
		parts = list(self)
		while True:
			i = 0
			remove = []
			replace = {}
			for part1, part2 in itertools.combinations(parts, 2):
				if part1.overlaps(part2):
					print >>sys.stderr, '{} overlaps {}'.format(part1, part2)
					i += 1
					if len(part1) >= len(part2):
						remove += [part2]
						replace[part1] = SeqParts([part1, part2]).merge()
					else:
						remove += [part1]
						replace[part2] = SeqParts([part1, part2]).merge()
			if i == 0:
				break
			remove = set(remove)
			new_parts = []
			for part in parts:
				if part in remove:
					continue
				if part in replace:
					part = replace[part]
				new_parts += [part]
			parts = new_parts
		seqparts = copy.deepcopy(self)
		seqparts.parts = parts
		return seqparts
	@lazyproperty
	def coord_map(self):
		d = {}
		last_start = 0	# 0-based
		for part in self:
			start = last_start
			end = last_start + len(part)
			d[part] = (start, end)
			last_start = end
		return d 
	def map_coord(self, exons):
		new_exons = []
		d_bins = {}
#		print >>sys.stderr, ' -> '.join(['{chrom}:{start}-{end}({strand})'.format(**exon.__dict__) \
 #                               for exon in exons])
		for old_exon in exons:
			for part in self:
				start, end = self.coord_map[part]	# 0-based
				exon = copy.deepcopy(old_exon)
				#exon.id = exon.chrom
				exon.chrom = part.seqid
				exon_start, exon_end = exon.start, exon.end
				if max(exon.start, exon.end) <= start or min(exon.start, exon.end) > end:
					# ---|----|---  part
					# --        --  exon	# no overlap
					continue
				elif exon.start > start and exon.end <= end: 
					# |-----|   # part  0-based
					#   ---     # exon  1-based 
					if part.strand == '-':
						exon.start = part.start + end - exon_end +1
						exon.end = part.start + end - exon_start +1
						exon.strand = '+' if exon.strand == '-' else '-'
					else:
						exon.start = part.start + exon_start - start
						exon.end = part.start + exon_end - start
				elif end >= exon.start > start and exon.end > end:
					# |-----|----	# part
					#     -----		# exon
					print >>sys.stderr, 'here'
					if part.strand == '-':
						exon.start = part.start+1	#part.end
						exon.end = part.start + end - exon_start +1
						exon.strand = '+' if exon.strand == '-' else '-'
					else:
						exon.start = part.start + exon_start - start
						exon.end = part.end
				elif  exon.start <= start and exon.end > end:
					# ---|----|---	part
					#   --------	exon
					if part.strand == '-':
						exon.start = part.end
						exon.end = part.start + 1
						exon.strand = '+' if exon.strand == '-' else '-'
					else:
						exon.start = part.start + 1
						exon.end = part.end
				elif exon.start <= start and start < exon.end <= end:
					# ----|-----|	part	
					#   ------		exon
					if part.strand == '-':
						exon.start = part.start + end - exon_end +1
						exon.end = part.end #part.start + 1
						exon.strand = '+' if exon.strand == '-' else '-'
					else:
						exon.start = part.start + 1
						exon.end = part.start + exon_end - start
				else:	# no overlap
					continue
				new_exons += [exon]
				try: d_bins[part] += [exon]
				except KeyError: d_bins[part] = [exon]
			#	print >>sys.stderr, part, exon_start, exon_end, start, end, exon
		new_parts = copy.deepcopy(self)
		parts = []
		for part in self:
			try: exons = GtfExons(d_bins[part])
			except KeyError: continue
			seqid, start, end, strand = exons.get_region()
			parts += [SeqPart(seqid, start, end, strand)]
		new_parts.parts = parts
		return GtfExons(new_exons), new_parts
	def to_exons(self, min_intron=20):
		source, score, frame = 'hmmsearch', '.', '.'
		type = 'exon'
		attributes = ''
		exons = []
		for part in self:
			paths = part.path.group_nodes(max_dist=min_intron)
			for path in paths:
			#	print >>sys.stderr, path
				first_hit, last_hit = path[0], path[-1]
				strand = first_hit.strand
				if strand == '-':
					start, end = last_hit.nucl_start, first_hit.nucl_end
				else:
					start, end = first_hit.nucl_start, last_hit.nucl_end
				chrom = first_hit.nucl_name
				line = [chrom, source, type, start, end, score, strand, frame, attributes]
				exon = GffLine(line)
				exons += [exon]
		return GtfExons(exons)

class SeqPart():
	'''Segment'''
	def __init__(self, seqid, start, end, strand, seq=None):	# 0-based
		self.seqid = seqid
		self.start, self.end, self.strand = start, end, strand
		self.seq = seq
	def __str__(self):
		return '{seqid}:{start}-{end}({strand})'.format(**self.__dict__)
	def __len__(self):
		return self.end-self.start
	def __hash__(self):
		return hash(self.key)
	@property
	def key(self):
		return (self.seqid, self.start, self.end, self.strand)
	def get_seq(self, d_seqs):
		return d_seqs[self.seqid][self.start:self.end]
	def overlaps(self, other):
		if self.seqid != other.seqid or self.strand != other.strand:
			return False
		return max(0, min(self.end, other.end)-max(self.start, other.start))

class HmmSearchDomHit:
	def __init__(self, line):
		self.line = line.strip().split()
		self.title = ['tname', 'tacc', 'tlen', 'qname', 'qacc', 'qlen',
					'full_evalue', 'full_score', 'full_bias',
					'domi', 'domn', 'cevalue', 'ievalue', 'score', 'bias',
					'hmmstart', 'hmmend',	# 1-based, 1-100
					'alnstart', 'alnend',
					'envstart', 'envend',
					'acc']
		self.ctype = [str, str, int, str, str, int,
					float, float, float,
					int, int, float, float, float, float,
					int, int, int, int, int, int,
					float]
		self.set_attr()
		self.hmm_name = self.qname
		self.hmm_length = self.qlen
	def set_attr(self):
		for key, value, type in zip(self.title, self.line, self.ctype):
			setattr(self, key, type(value))
	def __hash__(self):
		return hash(self.key)
	def __eq__(self, other):
		if self.key == other.key:
			return True
		return False
	@property
	def edit_score(self):
		return self.full_score * (1.0*(self.qlen- abs(self.tlen-self.qlen)) / self.qlen)
	@property
	def hmmcov(self):
		return round(1e2*(self.hmmend - self.hmmstart + 1) / self.hmm_length, 1)
	@property
	def hmmcovlen(self):
		return self.hmmend - self.hmmstart + 1
	@property
	def key(self):
		return (self.tname, self.qname, self.hmmstart, self.hmmend, self.alnstart, self.alnend)
	def contained(self, other):
		if other.hmmstart >= self.hmmstart and other.hmmend <= self.hmmend:
			return 1.0*other.hmmcovlen / self.hmmcovlen
		return 2
	def link_hits(self, other, max_dist=60, min_ovl=-60, min_flank_dist=2):
		'''--->	 self
			  --->  other
		'''
		# contains
		# ------->	self	 ---->
		#   ---->   other  -------->
		if other.hmmstart-10 >= self.hmmstart and other.hmmend+10 <= self.hmmend:
			return False
		if self.hmmstart-10 >= other.hmmstart and self.hmmend+10 <= other.hmmend:
			return False
		left_dist = other.hmmstart - self.hmmstart + 1
		mid_dist  = other.hmmstart - self.hmmend + 1
		right_dist = other.hmmend - self.hmmend + 1
		#print self.short, other.short, left_dist, mid_dist, right_dist
		#print >> sys.stderr, self.short, other.short, left_dist, mid_dist, right_dist
		if min_ovl <= mid_dist <= max_dist:
#			if min(left_dist, right_dist) < min_flank_dist:
			if left_dist >= min_flank_dist:
				return mid_dist
		# rRNA algnment
		if self.nucl_name == other.nucl_name and self.strand == other.strand:
			left_dist = other.nucl_alnstart - self.nucl_alnstart + 1
			mid_dist  = other.nucl_alnstart - self.nucl_alnend + 1
			right_dist = other.nucl_alnend - self.nucl_alnend + 1
		#	print >> sys.stderr, self.short, other.short, (self.alnstart,self.alnend), (other.alnstart,other.alnend), left_dist, mid_dist, right_dist
			if min_ovl*3 <= mid_dist <= max_dist*3:
				if left_dist >= min_flank_dist:
					return mid_dist
#		print >> sys.stderr, self.short, other.short, left_dist, mid_dist, right_dist
		return False
			
	def split_name(self):
		if not 'nucl_frame' in self.__dict__:
			self.nucl_name, self.nucl_frame = self.split_translated_id(self.tname)	# t
	def convert_to_nucl_coord(self, d_length, seq_type='prot'):
		self.split_name()
		nucl_length = d_length[self.nucl_name]
		self.strand, self.frame = self.parse_frame(self.nucl_frame)
		if seq_type == 'prot':
			if self.strand == '-':
				self.nucl_alnstart = nucl_length - (self.alnend * 3 + self.frame) + 1
				self.nucl_envstart = nucl_length - (self.envend * 3 + self.frame) + 1
				self.nucl_alnend   = nucl_length - ((self.alnstart-1) * 3 + self.frame)
				self.nucl_envend   = nucl_length - ((self.envstart-1) * 3 + self.frame)
			else: # +
				self.nucl_alnstart = (self.alnstart-1) * 3  + self.frame + 1
				self.nucl_envstart = (self.envstart-1) * 3  + self.frame + 1
				self.nucl_alnend   = self.alnend * 3 + self.frame
				self.nucl_envend   = self.envend * 3 + self.frame
		else:
			if self.strand == '-':
				self.nucl_alnstart = nucl_length - self.alnend + 1
				self.nucl_envstart = nucl_length - self.envend + 1
				self.nucl_alnend   = nucl_length - self.alnstart + 1
				self.nucl_envend   = nucl_length - self.envstart + 1
			else:
				self.nucl_alnstart = self.alnstart
				self.nucl_envstart = self.envstart
				self.nucl_alnend   = self.alnend
				self.nucl_envend   = self.envend
		self.nucl_length = nucl_length
		self.nucl_start, self.nucl_end = self.nucl_alnstart, self.nucl_alnend
	@property
	def nucl_hit(self):
		return (self.nucl_name, self.nucl_length, self.nucl_alnstart, self.nucl_alnend,
				self.strand, self.frame, 
				self.hmm_name, self.hmm_length, self.hmmstart, self.hmmend)
	@property
	def short(self):
		return (self.nucl_alnstart, self.nucl_alnend, self.strand, self.hmmstart, self.hmmend)
	@property
	def hmmcoord(self):
		return (self.hmmstart, self.hmmend)
	def __str__(self):
		return str(self.short)
	def get_hmm_distance(self, other):
		if self.hmmstart > other.hmmstart: # o--> s-->
			return self.hmmend - other.hmmstart
		else:	# s--> o-->
			return other.hmmend - self.hmmstart
	def get_distance(self, other):
		if self.nucl_name != other.nucl_name:
			return self.get_hmm_distance(other) * 3e6
		if self.strand == other.strand:
			if self.strand == '-':  # <-- <--
				if self.nucl_alnstart > other.nucl_alnstart:  # <--o <--s
					return self.nucl_alnend - other.nucl_alnstart
				else:										  # <--s <--o
					return self.nucl_alnstart + self.nucl_length - other.nucl_alnend
			else:
				if self.nucl_alnstart < other.nucl_alnstart:  # s--> o-->
					return other.nucl_alnend - self.nucl_alnstart
				else:										  # o--> s-->
					return self.nucl_length - self.nucl_alnend + other.nucl_alnstart
		else:  # s--> <--o, <--o s-->
			return self.get_hmm_distance(other) * 1e6 #self.nucl_length #min(self.nucl_alnstart + other.nucl_alnstart, 
					#self.nucl_length + other.nucl_length - self.nucl_alnend - other.nucl_alnend)
	def to_gff(self):
		line = [self.nucl_name, 'HMMsearch', 'protein_match', 
				self.nucl_alnstart, self.nucl_alnend, self.score, self.strand, 0,
				]
	def split_translated_id(self, name):
		name = name.split('|')
		nucl_name = '|'.join(name[:-1])
		nucl_frame = name[-1]
		return nucl_name, nucl_frame
	def parse_frame(self, string):	# frame=0-2
		if string.startswith('rev'):
			strand = '-'
		else:
			strand = '+'
		try:
			frame = int(string[-1]) -1
		except ValueError: 
			frame = '.'
		return strand, frame
class DiGraph(nx.DiGraph):
	def __init__(self, data=None, **attr):
		super(DiGraph, self).__init__(data, **attr)
	@property
	def sources(self):
		return [ node for node in self.nodes() if not self.predecessors(node)]
	@property
	def targets(self):
		return [ node for node in self.nodes() if not self.successors(node)]

class HmmStructueGraph(DiGraph):
	def __init__(self, data=None, **attr):
		super(HmmStructueGraph, self).__init__(data, **attr)
	def build_from_hmm_hits(self, hits):	# HmmSearchDomHit
		for i, hit in enumerate(hits):
			self.add_node(hit)
	def link_nodes(self):
		for hit1, hit2 in itertools.combinations(self.nodes(), 2):
			dist_12 = hit1.link_hits(hit2)
			dist_21 = hit2.link_hits(hit1)
			#print >>sys.stderr, hit1.short, hit2.short, dist_12, dist_21
			if not dist_12 is False:
				self.add_edge(hit1, hit2, dist=dist_12)
			if not dist_21 is False:
				self.add_edge(hit2, hit1, dist=dist_21)
	def prune_contained(self, min_cov=0.6):
		tobe_removed = set([])
		for hit1, hit2 in itertools.combinations(self.nodes(), 2):
			if hit1.contained(hit2) < min_cov:
				print >>sys.stderr, hit1.short, '-',  hit2.short
				tobe_removed.add(hit2)
			elif hit2.contained(hit1) < min_cov:
				print >>sys.stderr, hit2.short, '-',  hit1.short
				tobe_removed.add(hit1)
		
		for node in tobe_removed:
			self.remove_node(node)
	def break_circle(self):
		for circle in nx.simple_cycles(self):
#			print circle
			sg = self.subgraph(circle)
			max_dist_edge = max(sg.edges(), key=lambda x: sg.get_edge_data(*x))
			n1, n2 = max_dist_edge
			print >>sys.stderr, 'break circle', n1.short, '-',  n2.short
			self.remove_edge(*max_dist_edge)
	def prune_graph(self):
		# remove bridged edge
		for node in self.nodes():
			predecessors = self.predecessors(node)
			successors = self.successors(node)
			for predecessor, successor in itertools.product(predecessors, successors):
				if self.has_edge(predecessor, successor):
					self.remove_edge(predecessor, successor)
					n1, n2 = predecessor, successor
					print >>sys.stderr, 'remove bridged edge:', n1.short, '-',  n2.short
		# remove distant edge
		for node in self.nodes():
			predecessors = set(self.predecessors(node))
			if len(predecessors) < 2:
				continue
			pred_successors = []
			products = []
			for predecessor in predecessors:
				successors = self.successors(predecessor)
				pred_successors += successors
				for successor in successors:
					products += [(predecessor, successor)]
			pred_successors = set(pred_successors)
			if len(pred_successors) < 2:
				continue
			length = min(len(set(predecessors)), len(set(pred_successors)))
#			print >>sys.stderr, length
			combinations = list(ExcludeProduct(products, length))
			distances = []
			for combination in combinations:
				sum_distance = 0
				for hit1, hit2 in combination:
					distance = hit1.get_distance(hit2)
					sum_distance += distance
				distances += [sum_distance]
#			print >>sys.stderr, products, distances, combinations
			if not distances:
				continue
			min_dist, combination = min(zip(distances, combinations), key=lambda x:x[0])
			if min_dist == inf:
				continue
			for edge in products:
				if edge in set(combination):
					continue
				self.remove_edge(*edge)
				n1, n2 = edge
				print >>sys.stderr, 'remove distant edge:', n1.short, '-',  n2.short
	def linearize_path(self):
		for source, target in itertools.product(self.sources, self.targets):
			for path in nx.all_simple_paths(self, source, target):
				yield HmmPath(path)
		for node in self.nodes():
			if not self.predecessors(node) and  not self.successors(node):
				yield HmmPath([node])
	def cluster_genes(self):
		pass
class ExcludeProduct():
	def __init__(self, product, cutoff=None):
		self.product = product
		self.cutoff = cutoff
	def __iter__(self):
		return iter(self._exlude())
	def _exlude(self):
		for combinations in itertools.combinations(self.product, self.cutoff):
			include = False
			for p1, p2 in itertools.combinations(combinations, 2):
				if self.is_not_exclude(p1, p2):
					include = True
					break
			if include:
				continue
			yield list(combinations)
	def is_not_exclude(self, p1, p2):
		for v1, v2 in zip(p1, p2):
			if v1 == v2:
				return True
		return False
class HmmCluster(object):
	def __int__(self, hmmout, d_length=None): # only for nucl
		self.hmmout = hmmout
		self.d_length = d_length
	def cluster_by_hmm(self):
		d_cluster = OrderedDict()
		for hit in HmmSearch(self.hmmout):
			hit.split_name()
			nucl_length = self.d_length[hit.nucl_name]
			hit.convert_to_nucl_coord(nucl_length)
			key = hit.qname
			try: d_cluster[key] += [hit]
			except KeyError: d_cluster[key] = [hit]
		return d_cluster
	def cluster_hits(self, hits, max_diff=50):
		dG = nx.DiGraph()
		for i, hit in enumerate(hits):
			pass
	def xff(self):
		if len(records) == 1:
			return records
		records = sorted(records, key=lambda x:x.hmmstart)
		best_idx, best_rc = self.maxscore(records)
		new_rcs = [best_rc]
		for rc in records[:best_idx][::-1]:	 # <<- left extand
			right_rc = new_rcs[0]
			mal_pos = right_rc.hmmstart - rc.hmmend
			if abs(mal_pos) <= max_mal:
				if mal_pos <= 0:
					diff = 1-mal_pos
					rc.hmmend -= diff
					rc.alnend -= diff
				new_rcs = [rc] + new_rcs
		for rc in records[best_idx:]:		   # ->> right extand
			left_rc = new_rcs[-1]
			mal_pos = rc.hmmstart - left_rc.hmmend
			if abs(mal_pos) <= max_mal:
				if mal_pos <= 0:
					diff = 1-mal_pos
					rc.hmmstart += diff
					rc.alnstart += diff
				new_rcs += [rc]
		return HmmClusterRecord(new_rcs)
	def maxscore(self, records):
		for i, rc in enumerate(records):
			if i == 0:
				best = (i, rc)
				continue
			if rc.score > best[1].score:
				best = (i, rc)
		return best

class HmmClusterRecord(object):
	def __init__(self, records):
		self.records = records
		self.score = sum([rc.score for rc in records])
		self.hmmstart = records[0].hmmstart
		self.hmmend = records[-1].hmmend
		self.alnstart = records[0].alnstart
		self.alnend = records[-1].alnend
		self.tlen = records[0].tlen
		self.hmmcov = round(1e2*(self.hmmend - self.hmmstart + 1) / self.tlen, 1)
		self.evalue = multi(*[rc.evalue for rc in records])

