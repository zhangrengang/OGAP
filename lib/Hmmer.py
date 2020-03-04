import sys
import itertools
import networkx as nx
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
	def get_hit_nuclseqs(self, inseq, outseq):
		nucl_names = []
		for hit in self:
			hit.split_name()
			nucl_names += [hit.nucl_name]
		get_records(inseq, outseq, nucl_names, type='fasta')

class HmmSearchDomHit(object):
	def __init__(self, line):
		self.line = line.strip().split()
		self.title = ['tname', 'tacc', 'tlen', 'qname', 'qacc', 'qlen'
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
	def hmmcov(self):
		return round(1e2*(self.hmmend - self.hmmstart + 1) / self.qlen, 1)
	@property
	def key(self):
		return (self.tname, self.qname, self.hmmstart, self.hmmend, self.alnstart, self.alnend)
	def link_hits(self, other, max_dist=25, min_ovl=-50, min_flank_dist=10):
		left_dist = other.hmmstart - self.hmmstart + 1
		mid_dist  = self.hmmend - other.hmmstart + 1
		right_dist = other.hmmend - self.hmmend + 1
		if min_ovl <= mid_dist <= max_dist:
			return mid_dist
			if min(left_dist, right_dist) < min_flank_dist:
				return False
		return False
			
	def split_name(self):
		if not 'nucl_frame' in self.__dict__:
			self.nucl_name, self.nucl_frame = self.split_translated_id(self.tname)	# t
	def convert_to_nucl_coord(self, nucl_length):
		self.split_name()
		self.strand, self.nucl_frame = self.parse_frame(self.nucl_frame)
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
		frame = int(string[-1]) -1
		return strand, frame

class StructueGraph(nx.DiGraph):
	def __init__(self):
		super(SeqGraph, self).__init__()
	def build_from_hmm_hits(self, hits):
		for i, hit in enumerate(hits):
			self.add_node(hit)
		for hit1, hit2 in itertools.combinations(hits, 2):
			dist_12 = hit1.link_hits(hit2)
			dist_21 = hit2.link_hits(hit1)
			if not dist_12 is False:
				self.add_edge(hit1, hit2, dist=dist_12)
			if not dist_21 is False:
				self.add_edge(hit2, hit1, dist=dist_21)
	def prune_graph(self):
		pass
	def cluster_genes(self):
		pass
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

