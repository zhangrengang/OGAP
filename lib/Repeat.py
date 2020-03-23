import sys, os
import networkx as nx
from Bio import SeqIO
from RunCmdsMP import run_cmd, logger


class RepeatPipeline():
	def __init__(self, genome, tmpdir='/dev/shm/tmp', 
				prefix=None,
				vmatch_opts='-p -d -l 50 -e 3'):
		self.genome = genome
		if prefix is None:
			prefix = os.path.basename(self.genome)
		#self.prefix = prefix
		self.tmpdir = tmpdir
		self.prefix = '{}/{}'.format(self.tmpdir, prefix)
		self.vmatch_opts = vmatch_opts
	def run_vmatch(self):
		cmd = 'mkvtree -db {} -indexname {} -dna -allout -pl'.format(
				self.genome, self.prefix)
		run_cmd(cmd, log=True)
		vmatch_out = '{}.vmatch'.format(self.prefix)
		cmd = 'vmatch {} {} > {}'.format(self.vmatch_opts, self.prefix, vmatch_out)
		run_cmd(cmd, log=True)
		
	def get_seq_no(self):
		return {i:rc.id for i, rc in enumerate(SeqIO.parse(self.genome, 'fasta'))}

class Vmatch:
	def __init__(self, match, d_no=None):
		self.matchfile = matchfile
		self.d_no = d_no
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.matchfile):
			if line.startswith('#'):
				continue
			line = line.strip().split()
			if not len(line) == 10:
				continue
			yield VmatchRecord(line, self.d_no)
	def cluster(self):
		G = VmatchGraph()
		for match in self:
			G.add_edge(match.left_segment, match.right_segment, match=match)
		return G
class VmatchGraph(nx.Graph):
	def __init__(self):
		super(VmatchGraph, self).__init__()
class VmatchRecord:
	def __init__(self, line, d_no=None):
		self.title = ['left_len', 'left_no', 'left_pos', 'type', 
					'right_len', 'rigth_no', 'right_pos',
					'distance', 'evalue', 'score', 'identity']
		self.type = [int, int, int, str,
					int, int, int,
					int, float, float, float]
		self.line = line
		self._set_attr(self.title, self.line, self.type)
		if d_no is not None:
			self.left_id, self.right_id = d_no[self.left_no], d_no[self.right_no]
		self.left_start, self.left_end = self.left_pos+1, self.left_pos+self.left_len
		self.right_start, self.right_end = self.right_pos+1, self.right_pos+self.right_len
		self.left_segment = VmatchSegemnt(self.left_id, self.left_start, self.left_end)
		self.right_segment = VmatchSegemnt(self.right_id, self.right_start, self.right_end)

	def _set_attr(self, keys, values, types):
		for key, value, type in zip(keys, values, types):
			setattr(self, key, type(value))

class VmatchSegemnt:
	def __init__(self, id, start, end):
		self.id, self.start, self.end = id, start, end
	def __hash__(self):
		return hash(self.key)
	def __eq__(self, other):
		if self.key == other.key:
			return True
		return False
	@property
	def key(self):
		return (self.id, self.start, self.end)
