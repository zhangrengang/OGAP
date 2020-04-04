import sys, re
from Gff import GffLine, GffExons

class PslParser():
	def __init__(self, psl):
		if isinstance(psl, file):
			self.psl = psl
		else:
			self.psl = open(psl)
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in self.psl:
			if not re.compile(r'\d').match(line):
				continue
			yield PslRecord(line)
	def to_exons(self, minintron=200, d_seqs=None, fout=None):
		hits = []
		for i, record in enumerate(self):
			exons = record.to_exons(minintron=minintron)
			exons.id = str(i)
			if fout is not None:
				seq = d_seqs[record.qname]
				exons.seq = exons.extract_seq(seq)
				print >>fout, '>{}\n{}'.format(exons.id, exons.seq)
			hits += [exons]
		return hits

class PslRecord():
	def __init__(self, line):
		line = line.strip().split('\t')
		match, mismatch, rep_match, Ns, \
			qgap_count, qgap_bases, tgap_count, tgap_bases, strand, \
			qname, qsize, qstart, qend, tname, tsize, tstart, tend, \
			block_count, block_sizes, qstarts, tstarts = line[:21]
		match, mismatch, rep_match, Ns, \
			qgap_count, qgap_bases, tgap_count, tgap_bases, \
			qsize, qstart, qend, tsize, tstart, tend, \
			block_count = map(int, [
		match, mismatch, rep_match, Ns, \
			qgap_count, qgap_bases, tgap_count, tgap_bases, \
			qsize, qstart, qend, tsize, tstart, tend, \
			block_count])

#		block_sizes = map(int, block_sizes.strip(',').split(','))
#		qstarts = map(int, qstarts.strip(',').split(','))
#		tstarts = map(int, tstarts.strip(',').split(','))
		self.match, self.mismatch, self.rep_match, self.Ns, \
			self.qgap_count, self.qgap_bases, \
			self.tgap_count, self.tgap_bases, self.strand, \
			self.qname, self.qsize, self.qstart, self.qend, \
			self.tname, self.tsize, self.tstart, self.tend, \
			self.block_count, self._block_sizes, self._qstarts, self._tstarts = \
		match, mismatch, rep_match, Ns, \
			qgap_count, qgap_bases, tgap_count, tgap_bases, strand, \
			qname, qsize, qstart, qend, tname, tsize, tstart, tend, \
			block_count, block_sizes, qstarts, tstarts
		self.qlclip = self.qstart
		self.qrclip = self.qsize - self.qend
		self.tlclip = self.tstart
		self.trclip = self.tsize - self.tend
	
	@property
	def line(self):
		line = [self.match, self.mismatch, self.rep_match, self.Ns,
            self.qgap_count, self.qgap_bases,
            self.tgap_count, self.tgap_bases, self.strand,
            self.qname, self.qsize, self.qstart, self.qend,
            self.tname, self.tsize, self.tstart, self.tend,
            self.block_count, self.block_sizes, self.qstarts, self.tstarts]
		line = map(str, line)
		return '\t'.join(line)
	@property
	def block_sizes(self):
		return map(int, self._block_sizes.strip(',').split(','))
	@property
	def qstarts(self):
		return map(int, self._qstarts.strip(',').split(','))
	@property
	def tstarts(self):
		return map(int, self._tstarts.strip(',').split(','))
	def to_exons(self, minintron=200):
		starts = self.qstarts
		name = self.qname
		score = self.qscore
		exons = []
		last_end = None
		for start, block_size in zip(starts, self.block_sizes):
			start = start + 1
			end = start + block_size
			line = [name, 'blat', 'exon', start, end, score, 
					self.strand, '.', '']
			if last_end is None or start - last_end > minintron:
				exons += [GffLine(line)]
			else:
				exons[-1].end = end
			last_end = end
		exons = GffExons(exons)
		return exons
	@property
	def qscore(self):
		return self.match - self.mismatch - self.qgap_count - self.qgap_bases - self.tgap_count - self.tgap_bases - self.block_count - self.tlclip - self.trclip

	@property
	def score(self):
		return self.match - self.mismatch - self.qgap_count - self.qgap_bases - self.tgap_count - self.tgap_bases - self.block_count - self.qlclip - self.qrclip
	@property
	def qcov(self):
		return
	@property
	def tcov(self):
		return
def test():
	psl = sys.stdin
	for rc in PslParser(psl):
		print rc.line, rc.score
if __name__ == '__main__':
	test()
