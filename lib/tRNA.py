import sys
from Gff import GffLine, GffExons

AA = {
    'Phe':'F',
    'Leu':'L',
    'Ile':'I',
    'Met':'M',
    'Val':'V',
    'Ser':'S',
    'Pro':'P',
    'Thr':'T',
    'Ala':'A',
    'Tyr':'Y',
    'His':'H',
    'Gln':'Q',
    'Asn':'N',
    'Lys':'K',
    'Asp':'D',
    'Glu':'E',
    'Cys':'C',
    'Arg':'R',
    'Gly':'G',
    'Trp':'W',
}
class tRNAscan():
	def __init__(self, output, struct=None, **kargs):
		self.output = output
		self.struct = struct
		self.kargs = kargs
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for i, line in enumerate(open(self.output)):
			if i < 3:
				continue
			yield tRNAscanRecord(line, **self.kargs)

class tRNAscanRecord():
	def __init__(self, line, min_intron=20):
		self.line = line.strip().split()
		self.title = ['target', 'n', 'start', 'end', 'type', 'anti_codon',
					'intron_start', 'intron_end', 'score']
		self.ctype = [str, int, int, int, str, str, 
					int, int, float]
		self.set_attr()
		self.min_intron = min_intron
		self.anti_codon = self.anti_codon.replace('T', 'U')
	def set_attr(self):
		for key, value, type in zip(self.title, self.line, self.ctype):
			setattr(self, key, type(value))
	@property
	def aa(self):
		try:
			return AA[self.type]
		except KeyError:
			print >>sys.stderr, 'Unknown tRNA product: {}, ignored'.format(self.type)
			return 'X'
	@property
	def name(self):
		return 'trn{}-{}'.format(self.aa, self.anti_codon)
	@property
	def product(self):
		return 'trn{}-{}'.format(self.aa, self.type)
	def update_name(self, name):
		name = name.split('-')
		name[1] = self.anti_codon
		return '-'.join(name)
	def is_trn(self, trn=None):
		if trn is None:
			trn = self.target
		trn_aa = str(trn).split('-')[0][-1].upper()
		if trn_aa == self.aa:
			return True
		return False
	def to_exons(self):
		self.strand = '+'
		start, end = self.start, self.end
		intron_start, intron_end = self.intron_start, self.intron_end
		if start > end:
			self.strand = '-'
			start, end = end, start
			intron_start, intron_end = intron_end, intron_start

		intron_length = intron_end - intron_start + 1
		intron = True
		if intron_length <= self.min_intron:
			intron = False

		if intron:
			regions = [(start, intron_start-1), (intron_end+1, end)]
		else:
			regions = [(start, end)]

		exons = []
		source = 'tRNAscan'
		frame = '.'
		type = 'exon'
		attributes = ''
		for start, end in regions:
			line = [self.target, source, type, start, end, self.score, 
					self.strand, frame, attributes]
			line = GffLine(line)
			exons += [line]
		exons = GffExons(exons)
		return exons
