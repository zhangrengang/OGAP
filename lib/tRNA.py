import sys
import re
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
	def __init__(self, output=None, struct=None, **kargs):
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
	def __init__(self, line=None, min_intron=20):
		if line is None:
			return
		self.line = line.strip().split()
		self.title = ['target', 'n', 'start', 'end', 'type', 'anti_codon',
					'intron_start', 'intron_end', 'score']
		self.ctype = [str, int, int, int, str, str, 
					int, int, float]
		self.set_attr()
		self.min_intron = min_intron
		self.anti_codon = self.set_anti_codon()
	def set_attr(self):
		self._set_attr(self.title, self.line, self.ctype)
	def _set_attr(self, keys, values, types):
		for key, value, type in zip(keys, values, types):
			setattr(self, key, type(value))
	def set_anti_codon(self):
		return self.anti_codon.replace('T', 'U')
	@property
	def aa(self):
		try:
			return AA[self.type]
		except KeyError:
			print >>sys.stderr, 'Unknown tRNA product: {}, ignored'.format(self.type)
			return 'XX'
	@property
	def name(self):
		return 'trn{}-{}'.format(self.aa, self.anti_codon)
	@property
	def product(self):
		return 'trn{}-{}'.format(self.aa, self.type)
	def update_name(self, name, anti_codon=None):
		name = str(name)
		prefix = self.get_prefix(name)
		aa = self.get_aa(name)
		if anti_codon is None:
			anti_codon = self.get_codon(name)
		#if all([prefix, aa, anti_codon]):
		#	return name
		if aa is None:
			aa = self.aa
		if anti_codon is None:
			anti_codon = self.anti_codon
		rename = '{}{}-{}'.format(prefix, aa, anti_codon)
		if self.is_cp(name):
			rename += '-cp'
		return rename

	def get_prefix(self, name):
		name = str(name)	# (trn)M-CAU-cp
		pattern = r'(trn.*?)[A-Z]\-?'
		match = re.compile(pattern).match(name)
		if match:
			return match.groups()[0]
		return 'trn'
	def get_aa(self, name):
		name = str(name)	# trn(M)-CAU-cp
		pattern1 = r'trn.*?([A-Z])\-?'
		match = re.compile(pattern1).match(name)
		if match:
			return match.groups()[0]
		for key in AA:		# tRNA-(Tyr)-GUA
			match = re.compile(key).search(name)
			if match:
				return AA[key]
		return None
	def get_codon(self, name):
		name = str(name)	# trnM-(CAU)-cp
		pattern = r'\-([atcgu]{3})\-?'
		match = re.compile(pattern, re.I).search(name)
		if match:
			return match.groups()[0] #.upper()
		return None
	def is_cp(self, name):
		name = str(name)	# trnM-CAU-(cp)
		pattern = r'\-(cp)'
		match = re.compile(pattern, re.I).search(name)
		if match:
			return True
		return False
	def is_trn(self, trn=None):
		if trn is None:
			trn = self.target
		trn_aa = self.get_aa(trn)
		trn_anti_codon = self.get_codon(trn)
		if trn_anti_codon is None and trn_aa == self.aa:
			return True
		elif trn_aa == self.aa and trn_anti_codon == self.anti_codon:
			return True
		elif trn_anti_codon == self.anti_codon:
			# trnI-cau vs trnM-cau
			print >>sys.stderr, 'same tRNA anti-codon: {} with different product: {} vs {}, retrained'.format(
							trn_anti_codon, trn_aa, self.aa)
			return True
		elif trn_aa == self.aa:
			print >>sys.stderr, 'same tRNA product: {} with different anti-codon: {} vs {}, ignored'.format(
							self.type, trn_anti_codon, self.anti_codon)
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
			line = [self.target, source, type, start, end, round(self.score,1), 
					self.strand, frame, attributes]
			line = GffLine(line)
			exons += [line]
		exons = GffExons(exons)
		return exons

class tRNAscanStructs():
	def __init__(self, struct, **kargs):
		self.struct = struct
	def __iter__(self):
		return self._parse()
	def _parse(self):
		lines = []
		for i, line in enumerate(open(self.struct)):
			if not line.strip():
				yield tRNAscanStruct(lines)
				lines = []
			lines += [line]
			
class tRNAscanStruct(tRNAscanRecord):
	def __init__(self, lines=None):
		self.lines = lines
		self._parse()
	def _parse(self):
		for line in self.lines:
			line = line.strip()
			if line.startswith('Type'):
				self._parse_type(line)
			elif line.startswith('Seq'):
				self._parse_seq(line)
			elif line.startswith('Str'):
				self._parse_struct(line)
			elif line.endswith('bp'):
				self._parse_trna(line)
	def _parse_trna(self, line):
		patten = r'(\S+)\.trna(\d+)\s+\((\d+)\-(\d+)\)\s+Length:\s+(\d+)\s+bp'
		match = re.compile(patten).match(line).groups()
		keys = ['target', 'n', 'start', 'end', 'length']
		types = [str, int, int, int, int]
		self._set_attr(keys, match, types)
	def _parse_type(self, line):
		patten = r'Type:\s(\w+)\s+Anticodon:\s+(\S+)\sat\s(\d+)\-(\d+)\s+\((\d+)\-(\d+)\)\s+Score:\s+(\d+\.?\d*)'
		match = re.compile(patten).match(line).groups()
		keys = ['type', 'anti_codon', 'codon_start', 'codon_end', 'codon_xstart', 'codon_xend', 'score']
		types = [str, str, int, int, int, int, float]
		self._set_attr(keys, match, types)
		self.anti_codon = self.set_anti_codon()
	def _parse_seq(self, line):
		seq = line.strip().split(' ')[1]
		self._set_attr(['seq'], [seq], [str])
	def _parse_struct(self, line):
		struct = line.strip().split(' ')[1]
		self._set_attr(['struct'], [struct], [str])
	@property
	def rna_seq(self):
		return self.seq.upper().replace('T', 'U')
	@property
	def rna_struct(self):
		return self.struct.replace('>', '(').replace('<', ')')
		
