import sys
import re
from Database import Database

class GenbankSummary:
	def __init__(self, summary):
		self.summary = summary
	def __iter__(self):
		return self._parse()
	def _parse(self):
		lines = []
		for line in open(self.summary):
			line = line.strip()
			if not line and lines:
				yield GenbankSummaryLines(lines)
				lines = []
			if line:
				lines += [line]
		if lines:
			yield GenbankSummaryLines(lines)
	def select_by_species(self, fout=sys.stdout):
		d_records = {}
		for rc in self:
		#	print >>sys.stderr, rc.__dict__
			try: d_records[rc.species] += [rc]
			except KeyError: d_records[rc.species] = [rc]
		print >>sys.stderr, d_records.keys()
		for species, records in sorted(d_records.items()):
			refseq_records = self.get_refseq(records)
			if refseq_records:
				records = refseq_records
			record = self.get_longest(records)
			species = species.replace(' ', '_')
			line = [species, record.accession, record.length, record.topology]
			tax = record.lineage
			line += tax
			line = map(str, line)
			print >>fout, '\t'.join(line)
	def get_refseq(self, records):
		return [rc for rc in records if rc.is_refseq]
	def get_longest(self, records):
		return max(records, key=lambda x:x.length)
class GenbankSummaryLines:
	def __init__(self, lines):
		self.lines = lines
		self._parse()
	def _parse(self):
		title = self.lines[0]
		self.number, self.title = \
				re.compile(r'(\d+)\.\s+(.*)').match(title).groups()
		self.species = self._parse_species()
		topology = self.lines[1]
		self.length, self.topology = \
			re.compile(r'([,\d]+)\s+bp\s+(\w+)\s+\S+').match(topology).groups()
		self.length = int(self.length.replace(',', ''))
		accession = self.lines[2]
		self.id = self.accession = \
			re.compile(r'(\S+)\s+\S+').match(accession).groups()[0]
	def _parse_species(self):
		title = self.title.split()
		s,e = 0,2
		if title[0].startswith('UNVERIFIED'):
			title = title[1:]
		if title[1] == 'x':
			e = 3
		if title[1] == 'sp.':
			e = 1
		if title[1] == 'aff.':
			title.pop(1)
		if title[1] == 'cf.':
			title.pop(1)
		return ' '.join(title[s:e])

	@property
	def is_refseq(self):
		if re.compile(r'N[A-Z]_\d+').match(self.accession):
			return True
		return False
	@property
	def lineage(self):
		try:
			return Database().get_taxonomy(self.species)
		except ValueError as e:
			print >>sys.stderr, e
			return []

def main():
	summary, outtsv = sys.argv[1], sys.stdout
	GenbankSummary(summary).select_by_species(outtsv)

if __name__ == '__main__':
	main()
