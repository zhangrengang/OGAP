# coding: utf8
import sys
import re
import copy
from collections import Counter, OrderedDict, namedtuple
import networkx as nx
from Bio.Data import CodonTable
from Bio.Seq import Seq
from lazy_property import LazyWritableProperty as lazyproperty
#try: from small_tools import open_file as open
try: from xopen import xopen as open
except ImportError: pass
try: from RunCmdsMP import logger
except ImportError: pass
try: from translate_seq import translate_seq as translate_cds
except ImportError: pass
try: from Region import Regions, Region, Position
except ImportError: pass

class GffLine(object):
	'''parse a line of standard gff'''
	def __init__(self, line, **kargs):
		if isinstance(line, (list, tuple)):
			pass
		elif isinstance(line, str):
			line = line.strip().split('\t')
		else:
			raise TypeError('Type {} of `{}` is invalid for gff line'.format(type(line), line))
		if not len(line) == 9:
			raise ValueError('Length ({}) `{}` of gff line is not equal to 9'.format(len(line), line))
		(self.chrom, self.source, self.type, self.start, self.end, 
			 self.score, self.strand, self.frame, self.raw_attributes) = line
		self.line = tuple(line)
		self.start, self.end = int(self.start), int(self.end)
		try: self.score = float(self.score)
		except: pass
		try: self.frame = int(self.frame)
		except: pass
		self.attributes_ = self.raw_attributes

	def __hash__(self):
		return hash(self.key)
	def __eq__(self, other):
		return self.key == other.key
	def __len__(self):
		return self.end - self.start + 1
	def __str__(self):
		#print self.attributes_
		line = [self.chrom, self.source, self.type, self.start, self.end, 
				self.score, self.strand, self.frame, self.attributes_]
		line = map(str, line)
		return '\t'.join(line)
	@lazyproperty
	def key(self):
		return self.line #(self.chrom, self.type, self.start, self.end, self.strand)
	@lazyproperty
	def coord(self):
		return '{chrom}:{start}-{end}({strand})'.format(**self.__dict__)
	@lazyproperty
	def attributes(self):
		if isinstance(self.raw_attributes, OrderedDict):
			return self.raw_attributes
		else:
			return self.parse_attr(self.raw_attributes)
	@lazyproperty
	def id(self):
		return self.ID
	@lazyproperty
	def parent(self):
		return self.Parent
	@lazyproperty
	def name(self):
		return self.Name
	@lazyproperty
	def ID(self):
		try: 
			return self.attributes['ID']
		except KeyError:
			return None
	@lazyproperty
	def Parent(self):
		try: 
			return self.attributes['Parent']
		except KeyError:
			return None
	@lazyproperty
	def Name(self):
		try:
			return self.attributes['Name']
		except KeyError:
			return None
	@property
	def keys(self):
		return self._keys_
	@lazyproperty
	def _keys_(self):
		return ['ID', 'Parent', 'Name']
	def parse_attr(self, attributes):
		#attributes = attributes.rstrip(';')
		#return OrderedDict(kv.split('=') for kv in attributes.split(';') if kv)
		d_attr = OrderedDict()
		for kv in attributes.split(';'):
			skv = kv.split('=', 1)
			if not len(skv) ==2:
				continue
			k, v = skv
			k, v = k.strip(), v.strip()
			d_attr[k] = v
		return d_attr
	def update_attr(self):
		for attr_key in self.keys:
			try:
				attr_value = getattr(self, attr_key)
				self.attributes[attr_key] = attr_value
			except: pass
	def format_attr(self):
		return self._format_attr()
	def _format_attr(self):
		return ';'.join(['{}={}'.format(k,v) for k, v in self.attributes.items() if v])
	def _write(self, fout=sys.stdout):
		#line = [self.chrom, self.source, self.type, self.start, self.end, 
		#		self.score, self.strand, self.frame, self._attributes]
		#line = map(str, line)
		#print self.attributes_
		print >>fout, str(self) #'\t'.join(line)
	def write(self, fout):
		self.update_attr()
		self.attributes_ = self.format_attr()
		return self._write(fout)
	@property
	def region(self):
		return Region(chrom=self.chrom, start=self.start, end=self.end)
	@lazyproperty
	def positions(self):	# list
		if self.strand == '-':
			return list(reversed(self.region.positions))
		return self.region.positions
	def extract_seq(self, chrom):
		seq = chrom[self.start-1:self.end]
		if self.strand == '-':
			seq = seq.reverse_complement()
		return seq
	def has_overlap(self, other):
		if not self.chrom == other.chrom:
			return False
		return max(0, min(self.end, other.end) - max(self.start, other.start))
	def overlaps(self, other):
		return self.region.overlaps(other.region)



class GtfLine(GffLine):
	'''standard gtf'''
	def __init__(self, line, **kargs):
		super(GtfLine, self).__init__(line)
	def parse_attr(self, attributes):
		return self._parse_attr(attributes)
	def _parse_attr(self, attributes):
		d_attr = OrderedDict()
		for kv in attributes.strip(';').split('; '):
			skv = kv.split(' ', 1)
			if not len(skv) ==2:
				continue
			k, v = skv
			k, v = k.strip(), v.strip('"').strip()
			d_attr[k] = v
		return d_attr
#		return OrderedDict(re.compile(r'(\S+) "?(.*?)"?[;$]').findall(attributes))
	def format_attr(self):
		return ' '.join(['{} "{}";'.format(k,v) for k, v in self.attributes.items()])
	@lazyproperty
	def keys(self):
		return ['gene_id', 'transcript_id']
	@lazyproperty
	def gene_id(self):
		return self.attributes['gene_id']
	@lazyproperty
	def transcript_id(self):
		return self.attributes['transcript_id']
	def to_gff(self):
		attributes = copy.deepcopy(self.attributes)
		keys = self.keys
		if self.type == 'gene':
			self.ID = self.gene_id
		elif self.type == 'transcript' or self.type.endswith('RNA'):
			self.ID = self.transcript_id
			self.Parent = self.gene_id
		else:
			self.ID = '{}:{}'.format(self.transcript_id, self.type.lower())
			self.Parent = self.transcript_id
		self.attributes = OrderedDict()
		self.keys = self._keys_
		self.update_attr()
		self.attributes.update(attributes)
		for key in keys:
			self.attributes[key] = None
		self.attributes_ = self._format_attr()
	def write_gff(self, fout):
		self.to_gff()
		#print self.__dict__, self.attributes_
		return self._write(fout)
	#def _write(self, fout=sys.stdout):
	#	print >>fout, str(self)

class AugustusGtfLine(GtfLine):
	'''augustus gff2/gtf output is not standard'''
	def __init__(self, line, **kargs):
		super(AugustusGtfLine, self).__init__(line)
	def parse_attr(self, attributes):
		if self.type == 'gene':
			gene_id = attributes
			return OrderedDict(gene_id=gene_id, ID=gene_id)
		elif self.type == 'transcript':
			transcript_id = attributes
			gene_id = transcript_id.split('.')[0]
			return OrderedDict(gene_id=gene_id, transcript_id=transcript_id, ID=transcript_id, Parent=gene_id)
		elif self.source != 'AUGUSTUS':
			return OrderedDict()
		else:
			return self._parse_attr(attributes)
			
class ExonerateGtfLine(GtfLine):
	'''exonerate gff output is not standard'''
	def __init__(self, line, **kargs):
		super(ExonerateGtfLine, self).__init__(line)

	def parse_attr(self, attributes):
		return OrderedDict(kv.split(' ') for kv in attributes.split(' ; '))

global HAS_PRINTED
global PRINT_LIMIT
HAS_PRINTED = 0
PRINT_LIMIT = 10

class GffLines(object):
	'''gff parser for lines'''
	def __init__(self, gff, parser=None, format=None, ):
		self.gff = gff
		if format is None and parser is None:
			if isinstance(self.gff, str):
				self.format = self.gff.split('.')[-1]	# gff3/gff or gtf
			else:
				self.format = 'gff3'
		else:
			self.format = format
		self.parser = parser
		if self.parser is None:
			if self.format in {'gff3', 'gff'}:
				self.parser = GffLine
			elif self.format in {'gtf'}:
				self.parser = GtfLine
			else:
				self.parser = GffLine	# default
				#raise ValueError('Unknown format: `{}` out of gff3 (gff3/gff) and gtf'.format(self.format))
		self.annotations = []
	def __iter__(self):
		return self._parse()
	def _parse(self):
		HAS_PRINTED = 0
		PRINT_LIMIT = 2
		if isinstance(self.gff, file):
			handle = self.gff
		elif isinstance(self.gff, str):
			handle = open(self.gff)
		else:
			raise TypeError('Type {} of `{}` is invalid for gff'.format(type(self.gff), self.gff))
		for line in handle:
			#print >>sys.stderr, line
			if line.startswith('#'):
				self.annotations += [line]
				continue
			line = line.strip().split('\t')
			if not line:
				continue
			if len(line) != 9 :
				if HAS_PRINTED < PRINT_LIMIT:
					print >>sys.stderr, '[WARN] length of {} is not 9'.format(line)
					HAS_PRINTED += 1
				continue
			try: yield self.parser(line)
			except ValueError:
				print >>sys.stderr, '[WARN] LINE {} can not parsed'.format(line)

class GtfLines(GffLines):
	def __init__(self, gff, parser=GtfLine):
		super(GtfLines, self).__init__(gff, parser)
class AugustusGtfLines(GtfLines):
	def __init__(self, gff, parser=AugustusGtfLine):
		super(AugustusGtfLines, self).__init__(gff, parser)

	def __iter__(self):
		return self.parse()
	def parse(self):
		for line in self._parse():
			if not line.source == 'AUGUSTUS':
				continue
			yield line

class ExonerateGtfLines(GffLines):
	def __init__(self, gff, parser=ExonerateGtfLine):
		super(ExonerateGtfLines, self).__init__(gff, parser)
			
class GffGenes(object):
	'''gff parser for genes'''
	def __init__(self, gff, parser=GffLines, **kargs):
		self.gff = gff
		self.parser = parser
		self.kargs = kargs
	def __iter__(self):
		return self._parse()
	def _parse(self):
		record = GffRecord(**self.kargs)
		i = 0
		j = 0
		ids = {} #set([])
		has_gene = False
		for line in self.parser(self.gff):
			i += 1
			j += 1
			id, parent = line.id, line.parent
#			if i < 10:
#				print >> sys.stderr, line
			if line.type == 'gene':
				has_gene = True
			if id in ids:	# duplicated ID of CDS and UTR
				lasti = ids[id]
				id = '{}.{}'.format(id, i)
				if line.type == 'gene' and i - lasti < 10:	 # trans-splicing in RefSeq
					continue
			ids[id] = j #.add(id)
			if (parent is None or line.type == 'gene' or (parent not in record and not has_gene) ) \
					and len(record.nodes()) > 0:
				yield record
				record = GffRecord()
				i = 1
				record.add_node(id, line=line, index=i)
			else:
				record.add_node(id, line=line, index=i)
				if parent is not None and line.type != 'gene':	# gene: parent: gene_id
					record.add_edge(parent, id)
		if len(record.nodes()) > 0:
			yield record
class PasaGffGenes(GffGenes):
	def __init__(self, gff, parser=GffLines):
		self.gff = gff
		self.parser = parser
	def _parse(self):
		record = GffRecord()
		ids = {}
		i, j = 0, 0
		for line in self.parser(self.gff):
			i += 1
			j += 1
			id, parent = line.id, line.parent
			if id in ids:
				id = '{}.{}'.format(id, i)
			ids[id] = i
			if line.type == 'gene' and len(record.nodes()) > 0:
				yield record
				record = GffRecord()
				i = 1
				ids = {}
				record.add_node(id, line=line, index=i)
			else:
				record.add_node(id, line=line, index=i)
				if line.type != 'gene':
					record.add_edge(parent, id)
		if len(record.nodes()) > 0:
			yield record
class GtfGenes(GffGenes):
	def __init__(self, gff, parser=GtfLines):
		self.gff = gff
		self.parser = parser
		self.annotations = []
	def _parse(self):
		lines = self.parser(self.gff)
		record = GffRecord()
		for i, line in enumerate(lines):
			if line.type == 'gene':
				parent = None
				try: id = line.gene_id
				except: id=None
			elif line.type == 'transcript' or line.type.endswith('RNA'):
		#		print line.attributes
				id, parent = line.transcript_id, line.gene_id
			else:
				id, parent = i, line.transcript_id
			if parent is None and len(record.nodes()) > 0:
				self.annotations = lines.annotations
				lines.annotations = []
				yield record
				record = GffRecord()
				record.add_node(id, line=line, index=i)
			else:
				record.add_node(id, line=line, index=i)
				if parent is not None:
					record.add_edge(parent, id)
		self.annotations = lines.annotations
		if len(record.nodes()) > 0:
			yield record
class StringTieGtfGenes(GtfGenes):
	def __init__(self, gff, parser=GtfLines):
		super(StringTieGtfGenes, self).__init__(gff, parser)
	def _parse(self):
		lines = self.parser(self.gff)
		record = GffRecord()
		last_gene = ''
		for i, line in enumerate(lines):
			if line.type == 'transcript':
				id, parent = line.transcript_id, line.gene_id
				gene = parent
			else:
				id, parent = i, line.transcript_id
			if last_gene and last_gene != gene:
				yield record
				record = GffRecord()
				record.add_node(id, line=line, index=i)
			else:
				record.add_node(id, line=line, index=i)
				record.add_edge(parent, id)
			last_gene = gene
		yield record

class PinfishGtfGenes(GtfGenes):
	def __init__(self, gff, parser=GtfLines):
		super(PinfishGtfGenes, self).__init__(gff, parser)
	def _parse(self):
		lines = self.parser(self.gff)
		record = GffRecord()
		for i, line in enumerate(lines):
			if line.type == 'mRNA':
				id, parent = line.transcript_id, None
				if len(record.nodes()) > 0:
					yield record
					record = GffRecord()
			else:
				id, parent = i, line.transcript_id
			record.add_node(id, line=line, index=i)
			if parent is not None:
				record.add_edge(parent, id)
		yield record
class AugustusGtfGenes(GtfGenes):
	def __init__(self, gff, parser=AugustusGtfLines):
		super(AugustusGtfGenes, self).__init__(gff, parser)
	def __iter__(self):
		return self.parse()
	def parse(self):
		for record in self._parse():
			record.annotations = AugustusGtfAnnotations(self.annotations)
			types = {line.type for line in record.lines}
			record.is_complete = False
			if 'start_codon' in types and 'stop_codon' in types:
				record.is_complete = True
			yield record
			
class AugustusGtfAnnotations():
	def __init__(self, lines):
		self._parse(lines)
	def _parse(self, lines):
		protein = False
		exons = False
		obeyed = False
		self.protein_sequence = []
		self.blocks = []
		self.supported, self.total_exons = 0, 0
		self.supported_P = 0
		self.supported_E = 0
		self.fully_obeyed = 0
		self.obeyed_P = 0
		self.obeyed_E = 0
		for line in lines:
			line = line.strip('# \n')
			if line.startswith('protein sequence'):
				self.protein_sequence += [line.split(' ')[-1][1:].strip(']')]
				if not line.endswith(']'):
					protein = True
			elif protein:
				self.protein_sequence += [line.strip(']')]
				if line.endswith(']'):
					protein = False
			if line.startswith('sequence of block'):
				self.blocks += [re.compile(r'(\d+)\s+\[([A-Z]+)\]\s+(\d+)').search(line).groups()]	
			if line.startswith('CDS exons'):
				self.supported, self.total_exons = map(int, line.split(' ')[-1].split('/'))
				exons = True
			elif exons:
				if line.startswith('P:'):
					self.supported_P = int(re.compile(r'^P:\s+(\d+)').search(line).groups()[0])
				elif line.startswith('E:'):
					self.supported_E = int(re.compile(r'^E:\s+(\d+)').search(line).groups()[0])
				else:
					exons = False
			if line.startswith('hint groups fully obeyed'):
				self.fully_obeyed = int(line.split(' ')[-1])
				obeyed = True
			elif obeyed:
				if line.startswith('P:'):
					self.obeyed_P = int(re.compile(r'^P:\s+(\d+)').search(line).groups()[0])
				elif line.startswith('E:'):
					self.obeyed_E = int(re.compile(r'^E:\s+(\d+)').search(line).groups()[0])
				else:
					obeyed = False
		self.protein_sequence = ''.join(self.protein_sequence)
class GffExons(object):
	def __init__(self, exons):
		self.exons = exons
		self.rna_type = 'mRNA'
	def __iter__(self):
		return iter(self.exons)
	def __len__(self):
		return len(self.exons)
	def __str__(self):
		return ' -> '.join(['{chrom}:{start}-{end}({strand})'.format(**exon.__dict__) \
								for exon in self if exon.type=='exon'])
	def count_exon(self):
		return len([1 for exon in self if exon.type=='exon'])
	def __getitem__(self, index):
		if isinstance(index, int):
			return self.exons[index]
		else:
			return self.__class__(path=self.exons[index])
	def __hash__(self):
		return hash(str(self))
	def __add__(self, other):
		exons = self.exons + other.exons
		return self.__class__(exons=exons)
	def overlaps(self, other):
		self_genes = self.filter('gene')
		other_genes = other.filter('gene')
		if not (self_genes and other_genes):
			return self.region.overlaps(other.region)
		has_overlaps = [self_gene.has_overlap(other_gene) \
					for self_gene, other_gene in zip(self_genes, other_genes)]
		if len(has_overlaps) == 0:
			return False
		return all(has_overlaps)
	
	@property
	def total_length(self):
		length = 0
		for exon in self.exons:
			length += exon.end - exon.start + 1
		return length
	def get_region(self):
		starts  = [exon.start  for exon in self]
		ends	= [exon.end	for exon in self]
		strands = [exon.strand for exon in self]
		chroms  = [exon.chrom  for exon in self]
		chrom  = chroms[0]
		strand = strands[0]
		start  = min(starts + ends)
		end	= max(starts + ends)
		return chrom, start, end, strand
	@lazyproperty
	def coord(self):
		return self.get_region()
	@lazyproperty
	def chroms(self):
		return {exon.chrom  for exon in self}
	@lazyproperty
	def region(self):
		chrom, start, end = self.coord[:3]
		return Region(chrom, start, end)
	@lazyproperty
	def positions(self):
		positions = []
		for line in self:
			positions += line.region.positions
		return positions
	@lazyproperty
	def chrom(self):
		return self.coord[0]
	@lazyproperty
	def start(self):
		return self.coord[1]
	@lazyproperty
	def end(self):
		return self.coord[2]
	@lazyproperty
	def strand(self):
		return self.coord[3]
	def extract_seq(self, d_seqs, type='exon'):
		'''exons must be in order'''
		seqs = []
		for line in self:	# exon
			if not line.type == type:
				continue
			if isinstance(d_seqs, dict):
				seq = d_seqs[line.chrom]	# support for extracting from different contigs
			else:
				seq = d_seqs
			exon_seq = seq[line.start-1:line.end]
			if isinstance(exon_seq, list):
				exon_seq = Seq(''.join(exon_seq))
			if line.strand == '-':
				exon_seq = exon_seq.reverse_complement()
			seqs += [str(exon_seq)]
		return ''.join(seqs)
	def to_positions(self, type='exon'):
		positions = []
		for line in self:   # exon
			if not line.type == type:
				continue
		#   if line.strand == '-':
		#	   positions += list(reversed(line.positions))
		#   else:
			positions += line.positions	# reverse if strand == '-'
		return positions

	def to_gff_record(self):
		pass
	def extend_gene(self, gene, gene_copy, source=None, rna_type='mRNA', pseudo=False):
		rna_type = 'mRNA' if rna_type == 'CDS' else rna_type
		if source is None:
			source = gene_copy.source
		gene_name = gene.name
		product = gene.product
		gene_id = gene_copy.id
		rna_id = '{}.t1'.format(gene_id)
		type = 'gene'
		score, frame = '.', '.'
		trans_splicing = False
		if len(gene_copy) > 1:
			trans_splicing = True
		chrom, start0, end0, strand0 = self.get_region()
		
		attributes0 = OrderedDict()
		for k, v in zip(('ID', 'Parent', 'Name', 'gene', 'product',
							 ),
						(gene_id, None, gene_name, gene_name, product,
							 )):
			attributes0[k] = v
		if trans_splicing:
			attributes0['exception'] = 'trans-splicing'
		if pseudo:
			attributes0['pseudo'] = 'true'
		lines = []
		for i, part in enumerate(gene_copy):
			_part = '{}/{}'.format(i+1, len(gene_copy))
			# gene
			chrom, start, end, strand = (part.seqid, part.start, 
						part.end, part.strand)
			if len(gene_copy) > 1:
				attributes0['part'] = _part
			_attributes0 = copy.deepcopy(attributes0)
			line = [chrom, source, type, start, end, score, strand, frame, _attributes0]
			line = GffLine(line)
			lines += [line]
			# RNA
			line = copy.deepcopy(line)
			attributes = copy.deepcopy(attributes0)
			attributes.update(ID=rna_id, Parent=gene_id)
			line.attributes = attributes
			line.type = rna_type
			lines += [line]
		# exons
		attributes = copy.deepcopy(attributes0)
		for key in ['Name', 'gene']:
			attributes[key] = None
		if len(gene_copy) > 1:
			attributes['part'] = None
		attributes.update(Parent=rna_id)
		for i, exon in enumerate(self):
			# CDS, etc.
			attributes = copy.deepcopy(attributes)
			line = copy.deepcopy(line)
			line.__dict__.update(**exon.__dict__)
			line.ID = '{}:{}{}'.format(rna_id, exon.type[0].lower(), i+1)
#			attributes.update(ID='{}:{}{}'.format(rna_id, exon.type[0].lower(), i+1))
			line.attributes = attributes
			line.Parent = rna_id

			line.__dict__.update(source=source) #, score=score)
			#print >>sys.stderr, line.attributes
			# exon
			if exon.type != 'exon':
				type = 'exon'
				exon_line = copy.deepcopy(line)
				exon_line.ID = '{}:{}{}'.format(rna_id, type[0], i+1)
#				exon_line.attributes.update(ID='{}:{}{}'.format(rna_id, type[0], i+1))
				exon_line.__dict__.update(type=type, frame=frame, score=score)
				lines += [exon_line]
			lines += [line]
		record = GffExons(lines)
		record.source = source
		record.chrom = chrom
		record.start = start0
		record.end = end0
		record.strand = strand0
		record.gene = record.name = gene_name
		record.id = gene.id
		record.product = product
		record.gene_id = gene_id
		record.rna_id = rna_id
		record.rna_type = rna_type
		record.trans_splicing = trans_splicing
		return record
	def to_tbl(self, fout, chrom=None, feat_type='gene', wgs=True, gene=True,
			transl_table=1, rna_type=None, locus_tag=None, note=None):
		if self.rna_type == 'repeat':
			exon = self[0]
			start, end = exon.start, exon.end
			line = [start, end]		# 1-based
			line += ['repeat_region']
			line = map(str, line)
			print >>fout, '\t'.join(line)
			try:
				line = ['', '', 'note', exon.id]
				print >>fout, '\t'.join(line)
			except AttributeError:
				pass
			try:
				line = ['', '', 'rpt_type', exon.attributes['rpt_type']]
				print >>fout, '\t'.join(line)
			except AttributeError:
				pass
			return None
		# if feat_type == 'gene': # otherwise repeat etc.
		# gene
		if rna_type is None:
			rna_type = 'CDS' if self.rna_type == 'mRNA' else self.rna_type
		if rna_type == 'CDS':
			exons = self.filter(rna_type)
		else:
			exons = self.filter('exon')
		if chrom is not None:
			exons.exons = [exon for exon in exons if exon.chrom == chrom]
		genes = self.filter('gene')
		if chrom is not None:
			genes.exons = [exon for exon in genes if exon.chrom == chrom]
		for i, exon in enumerate(genes):
			start, end = exon.start, exon.end
			start, end = (end, start) if exon.strand == '-' else (start, end)
			line = [start, end]		# 1-based
			if i == 0:
				line += ['gene']
			line = map(str, line)
			print >>fout, '\t'.join(line)
		try:
			line = ['', '', 'gene', self.gene]
			print >>fout, '\t'.join(line)
		except AttributeError:
			pass
		try:
			if self.trans_splicing:
				line = ['', '', 'exception', 'trans-splicing']
				print >>fout, '\t'.join(line)
		except AttributeError:	# no gene name
			pass
		if locus_tag is not None:
			line = ['', '', 'locus_tag', locus_tag]
			print >>fout, '\t'.join(line)

		if note:
			line = ['', '', 'note', note]
			line = map(str, line)
			print >>fout, '\t'.join(line)

		if getattr(self, 'pseudo', None):
			line = ['', '', 'pseudo']
			print >>fout, '\t'.join(line)
			return

		# shared lines by mRNA and CDS
		lines = []
		try:
			product = self.product
		except AttributeError:  # no product
			if rna_type == 'CDS':
				product = 'hypothetical protein'
			else:
				product = rna_type
		lines += [ ['', '', 'product', product] ]
		if locus_tag is not None:
			transcript_id = '{}|{}|{}'.format('gnl', locus_tag.split('_')[0], self.rna_id)
			if wgs: # nesscery to submit to WGS, but not for GenBank
				lines += [ ['', '', 'transcript_id', transcript_id] ]
			if rna_type == 'CDS':
				lines += [ ['', '', 'protein_id', transcript_id+':cds'] ]
		# mRNA
		if locus_tag is not None and rna_type == 'CDS' and wgs:
			for i, exon in enumerate(exons):
				start, end = exon.start, exon.end
				start, end = (end, start) if exon.strand == '-' else (start, end)
				line = [start, end]
				if i == 0:
					line += ['mRNA']
				line = map(str, line)
				print >>fout, '\t'.join(line)
			for line in lines:
				print >>fout, '\t'.join(line)
		# CDS		
		for i, exon in enumerate(exons):
			start, end = exon.start, exon.end
			start, end = (end, start) if exon.strand == '-' else (start, end)
			line = [start, end]
			if i == 0:
				line += [rna_type]
				if rna_type == 'CDS':
					#print >>sys.stderr, vars(exon)
					codon_start = exon.frame + 1
			line = map(str, line)
			print >>fout, '\t'.join(line)

		for line in lines:
			print >>fout, '\t'.join(line)
		if rna_type == 'CDS' and transl_table != 1:
			line = ['', '', 'transl_table', transl_table]
			line = map(str, line)
			print >>fout, '\t'.join(line)
		if rna_type == 'CDS':
			line = ['', '', 'codon_start', codon_start]
			line = map(str, line)
			print >>fout, '\t'.join(line)

	def reverse(self, length):
		exons = []
		for exon in self:
			exon.start, exon.end = length - exon.end + 1, length - exon.start + 1
			if exon.strand == '-':
				exon.strand = '+'
			else:
				exon.strand = '-'
			exons += [exon]
		return GffExons(exons)
	def filter(self, type='exon'):
		if isinstance(type, str):
			types = {type}
		else:
			types = set(type)
		return GffExons([line for line in self if line.type in types])
	def translate_cds(self, seq, **kargs):
		return translate_cds(Seq(seq), **kargs)
	def write(self, fout):
		for exon in self:
			#print >>sys.stderr, exon.__dict__
			exon.write(fout)
	def update_frame(self):
		'''add frame'''
		accumulated_length = 0
		for exon in self:
			length = exon.end - exon.start + 1
			_frame = accumulated_length % 3
			frame = 3 - _frame if _frame > 0 else _frame
			exon.frame = frame
			accumulated_length += length

class GtfExons(GffExons):
	def __init__(self, exons):
		super(GtfExons, self).__init__(exons)
	def write_gff(self, fout):
		for exon in self:
			#print >>sys.stderr, exon.__dict__, exon.keys
			#exon.keys = ['ID', 'Parent', 'Name']
			exon.write_gff(fout)
class PinfishGtfExons(GtfExons):	# 202-4-29 for shuiqingshu add_mRNA
	def __init__(self, exons):
		super(PinfishGtfExons, self).__init__(exons)
	def extend_cds(self):
		'''only work when mapping with cds sequences'''
		exons = self.filter('exon')
		exons = copy.deepcopy(exons)
		exons.update_frame()
		for exon in exons:
			exon.type = 'CDS'
		return self + exons
class ExonerateGtfExons(GtfExons):
	def __init__(self, exons):
		super(ExonerateGtfExons, self).__init__(exons)
	@property
	def is_compelte(self):
		if self.start_codon_found and self.stop_codon_found:
			return True
		return False
	def get_seq(self, seq, start, end):
		seq = seq[start:end]
		if isinstance(seq, list):
			seq = Seq(''.join(seq))
		return seq
	def trace_start_codon(self, seq, transl_table=1, max_dist=150):
		first_exon = self.exons[0]
		strand = first_exon.strand
		# trace start codon: min distance
		tss = first_exon.end if strand == '-' else first_exon.start-1
		# --> right
		upper = len(seq) if strand == '-' else first_exon.end
		#print >>sys.stderr, tss, upper, strand, first_exon.start, first_exon.end
		i = tss
		for i in range(tss, upper, 3):
			codon = self.get_seq(seq, i, i+3) #seq[i:i+3]
			if self.is_start_codon(codon, strand, transl_table):
				dist1 = i - tss
				break
			if strand == '-' and self.is_stop_codon(codon, strand, transl_table):
				dist1 = max_dist + 1
				break
		else:
			dist1 = max_dist + 1
		if strand == '+' and i >= first_exon.end:	# exceed exon
			dist1 = max_dist + 1
		# <-- left
		lower = first_exon.start if strand == '-' else 0
		j = tss
		for j in range(tss, lower, -3):
			codon = self.get_seq(seq, j-3,j) #seq[j-3:j]
			if self.is_start_codon(codon, strand, transl_table):
				dist2 = tss - j + 3
				break
			if strand != '-' and self.is_stop_codon(codon, strand, transl_table):
				dist2 = max_dist + 1
				break
		else:
			dist2 = max_dist + 1
		if strand == '-' and i <= first_exon.start:	# exceed exon
			dist2 = max_dist + 1
		start_codon_found = True
		if min(dist1, dist2) > max_dist:
			start_codon_found = False
		if start_codon_found:
			if dist1 < dist2:
				if strand == '-':
					first_exon.end = i + 3
				else:
					first_exon.start = i + 1
			else:
				if strand == '-':
					first_exon.end = j
				else:
					first_exon.start = j -3 + 1
		self.has_start_codon = start_codon_found
	def trace_stop_codon(self, seq, transl_table=1, max_dist=150):
		last_exon = self.exons[-1]
		strand = last_exon.strand
		# trace stop codon
		tes = last_exon.start-1 if strand == '-' else last_exon.end
		# --> right
		upper = last_exon.end if strand == '-' else len(seq)
		for i in range(tes, upper, 3):
			codon = self.get_seq(seq, i, i+3) #seq[i:i+3]
			if self.is_stop_codon(codon, strand, transl_table):
				dist1 = i - tes
				break
		else:
			dist1 = max_dist + 1
		# <-- left
		lower = 0 if strand == '-' else last_exon.start
		for j in range(tes, lower, -3):
			codon = self.get_seq(seq, j-3, j) #seq[j-3:j]
			if self.is_stop_codon(codon, strand, transl_table):
				dist2 = tes - j + 3
				break
		else:
			dist2 = max_dist + 1
			
		stop_codon_found = True
		if min(dist1, dist2) > max_dist:
			stop_codon_found = False
		if stop_codon_found:
			if dist1 < dist2:	# -->
				if strand == '-':
					last_exon.start = i + 1
				else:
					last_exon.end = i + 3
			else:	# <--
				if strand == '-':
					last_exon.start = j - 2 #3
				else:
					last_exon.end = j
		self.has_stop_codon = stop_codon_found
	def truncate_exons(self, seq, transl_table=1, max_dist=150):
		last_truncated = False
		for exon in self:
			strand = exon.strand
			if strand == '-':
				start, stop, step = exon.end - exon.frame, exon.start+3, -3
				if last_truncated:
					exon.end -= exon.frame
				for i in range(start, stop, step):
					codon = self.get_seq(seq, i-3, i) #seq[i-3:i]
					if self.is_stop_codon(codon, strand, transl_table):
						exon.start = i
						last_truncated = True
						break
				else:
					last_truncated = False
					
			else:	# +
				start, stop, step = exon.start-1+exon.frame, exon.end-3, 3
				if last_truncated:
					exon.start += exon.frame
					exon.frame = 0
				for i in range(start, stop, step):
					codon = self.get_seq(seq, i, i+3) #seq[i:i+3]
					if self.is_stop_codon(codon, strand, transl_table):
						exon.end = i	# 1-based
						last_truncated = True
						break
				else:
					last_truncated = False
		if last_truncated:
			if strand == '-':
				exon.start -= 3
			else:
				exon.end += 3
	def is_stop_codon(self, codon, strand, transl_table=1):
		stop_condons = CodonTable.unambiguous_dna_by_id[transl_table].stop_codons
		return self.is_codon_in_set(codon, strand, stop_condons)
	def is_start_codon(self, codon, strand, transl_table=1):
		start_codons = CodonTable.unambiguous_dna_by_id[transl_table].start_codons
		return self.is_codon_in_set(codon, strand, start_codons)
	def is_codon_in_set(self, codon, strand, codons):
		if strand == '-':
			codon = codon.reverse_complement()
		if str(codon.upper()) in set(codons):
			return True
		return False
	
	def to_gff(self, fout):
		first_exon, last_exon = self.exons[0], self.exons[-1]
		chrom, start, end, strand = (first_exon.chrom, first_exon.start, 
									last_exon.end, first_exon.strand)
		source= 'Exonerate'
		gene_id = '{}_{}-{}'.format(chrom, start, end)
		transcript_id = gene_id + '.1'
		score, frame = self.score, '.'
		gene_attr = 'ID={};Target={}'.format(gene_id, self.hit)
		tran_attr = 'ID={};Parent={}'.format(transcript_id,gene_id)
		for type, attribute in zip(['gene', 'mRNA'], [gene_attr, tran_attr]):
			line = [chrom, source, type, start, end, score, strand, frame, attribute]
			line = map(str, line)
			print >>fout, '\t'.join(line)
		for exon in self:	# no intron
			cds = copy.deepcopy(exon)
			cds.type = 'CDS'
			cds.source = source
			cds.score = '.'
			cds.gene_id = gene_id
			cds.transcript_id = transcript_id
			cds.write_gff(fout)		# gtf
			
	def to_augustus_gtf(self, fout, index=1):
		first_exon, last_exon = self.exons[0], self.exons[-1]
		chrom, start, end, strand = (first_exon.chrom, first_exon.start, 
									last_exon.end, first_exon.strand)
		source= 'AUGUSTUS'
		gene_id = 'g{}'.format(index)
		transcript_id = gene_id + '.t1'
		score, frame = '.', '.'
		for type, attribute in zip(['gene', 'transcript'], [gene_id, transcript_id]):
			line = [chrom, source, type, start, end, score, strand, frame, attribute]
			line = map(str, line)
			print >>fout, '\t'.join(line)
		if self.has_start_codon:
			if strand == '-':
				start_codon = copy.deepcopy(last_exon)
				start_codon.start = start_codon.end - 2
			else:
				start_codon = copy.deepcopy(first_exon)
				start_codon.end = start_codon.start + 2
			start_codon.source = source
			start_codon.type = 'start_codon'
			start_codon.gene_id = gene_id
			start_codon.transcript_id = transcript_id
			start_codon.write(fout)
		for exon in self:	# no intron
			cds = copy.deepcopy(exon)
			cds.type = 'CDS'
			cds.source = source
			cds.score = self.score
			cds.gene_id = gene_id
			cds.transcript_id = transcript_id
			cds.attributes['score'] = self.score
			#print >>sys.stderr, cds.keys
			cds.write(fout)		# gtf
		if self.has_stop_codon:
			if strand == '-':
				stop_codon = copy.deepcopy(last_exon)
				stop_codon.end = start_codon.start + 2
			else:
				stop_codon = copy.deepcopy(first_exon)
				stop_codon.start = stop_codon.end - 2
			stop_codon.source = source
			stop_codon.type = 'stop_codon'
			stop_codon.gene_id = gene_id
			stop_codon.transcript_id = transcript_id
			stop_codon.write(fout)
		print >> fout, '# coding sequence = [{}]'.format(self.cds_seq)
		print >> fout, '# protein sequence = [{}]'.format(self.pep_seq)
		print >> fout, '# hit = {}::{}'.format(self.hit, self.score)
		
class ExonerateGffGenes(GffGenes):	# each alignment
	def __init__(self, gff=None, parser=ExonerateGtfLines):
		self.gff = gff
		self.parser = parser
	def _parse(self):
		record = GffRecord()
		record.has_frameshift = False
		last_gene = None
		for i, line in enumerate(self.parser(self.gff)):
			id, parent = i, last_gene
			if line.type == 'gene':
				if len(record.nodes()) == 0:
					record.score = line.score
				last_gene = i
				parent = None
			if line.type == 'exon' and 'frameshifts' in line.attributes:
				record.has_frameshift = True
			if parent is None and len(record.nodes()) > 0:
				yield record
				record = GffRecord()
				record.score = line.score
				record.has_frameshift = False
				record.add_node(id, line=line, index=i)
			else:
				record.add_node(id, line=line, index=i)
				if parent is not None:
					record.add_edge(parent, id)
		if len(record.nodes()) > 0:
			yield record
	def get_best_hit(self):
		records = [record for record in self if not record.has_frameshift]
		return max(records, key=lambda x: x.score)
	def get_gene_gtf(self, d_seqs, fout=None, transl_table=1, max_dist=150):
		#best_hit = self.get_best_hit()
		#seq = d_seqs[best_hit.chrom]
		hits = []
		for i, record in enumerate(self):
			if record.has_frameshift:	# exclude frameshift
				continue
			seq = d_seqs[record.chrom]
			exons = self.fit_structure(record, seq,
						transl_table=transl_table, max_dist=max_dist)
			if fout is not None:
				exons.to_augustus_gtf(fout, index=i)
			hits += [exons]
		return hits
	def to_exons(self, d_seqs=None, fout=None):
		hits = []
		for i, record in enumerate(self):
			exons = []
			for line in record.lines:
				if line.type == 'exon':
					line.score = record.score
					exons += [line]
			exons = ExonerateGtfExons(exons)
			try: exons.score = int(record.score)
			except: exons.score = 0
			try: exons.hit = record.gene.attributes['sequence']
			except : exons.hit = None
			exons.id = str(i)
			if fout is not None:
				seq = d_seqs[record.chrom]
				exons.seq = exons.extract_seq(seq)
				print >>fout, '>{}\n{}'.format(exons.id, exons.seq)
			hits += [exons]
		return hits
	def cluster(self):
		hits = self.to_exons()
		if not hits:
			return []
		hits = sorted(hits, key=lambda x: (x.chrom, x.strand, x.start, x.end))
#		for x in hits:
#			print (x.chrom, x.strand, x.start, x.end)
		clster = [[hits[0]]]
		for hit in hits[1:]:
#			print hit
			if hit.overlaps(clster[-1][-1]):
				clster[-1] += [hit]
			else:
				clster += [[hit]]
		return clster
	def get_best_loci(self):
		cluster = self.cluster()
		for clust in cluster:
			yield max(clust, key=lambda x: x.score)
	def fit_structure(self, record, seq, **kargs):
		#print >>sys.stderr, record.score
		exons = [line for line in record.lines if line.type == 'exon']
		exons = ExonerateGtfExons(exons)
		try: exons.hit = record.gene.attributes['sequence']
		except : exons.hit = None
		try: exons.score = int(record.score)
		except: exons.score = 0
		exons.update_frame()
#		seq = list(seq)
		exons.trace_start_codon(seq, **kargs)
		try:
			assert exons.total_length % 3 == 0
		except AssertionError:
			print >>sys.stderr, exons.total_length, exons.score, exons
			
		exons.trace_stop_codon(seq, **kargs)
		exons.cds_seq = exons.extract_seq(seq)
		exons.pep_seq = exons.translate_cds(exons.cds_seq, **kargs)
		try:
			assert exons.total_length % 3 == 0
		except AssertionError:
			print >>sys.stderr, 'not 3 CDS', exons.total_length, exons.score, exons
			
		if '*' in set(exons.pep_seq):
			exons.truncate_exons(seq, **kargs)
			exons.cds_seq = exons.extract_seq(seq)
			exons.pep_seq = exons.translate_cds(exons.cds_seq, **kargs)
		try:
			assert exons.total_length % 3 == 0
		except AssertionError:
			print >>sys.stderr, 'not 3 CDS', exons.total_length, exons.score, exons
		return exons

	def to_hints(self, fout, src='P', pri=4, source='exonerate', intron_type = 'intronpart'):
		d_type = {'P': ('cds', 'CDSpart'), 'E': ('exon', 'exonpart')}
		target_type, hint_type = d_type[src]
		attributes = OrderedDict(src=src, pri=pri)
		grps = set([])
		for i, record in enumerate(self):
			last_start = None
			for line in record.lines: # lines
				if line.type == 'gene':
					grp = line.attributes['sequence']
					if grp in grps:
						grp = '{}|{}'.format(line.attributes['sequence'], i) # multi-hits
					attributes.update(grp=grp)
					grps.add(grp)
					score = line.score
				if line.type == target_type:
					# intron
					if last_start is not None and intron_type is not None:
						end = line.start - 1
						hint_line = [line.chrom, source, intron_type, last_start, end,
							line.score, line.strand, line.frame, attributes]
						hint_line = GffLine(hint_line)
						hint_line.write(fout)
					# exon
					hint_line = [line.chrom, source, hint_type, line.start, line.end, 
						score, line.strand, line.frame, attributes]
					hint_line = GffLine(hint_line)
					hint_line.write(fout)
					last_start = line.end + 1

					
class GffRecord(nx.DiGraph):
	'''DiGraph: gene-mRNA-exon-CDS'''
	def __init__(self, graph=None, flank=5000):	# graph is a nx.DiGraph object
		super(GffRecord, self).__init__()
		self.flank = flank
	#def __iter__(self):
	#	return iter(self.lines)
	def __str__(self):
		return str(self.region)
	@lazyproperty
	def id(self):
		'''use the init gene id as record id'''
		for node in self.non_empty_nodes:
			predecessors = self.predecessors(node)
			if not predecessors or not (set(predecessors)-self.empty_nodes):
				return node
	# @lazyproperty
	# def source_record(self):	# gene
		# return self.node[self.id]['line']
	# @lazyproperty
	# def start(self):
		# return self.source_record.start
	# @lazyproperty
	# def end(self):
		# return self.source_record.end
	@lazyproperty
	def length(self):
		return self.end - self.start + 1
	def __len__(self):
		return self.end - self.start + 1
	def __hash__(self):
		return hash(tuple(self.lines))
	def __eq__(self, other):
		return hash(self) == hash(other)
	@property
	def lines(self):
		return [self.node[node]['line'] for node in self.sort_nodes()]
	@property
	def rna_type(self):
		for line in self.lines:
			if line.type.find('RNA') >=0:
				return line.type
		for line in self.lines:
			if line.type not in {'gene', 'exon', 'CDS'}:
				return line.type
		for line in self.lines:
			if line.type == 'CDS':
				return 'mRNA'
	def to_exons(self):
		return GffExons(self.lines)
	@property
	def rnas(self):
		return GffRNARecords(self)	
	def sort_nodes(self):
		nodes = [node for node in self.nodes() if 'index' in self.node[node]]	# only lines
		return sorted(nodes, key=lambda x: self.node[x]['index'])
	def recur_remove_node(self, node):
		last_nodes = [node]
		nodes_to_del = [] + last_nodes
		while True:
			successors = []
			for last_node in last_nodes:
				successors += self.successors(last_node)
			if not successors:
				break
			nodes_to_del += successors
			last_nodes = successors
		nodes_to_del = set(nodes_to_del)
		removed_nodes = 0
		for node in nodes_to_del:
			self.remove_node(node)
			removed_nodes += 1
		if len(self.nodes()) == 1: # only gene line
			self.remove_node(self.nodes()[0])
			removed_nodes += 1
		return removed_nodes
	def write(self, fout=sys.stdout):
		for line in self.lines:
			#line = self.get_node_feature(node)
			line.write(fout)
	@lazyproperty
	def type_count(self):
		return Counter([self.get_node_feature(node).type for node in self.nodes()])
	def count_type(self, *types):
		types = set(types)
		return len([node for node in self.lines if node.type in types])
	@lazyproperty
	def nexons(self):
		return self.count_type('exon')
	@property
	def is_coding(self):
		if 'CDS' in self.feature_regions:
			return True
		else:
			return False
	@lazyproperty
	def gene(self):		# GffLine
		return self.get_node_feature(self.id)
	@lazyproperty
	def strand(self):
		return self.gene.strand
	@lazyproperty
	def chrom(self):
		return self.gene.chrom
	@lazyproperty
	def start(self):
		return self.gene.start  #Position(chrom=self.chrom, pos=self.gene.start)
	@lazyproperty
	def end(self):
		return self.gene.end	#Position(chrom=self.chrom, pos=self.gene.end)
	@lazyproperty
	def type(self):
		return self.gene.type
	@lazyproperty
	def region(self):
		return self.gene.region
	def overlaps(self, other):
		return self.region.overlaps(other.region)

	@property
	def sorted_nodes(self):
		return self.sort_features(self.node())

	def sort_features(self, fetures):
		return sorted(fetures, key=lambda x: self.get_node_index(x))
	def get_node_feature(self, node):
		try: return self.node[node]['line']
		except KeyError as e:
			print >>sys.stderr, 'get_node_feature-KeyError', node, self.nodes()
			raise KeyError(e)
	@lazyproperty
	def empty_nodes(self):
		return {node for node in self.nodes() if 'line' not in self.node[node]}
	@lazyproperty
	def non_empty_nodes(self):
		return set(self.nodes()) - self.empty_nodes

	def get_node_index(self, node):
		return self.node[node]['index']

	@property
	def streams(self):
		flank = self.flank
		gene_region = self.region
		upstream = Region(chrom=gene_region.chrom, start=gene_region.start-flank-1, end=gene_region.start-1)
		downsteam = Region(chrom=gene_region.chrom, start=gene_region.end+1, end=gene_region.end+1+flank)
		if self.strand == '-':
			upstream, downsteam = downsteam, upstream
		upstream, downsteam = [upstream], [downsteam]
		return upstream, downsteam
	@lazyproperty
	def features(self):
		d_features = {}
		for feature in self.lines:
			#feature = self.get_node_feature(node)
			try: d_features[feature.type] += [feature]
			except KeyError: d_features[feature.type] = [feature]
		return d_features
	@lazyproperty
	def feature_length(self):
		d_length = {}
		for feature in self.lines:
			try: d_length[feature.type] += len(feature)
			except KeyError: d_length[feature.type] = len(feature)
		return d_length
	@lazyproperty
	def feature_count(self):
		return Counter([self.get_node_feature(node).type for node in self.nodes()])
	def share_feat(self, other, feat='exon'):
		self_feat = self.features[feat]
		other_feat = other.features[feat]
		return set(self_feat) & set(other_feat)
	@lazyproperty
	def feature_regions(self):
		d_regions = {}
		for type, features in self.features.items():
			regions = [feature.region for feature in features]
			d_regions[type] = regions
		return d_regions

	def get_longest_rna(self, type='CDS'):
		rnas = GffRNARecords(self)
#		rnas = list(GffRNARecords(self))
#		if len(rnas) == 1:
#			return self
#		print self.id, map(lambda x:x.feature_length, rnas)
		longest = max(rnas, key=lambda x:x.feature_length.get(type, 0))
		return longest
	@property
	def UTRs(self):
		exons, cds = self.feature_regions['exon'], self.feature_regions['CDS']
		utrs = Regions.subtract(exons, cds)
		cds_start = min(cds, key=lambda x: x.start).start
		cds_end = max(cds, key=lambda x: x.end).end
		utr5, utr3 = [], []
		for utr in utrs:
			utr.chrom = self.chrom
			if utr.start < cds_start:
				utr5 += [utr]
			elif utr.end > cds_end:
				utr3 += [utr]
			else:
				raise ValueError('never')
		if self.strand == '-':
			utr5, utr3 = utr3, utr5
		return utr5, utr3
	def get_coords(self, type='exon'):
		feats = self.feature_regions.get(type, [])
		reverse = 1 if self.strand == '-' else 0
		feats = sorted(feats, key=lambda x:x.start, reverse=reverse)
		return ','.join('{}-{}'.format(feat.start, feat.end) for feat in feats)

	def classify_positons(self, d_positon=OrderedDict(), flank=5000, type='exon'):
		'''codon1-3, utr, intron, upstream, downstream'''
		self.flank = flank
		if not self.is_coding:
			return d_positon
		for RNARecord in GffRNARecords(self):
			feat_regions = RNARecord.feature_regions
			utr5, utr3 = RNARecord.UTRs
			upstream, downsteam = RNARecord.streams
			exons = Regions(feat_regions[type])
			introns = exons.get_introns()
			splicing = exons.get_splicing()
			cds = self.features['CDS']
			for postype, regions in zip(
					['upstream', 'utr5', 'intron', 'utr3', 'downstream', 'splicing', 'cds'], 
					[upstream, utr5, introns, utr3, downsteam, splicing, cds]):
				for region in regions:
					for pos in region.positions:
						try: d_positon[pos].add(postype)
						except KeyError: d_positon[pos] = {postype}
			RNARecord.classify_codons(d_positon)
		return d_positon
	def classify_regions(self, flank=3000, type='CDS'):
		self.flank = flank
#		record = self.get_longest_rna()
		utr5, utr3 = self.UTRs
		exons = Regions(self.feature_regions[type])
		introns = exons.get_introns()
		splicing = exons.get_splicing()
		cds = self.feature_regions['CDS']
		region = exons.to_region()
		upstream, downsteam = region.get_streams()
		if self.strand == '-':
			upstream, downsteam =  downsteam, upstream
		for postype, regions in zip(
			['upstream', 'utr5', 'intron', 'utr3', 'downstream', 'splicing', 'cds'],
			[[upstream,], utr5, introns, utr3, [downsteam,], splicing, cds]):
			yield postype, regions
#			print regions
#			for region in regions:
#				yield postype, region

	def classify_codons(self, d_positon=OrderedDict()):
		d_condon = {0:'condon1', 1:'condon2', 2:'condon3'}
		for cds in self.features['CDS']:
			for pos in cds.region.positions:
				if cds.strand == '-':
					condon = (cds.end-pos.pos-cds.frame) % 3
				else:
					condon = (pos.pos-cds.start-cds.frame) % 3
				postype = 'codon' + str(condon+1)
				try: d_positon[pos].add(postype)
				except KeyError: d_positon[pos] = {postype}
		return d_positon
	def classify_fold(self, d_seqs, d_positon=OrderedDict()):
		if not self.is_coding:
			return d_positon
		seq = d_seqs[self.chrom]
		for RNARecord in GffRNARecords(self):
			try: cds = RNARecord.features['CDS']
			except KeyError:
				print >> sys.stderr, 'RNA {} is non-coding in a coding gene'.format(RNARecord.id)
				continue
			reverse = 1 if RNARecord.strand == '-' else 0
			cds = sorted(cds, key=lambda x:x.start, reverse=reverse)
			exons = GffExons(cds)
		#	exons = RNARecord.to_exons()	# sorted
			frame = cds[0].frame
			cds_seq = exons.extract_seq(seq, 'CDS')
			cds_pos = exons.to_positions('CDS')
			#if RNARecord.id == 'Acyan01G0001200.1':
			#   print >>sys.stderr, cds_pos, cds_seq, str(Seq(cds_seq).translate())
			assert len(cds_seq) == len(cds_pos)
			if len(cds_seq) % 3 != 0:
				print >> sys.stderr, '[WARN] CDS of {} is not Multiple of 3 with frame {}; discarded..'.format(RNARecord.id, frame)
				continue
			else:
				try: assert frame == 0, '[WARN] {}: frame ({}) of the first CDS region is not 0'.format(RNARecord.id, frame)
				except AssertionError: pass
			for i in range(frame, len(cds_pos), 3):
				codon = cds_seq[i:i+3]
				cod_pos = cds_pos[i:i+3]
				for j, pos in zip(range(3), cod_pos):
					postype = fold_codon(codon, j)
					if postype is None:
						continue
					postype = 'codon{}-{}'.format(j+1, postype)
					try: d_positon[pos].add(postype)
					except KeyError: d_positon[pos] = {postype}
		return d_positon
	def classify_fold2(self, d_seqs, d_positon=OrderedDict()):
		'''exclude splicing sites'''
		if not self.is_coding:
			return d_positon
		seq = d_seqs[self.chrom]
		for RNARecord in GffRNARecords(self):
			for cds in RNARecord.features['CDS']:
				frame = cds.frame
				cds_seq = cds.extract_seq(seq)	# coding strand
				cds_pos = cds.positions			# coding strand
				cds_seq = cds_seq[frame:]
				cds_pos = cds_pos[frame:]
				assert len(cds_seq) == len(cds_pos)
				for i in range(0, len(cds_pos), 3):
					codon = cds_seq[i:i+3]
					cod_pos = cds_pos[i:i+3]
					for j, pos in zip(range(3), cod_pos):
						postype = fold_codon(codon, j)
						if postype is None:
							continue
						try: d_positon[pos].add(postype)
						except KeyError: d_positon[pos] = {postype}
		return d_positon

def fold_codon(codon, index):
	if len(codon) != 3:
		return None
	my_codon = list(codon)
	aas = []
	for base in 'ATCG':
		my_codon[index] = base
		aa = Seq(''.join(my_codon)).translate()
		aa = str(aa)
		aas += [aa]
	aas = set(aas)
	if len(aas) == 4:
		return 'fold0'
	if len(aas) == 1:
	#   print >> sys.stderr, codon, index
		return 'fold4'
	return None
		
class GffRNARecords():
	def __init__(self, GeneRecord):
		self._GeneRecord = GeneRecord
	def __iter__(self):
		return self._parse()
	def _parse(self):
		self.GeneRecord = copy.deepcopy(self._GeneRecord)
		if self.GeneRecord.gene.type == 'gene':
			if len(self.GeneRecord.lines) == 1:
				pass
			else:
				self.GeneRecord.remove_node(self.GeneRecord.id)
		for cmpt in nx.connected_components(self.GeneRecord.to_undirected()):
			yield self.GeneRecord.subgraph(cmpt)
	def to_gene(self, rnas):
		nodes = [self._GeneRecord.id]
		for rna in rnas:
			nodes += rna.nodes()
		remove_nodes = set(self._GeneRecord.nodes()) - set(nodes)
		if remove_nodes:
			self._GeneRecord.remove_nodes_from(remove_nodes)
		return self._GeneRecord

def main():
	pass
if __name__ == '__main__':
	main()
