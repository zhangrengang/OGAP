import sys, os
import sqlite3
import json
from collections import OrderedDict
from lazy_property import LazyWritableProperty as lazyproperty
from RunCmdsMP import run_cmd
from small_tools import open_file as open

rootdir = os.path.dirname(os.path.realpath(__file__))
dbdir = '{}/../db'.format(rootdir)
DB = '{}/taxonomy.json.gz'.format(dbdir)

class Taxonomy():
	def __init__(self, spname=None, taxid=None,
				jsonfile=DB, load=True, 
				dbfile=None):
		self.jsonfile = jsonfile
		if os.path.exists(jsonfile) and load:
			self.load_db()
		elif dbfile is None or not os.path.exists(jsonfile):
			homedir = os.environ['HOME']
			#homedir = '/export/tmp'
			dbfile = '{}/.etetoolkit/taxa.sqlite'.format(homedir)
			if not os.path.exists(dbfile):
				cmd = 'ete3 ncbiquery --search 9606 --info'
				run_cmd(cmd, log=True)
			if load:
				self.dbfile = dbfile
				self.dump_db()
				self.load_db()
			assert os.path.exists(dbfile), dbfile
		self.taxid = taxid
		self.spname = spname
		self.dbfile = dbfile

	def load_db(self):
		print >>sys.stderr, 'loading from {}'.format(self.jsonfile)
		with open(self.jsonfile, 'r') as fin:
			self.db = json.load(fin)
	def dump_db(self):
		d_results = {}
		d_names = {}
		for table, obj in zip(['species', 'synonym'], [Species, Synonym]):
			query = 'SELECT * FROM {};'.format(table)
			print >>sys.stderr, query
			conn = sqlite3.connect(self.dbfile)
			curs = conn.cursor()
			curs.execute(query)
			for result in curs.fetchall():
				result = obj(result)
				if result.taxid not in d_results:
					d_results[result.taxid] = result
				d_names[result.spname] = result.taxid
		d_species = {}
		print >>sys.stderr, 'collecting taxonomy'
		for spname, taxid in d_names.items():
			track = map(int, d_results[taxid].track.split(','))
			ranks = [d_results[tid].rank for tid in track]
			if not 'species' in set(ranks):
				continue
			taxonomy = [d_results[tid].spname for tid in track]
			if "Eukaryota" not in set(taxonomy):
				continue
			ranks = ['' if rank == 'no rank' else rank for rank in ranks]
			taxonomy = taxonomy[::-1]	# reverse
			ranks = ranks[::-1]
			d_species[spname.lower()] = (taxid, taxonomy, ranks)
				
		print >>sys.stderr, 'dumping into {}'.format(self.jsonfile)
		with open(self.jsonfile, 'w') as fout:
			json.dump(d_species, fout)

	def query_db(self, taxid=None, spname=None):
		'''TABLE species (taxid INT PRIMARY KEY, parent INT, spname VARCHAR(50) COLLATE NOCASE, common VARCHAR(50) COLLATE NOCASE, rank VARCHAR(50), track TEXT)'''
		column = '*'
		query = 'SELECT {column} FROM species where {key}={value};'
		if taxid is not None:
			key, value = 'taxid', taxid
		elif spname is not None:
			key, value = 'spname', spname
		query = query.format(column=column, key=key, value=repr(value))
		conn = sqlite3.connect(self.dbfile)
		curs = conn.cursor()
		#print >>sys.stderr, query
		curs.execute(query)
		result = curs.fetchone()
		if result is not None:
			return Species(result)
		else:	# search in synonym
			if taxid is None and spname is not None:
				# synonym: taxid INT,spname VARCHAR(50)
				query = 'SELECT {column} FROM synonym where spname={spname};'
				query = query.format(column=column, spname=repr(spname))
				conn = sqlite3.connect(self.dbfile)
				curs = conn.cursor()
				curs.execute(query)
				result = curs.fetchone()
				if result is not None:
					taxid = Synonym(result).taxid
					return self.query_db(taxid=taxid)
	@lazyproperty
	def is_family(self):
		return self.is_rank('family')
	@lazyproperty
	def is_order(self):
		return self.is_rank('order')
	def is_rank(self, rank, spname=None, taxid=None):
		if taxid is None:
			taxid = self.taxid
		if spname is None:
			spname = self.spname
		rrank = self.query_rank(taxid=taxid, spname=spname)
		return rrank == rank
	def get_rank(self, rank):
		track = self.get_track()
		if not track:
			return 
		for taxid in track[1:]:
			if self.is_rank(rank, taxid=taxid):
				return taxid
	def get_ranks(self):
		track = self.get_track()
		d = OrderedDict()
		for taxid in track:
			rank = self.query_rank(taxid=taxid)
			if rank == 'no rank':
				continue
			d[rank] = taxid
		return d
	def query_rank(self, taxid=None, spname=None):
		result = self.query_db(taxid=taxid, spname=spname)
		rank = result.rank
		return rank
	def get_track(self):
		result = self.query_db(taxid=self.taxid, spname=self.spname)
		try:
			track = map(int, result.track.split(','))
		except AttributeError:
			print >> sys.stderr, '{} queried None'.format((self.taxid, self.spname))
			track = []
		return track
	def get_taxonomy(self):
		track = self.get_track()
		if not track:
			return
		taxonomy = [self.query_db(taxid=taxid).spname for taxid in track]
		return taxonomy

	@lazyproperty
	def track(self):
		track = self.get_track()
		if not track:
			return []
		return [self.query_db(taxid=taxid) for taxid in track]
	@lazyproperty
	def taxonomy(self):	# revesed
		return [sp.spname for sp in self.track][::-1]
	@lazyproperty
	def rank(self):
		return [sp.rank for sp in self.track][::-1]
	@lazyproperty
	def ranks(self):
		return OrderedDict([(rank, taxon) for rank, taxon in zip(self.rank, self.taxonomy)])

class Species():
	def __init__(self, record=None):
		'''a record is fetchone from TABLE species'''
		self.record = record
		if record is None:
			self = None
		else:
			(self.taxid, self.parent, self.spname, 
				self.common, self.rank, self.track) = record
class Synonym():
	def __init__(self, record=None):
		if record is None:
			self = None
		else:
			self.taxid, self.spname = record	
def get_info(lstfile, fout=sys.stdout):
#	d = Taxonomy().db
	d_genus = {}
	for line in open(lstfile):
		sp0 = line.split()[0]
		sp = sp0.replace('_', ' ') #.lower()
		g = sp.split()[0]
		if g in d_genus:
			rank = d_genus[g]
		else:
			db = Taxonomy(spname=g, load=False)
			d_tax = db.ranks
			#taxid, taxonomy, ranks = d[sp]
			#d_tax = dict(zip(ranks, taxonomy))
			try: d_genus[g] = rank = [d_tax['order'], d_tax['family'], ','.join(db.taxonomy)]
			except KeyError:rank = ['error']
		line = [sp0] + rank
		print >>fout, '\t'.join(line)

def main():
	get_info(sys.argv[1])
	return
	#Taxonomy(jsonfile='db/taxonomy.json')
	d = Taxonomy(jsonfile='db/taxonomy.json.gz').db
	print len(d)
	for k in d:
		print d[k]
		break
if __name__ == '__main__':
	main()
