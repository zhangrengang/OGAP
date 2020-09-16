import sys
import re
import csv
from collections import OrderedDict, Counter
ITEMTYPE = OrderedDict([
	('Run', str),
	('ReleaseDate', str),
	('LoadDate', str),
	('spots', int),
	('bases', int),
	('spots_with_mates', int),
	('avgLength', int),
	('size_MB', int),
	('AssemblyName', str),
	('download_path', str),
	('Experiment', str),
	('LibraryName', str),
	('LibraryStrategy', str),
	('LibrarySelection', str),
	('LibrarySource', str),
	('LibraryLayout', str),
	('InsertSize', int),
	('InsertDev', int),
	('Platform', str),
	('Model', str),
	('SRAStudy', str),
	('BioProject', str),
	('Study_Pubmed_id', int),
	('ProjectID', int),
	('Sample', str),
	('BioSample', str),
	('SampleType', str),
	('TaxID', int),
	('ScientificName', str),
	('SampleName', str),
	('g1k_pop_code', str),
	('source', str),
	('g1k_analysis_group', str),
	('Subject_ID', str),
	('Sex', str),
	('Disease', str),
	('Tumor', str),
	('Affection_Status', str),
	('Analyte_Type', str),
	('Histological_Type', str),
	('Body_Site', str),
	('CenterName', str),
	('Submission', str),
	('dbgap_study_accession', str),
	('Consent', str),
	('RunHash', str),
	('ReadHash', str),
])

class SraRunInfo(object):
	def __init__(self, run_info=None):
		self.run_info = run_info
	def __iter__(self):
		return self.parse()
	def parse(self):
		with open(self.run_info) as fp:
			for line in csv.reader(fp):
		#		print line
				if not line:
					continue
				if line[0] == ITEMTYPE.keys()[0]:
					continue
				yield RunInfoRecord(line)
	def filter(self, critera):
		string = re.compile(r'(\w+)\s*([>=<]+)').sub(r'record.\1 \2', critera)
#		print string
		for record in self.parse():
			result = eval(string)
#			print result
#			print record.dict
#			print record.LibrarySelection, record.avgLength, record.bases, record.Sample
			if result is True:
				yield record
	def write(self, records=None, fout=sys.stdout):
		writer = csv.writer(fout)
		writer.writerow(ITEMTYPE.keys())
		for record in records:
			writer.writerow(record.line)
	def select(self, columns, records=None, fout=sys.stdout):
		if records is None:
			records = self.parse()
#		print >> fout, '\t'.join(columns)
		for record in records:
			line = [record.dict[col] for col in columns]
#			line = map(str, line)
#			print >> fout, '\t'.join(line)
			yield RunInfoRecord2(columns, line)
	def add_taxonomy(self, fout=sys.stdout):
		from Taxonomy import Taxonomy
		from Database import Database
		taxonomy_dbfile = Database().taxonomy_dbfile
		d_species = Taxonomy(jsonfile=taxonomy_dbfile).db
		records = []
		for record in self:
			sp = record.ScientificName.lower()
			try: taxid, record.taxonomy, record.ranks = d_species[sp]
			except KeyError: record.taxonomy = []
			record.line += record.taxonomy
			records += [record]
		self.write(records)
class RunInfoRecord(object):
	def __init__(self, line):
		self.line = line
		self.dict = OrderedDict()
		for (key, func), value in zip(ITEMTYPE.items(), line):
			try:
				value = func(value)
			except ValueError:
				value = None
			setattr(self, key, value)
			self.dict[key] = value
	def write(self, fout):
		print >>fout, '\t'.join(self.line)
	def write_head(self, fout):
		print >>fout, '\t'.join(ITEMTYPE.keys())
class RunInfoRecord2(object):
	def __init__(self, keys, values):
		self.values = values
		self.dict = OrderedDict()
		for key, value in zip(keys, values):
			self.dict[key] = value
			setattr(self, key, value)
	def __str__(self):
		self.values = map(str, self.values)
		return '\t'.join(self.values)
def add_taxonomy(run_info, fout=sys.stdout):
	SraRunInfo(run_info).add_taxonomy(fout=fout)

def stats_by_bioproject(run_info, fout=sys.stdout):
	import numpy as np
	d_porject = {}
	for run in SraRunInfo(run_info):
		proj = run.dict['BioProject']
		try: d_porject[proj].append(run)
		except KeyError: d_porject[proj] = [run]

	columns = ['BioProject', 'Runs', 'BioSamples', 'Taxas', 'Taxa', 'LibrarySelection', 'Platforms', 'Platform',
			   'mean_bases', 'median_bases', 'mean_avgLength', 'median_avgLength']
	print >>fout, '\t'.join(columns)
	for proj, runs in sorted(d_porject.items(), key=lambda x:len(x[1]), reverse=1): # sort by RUN number
		sra_num = len(runs)
		sample_num = len({run.dict['BioSample'] for run in runs})
		taxa = [run.dict['ScientificName'] for run in runs]
		taxa_count = Counter(taxa)
		taxa_num = len(taxa_count)
		taxa_count = sorted(taxa_count.items(), key=lambda x:x[1], reverse=1)
		selection = [run.dict['LibrarySelection'] for run in runs]
		selc_count = Counter(selection)
		selc_count = sorted(selc_count.items(), key=lambda x:x[1], reverse=1)
		bases = [run.dict['bases'] for run in runs]
		bases_mean = np.mean(bases)
		bases_median = np.median(bases)
		avgLengths = [run.dict['avgLength'] for run in runs]
		len_mean = np.mean(avgLengths)
		len_median = np.median(avgLengths)
		platforms = [run.dict['Platform'] for run in runs]
		pf_count = Counter(platforms)
		pf_num = len(pf_count)
		pf_count = sorted(pf_count.items(), key=lambda x:x[1], reverse=1)
		line = [proj, sra_num, sample_num, taxa_num, taxa_count, selc_count, pf_num, pf_count, bases_mean, bases_median, len_mean, len_median]
		line = map(str, line)
		print >> fout, '\t'.join(line)
def filter_and_select_runinfo(
		run_info, 
		fout=sys.stdout,
		critera='LibrarySelection == "RANDOM" and avgLength >= 72 and bases >= 100e6 and SampleType == "simple" and TaxID > 1',
		columns=['Run', 'BioSample', 'TaxID', 'ScientificName', 				# 4
				 'Platform', 'Model', 'spots', 'bases', 'avgLength', 			# 9
				 'BioProject', 'ReleaseDate', 'LoadDate',						# 12
				 'LibraryLayout', 'SRAStudy', 'SampleName', 'LibraryName',],	# 16
		):
	filter_out = SraRunInfo(run_info).filter(critera)
	#print filter_out
	print >> fout, '\t'.join(columns)
	lines = []
	for line in SraRunInfo().select(columns, records=filter_out, fout=fout):
		print >>fout, line
		lines.append(line)
	return lines
def get_singletons(run_info, fout=sys.stdout):
	lines = filter_and_select_runinfo(run_info, fout=sys.stderr)
	taxa = [line.TaxID for line in lines]
	taxa_count = Counter(taxa)
	singletons = {tax for tax,count in taxa_count.items() if count == 1} # set
	for line in lines:
		if line.TaxID in singletons:
			print >>fout, line
	return singletons
def get_last(run_info, fout=sys.stdout):
	lines = filter_and_select_runinfo(run_info, fout=sys.stderr)
	d_last = {}
	for line in lines:
		key = (line.ScientificName, line.Platform)
		if key in d_last and line.BioSample != d_last[key].BioSample:
			continue
		d_last[key] = line
		print >>fout, line
if __name__ == '__main__':
	subcmd = sys.argv[1]
	run_info = sys.argv[2]
	if subcmd == 'stats_bioproject':
		stats_by_bioproject(run_info)
	elif subcmd == 'get_singletons':
		#filter_and_select_runinfo(run_info)
		get_singletons(run_info)
	elif subcmd == 'get_last':
		#filter_and_select_runinfo(run_info)
		get_last(run_info)
	elif subcmd == 'add_taxonomy':
		add_taxonomy(run_info)
	else:
		raise ValueError
