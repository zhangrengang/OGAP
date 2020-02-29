import sys, os
import time
import shutil
import glob
import argparse
from collections import Counter
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Genbank import GenbankParser
from OrthoFinder import OrthoFinder, OrthoMCLGroupRecord
from RunCmdsMP import logger, run_cmd, run_job
from small_tools import mkdirs, rmdirs

AUGUSTUS_CONFIG_PATH=os.environ['AUGUSTUS_CONFIG_PATH']

class DB():
	def __init__(self, organ, taxon=None, dbdir=None,
				gbfiles=None, custom=None, 
				tmpdir='/dev/shm/tmp', 
				ncpu=16,
				**kargs):
		if custom is None:
			if taxon is None:
				ogtype = '{}-custom'.format(dbdir, organ)	# db name and prefix
			else:
				ogtype = '{}-{}'.format(organ, taxon)  # e.g. 'mt-Lamiales'
		else:
			ogtype = '{}-{}'.format(organ, custom)
		self.ogtype = ogtype
		self.organ = organ
		self.dbdir = dbdir
		self.taxon = taxon
		self.ncpu = ncpu
		self.gbfiles = gbfiles

		# folders
		outdir = '{}/{}'.format(dbdir, ogtype)
		tmpdir = '{}/{}'.format(tmpdir, ogtype)
		self.outdir = outdir = os.path.realpath(outdir)
		self.tmpdir = tmpdir = os.path.realpath(tmpdir)
		self.pepdir = '{}/prot'.format(tmpdir)
		self.rnadir = '{}/rrna'.format(tmpdir)
		self.alndir = '{}/alns'.format(tmpdir)
		self.trndir = '{}/trng'.format(tmpdir)
		self.msadir = '{}/msas'.format(outdir)
		self.hmmdir = '{}/hmms'.format(outdir)
		self.pfldir = '{}/prfl'.format(outdir)
		self.tnddir = '{}/agst'.format(outdir)

		# files
		self.vesion_file = '{}/version'.format(self.outdir)
		self.cds_file = '{}/cds.info'.format(self.outdir)
		self.rna_file = '{}/rna.info'.format(self.outdir)
		self.species_file = '{}/species.info'.format(self.outdir)
		self.taxonomy_file = '{}/taxonomy.info'.format(self.outdir)

	def makedb(self):
		rmdirs(self.pepdir, self.rnadir)  # for OrthoFinder
		mkdirs(self.outdir, self.pepdir, self.rnadir, self.alndir, \
				self.msadir, self.hmmdir, self.pfldir, self.trndir, self.tnddir)
		logger.info('parsing {}'.format(self.gbfiles))
		gb = GenbankParser(self.gbfiles)
		version = time.strftime("%Y-%m-%d", time.gmtime(os.path.getmtime(self.gbfiles[0])))
		cmd = 'echo RefSeq {} > {}/version'.format(version, self.outdir)
		run_cmd(cmd)
		ofbin = '/share/home/app/build/OrthoFinder-2.3.1/OrthoFinder-2.3.1-beta/orthofinder/orthofinder.py'
		for feat_type, seqdir, prefix, of_opts, build_pfl, min_ratio, cdhit, out_info in zip(
				['protein', ('rRNA','tRNA')],
				[self.pepdir, self.rnadir],
				[self.cds_file, self.rna_file],
				['', '-S blastn'],
				[True, False],
				[0.1, 0.05],
				['cd-hit', 'cd-hit-est'],
				[True, False],
				):
			self.cdhit = cdhit
			logger.info('writing {} sequences into {}'.format(feat_type, seqdir))
			mkdirs(seqdir)
			records = list(gb.filter_by_taxon(taxon=self.taxon))	# read into MEM
			features = gb.write_seqs_by_taxon(seqdir, feat_type=feat_type, records=records, taxon=self.taxon)
			gb.ntaxa = len(set(gb.taxa))	# some taxa have >1 records
			if out_info:
				logger.info('{} taxa: {}'.format(gb.ntaxa, sorted(set(gb.taxa))))
				with open(self.species_file, 'w') as fout:
					print >> fout, '\t'.join(['#accession', 'organism', 'taxonomy'])
					for taxon, acc, taxonomy in sorted(zip(gb.taxa, gb.accs, gb.taxonomy)):
						taxonomy = ';'.join(taxonomy)
						print >> fout, '\t'.join([acc, taxon, taxonomy])
				info_file = '{}/taxonomy.info'.format(self.outdir)
				taxa = [';'.join(tax) for tax in gb.tax]
				taxon, taxon_count = self.get_most_common(taxa)
				logger.info('get taxon `{}` from {}'.format(taxon, taxon_count))
				tables = [feat.transl_table for feat in features]
				table, table_count = self.get_most_common(tables)
				logger.info('get transl_table `{}` from {}'.format(table, table_count))
				title = ["#organelle", "taxon", "transl_table", "taxonomy"]
				line = [self.organ, self.taxon, table, taxon]
				line = map(str, line)
				with open(self.taxonomy_file, 'w') as fout:
					print >> fout, '\t'.join(title)
					print >> fout, '\t'.join(line)

			logger.info('running OrthoFinder to cluster genes')
			cmd = '{} -f {} -M msa -T fasttree -t {} {}'.format(ofbin, seqdir, self.ncpu, of_opts)
			run_cmd(cmd, log=True)

			logger.info('parsing OrthoFinder orthologs')
			orfdir = glob.glob('{}/OrthoFinder/*/'.format(seqdir))[0]
			orf = OrthoFinder(orfdir)
			logger.info('parsing orthologs group into {}'.format(self.alndir))
			min_count = gb.ntaxa * min_ratio
			self.min_count = 3 if min_count < 3 else min_count
			info_file = prefix + '.raw'

			# parse orthologs, bin sequences by orthologs
			with open(info_file, 'w') as fout:
				gene_names, seqfiles, groups = \
						self.parse_orthologs(orf, features, self.alndir, fout) # names
			if not os.path.exists(prefix):
				os.rename(info_file, prefix)
			logger.info('extrat {} genes: {}'.format(len(gene_names), gene_names))		# gene id

			logger.info('aligning sequences')
			alnfiles = self.align_seqs(gene_names, seqfiles, outdir=self.msadir)
			logger.info('building hmm profiles')
			hmmfiles = self.hmmbuild_alns(gene_names, alnfiles, outdir=self.hmmdir)
			if build_pfl:
				logger.info('building augustus profiles')
				pflfiles = self.msa2prfl(gene_names, alnfiles, outdir=self.pfldir)
				logger.info('training augustus models')
				# train augustus
				d_features = {feat.id: feat for feat in features}
				for gene_name, group in zip(gene_names, groups):
					training_features = {d_features[gene] for gene in group}
					train_db = '{}/train-{}.gb'.format(self.trndir, gene_name)
					with open(train_db, 'w') as fout:
						gene_number = self.get_augustus_training_set(records, training_features, fout)
					sepcies = '{}-{}'.format(self.ogtype, gene_name)
					self.train_augustus(train_set=train_db, species=sepcies, outdir=self.tnddir,
						gene_number=gene_number, ncpu=self.ncpu)

	def parse_orthologs(self, orf, features, outdir, fout):
		d_seqs = {feat.id: feat.seq for feat in features}
		d_names = {feat.id: feat.gene for feat in features}
		#print >>sys.stderr, d_names
		d_products = {feat.id: feat.product for feat in features}
		d_types = {feat.id: feat.type for feat in features}
		d_nexon = {feat.id: feat.nexon for feat in features}
		d_trans_splicing = {feat.id: feat.trans_splicing for feat in features}

		logger.info('using min_count={} to get gene name'.format(self.min_count))
		line = ['#id', 'name', 'product', 'seq_type', 
				'exon_number', 'trans_splicing', 
				'gene_count', 'species_count', 'names', 'species']
		print >>fout, '\t'.join(line)
		d_genes = {}
		lines = []
		gene_names, outseqs, groups = [], [], []
		for i, group in enumerate(orf.get_orthologs_cluster()):		# group by feat.id
			gene_name, name_count = self.get_gene_name(group, d_names,)
			tmpfix = '{}/{}.para'.format(outdir, gene_name)
			group = self.filter_paralogs(group, d_seqs, tmpfix)
			# name
			gene_name, name_count = self.get_gene_name(group, d_names,)
			logger.info('get gene name `{}` from {}'.format(gene_name, name_count))
			max_count = name_count[gene_name]
			if max_count < self.min_count:
				continue
			gene_id = gene_name
			if gene_id in d_genes:   # for duplicated gene names
				gene_id = '{}-{}'.format(gene_name, i)
			d_genes[gene_id] = (gene_name, group)

			# product
			product, product_count = self.get_gene_name(group, d_products)
			logger.info('get product `{}` from {}'.format(product, product_count))
			top_names = name_count.most_common(25)
			top_names = ','.join(['{}:{}'.format(name, count) for name, count in top_names])
			# count
			gene_count = len(group)
			species = [gene.split('|')[0] for gene in group]
			species_count = Counter(species)
			species = ','.join(['{}:{}'.format(name, count) \
								for name, count in sorted(species_count.items())])
			species_count = len(species_count)
			# exon number and trans_splicing
			nexon, nexon_count = self.get_gene_name(group, d_nexon)
			trans_splicing, ts_count = self.get_gene_name(group, d_trans_splicing)
			logger.info('get exon number `{}` from {}'.format(nexon, nexon_count))
			logger.info('get trans_splicing `{}` from {}'.format(trans_splicing, ts_count))

			# type
			seqtype, _ = self.get_gene_name(group, d_types)
			# write info
			line = [gene_id, gene_name, product, seqtype, 
					nexon, trans_splicing,
					gene_count, species_count, top_names, species]
			lines += [line]
			# write sequences
			outseq = '{}/{}.fa'.format(outdir, gene_name)
			self.write_orthologs(gene_name, group, d_seqs, outseq)	# by group
			outseqs += [outseq]
			gene_names += [gene_id]
			groups += [group]

		for line in sorted(lines):
			line = map(str, line)
			print >>fout, '\t'.join(line)
		return gene_names, outseqs, groups

	def train_augustus(self, train_set, species, gene_number, outdir, ncpu, transl_table=1):
		gene_number = gene_number/3
		kfold = 8 if ncpu > 8 else ncpu
		spdir = '{}/species/'.format(AUGUSTUS_CONFIG_PATH,)
		cmd = '''# train
(cd {spdir} && rm {species} -r)
randomSplit.pl {train_set} {gene_number}
new_species.pl --species={species}
etraining --species={species} {train_set}.train --translation_table={transl_table}
augustus --species={species} {train_set}.test --translation_table={transl_table} &> {train_set}.train_test1.out
# optimize
#optimize_augustus.pl --species={species} --cpus={ncpu} --kfold={kfold} --rounds=5 {train_set}.train && rm -r tmp_opt_{species}
#etraining --species={species} {train_set}.train
#augustus --species={species} {train_set}.test &> {train_set}.train_test2.out
(cd {spdir} && tar cvzf {outdir}/{species}.tgz {species})
'''.format(train_set=train_set, gene_number=gene_number, 
			species=species, ncpu=ncpu, kfold=kfold,
			outdir=outdir, spdir=spdir,
			transl_table=transl_table)
		run_cmd(cmd, log=True)
		
	def get_augustus_training_set(self, records, training_features, fout, flank=1000):
		'''get small DNA'''
		gene_number = 0
		for rc in records:
			record = rc.record
			indices = []
			trans_splicing_records = []
			for feat in rc.get_features('protein'):
				if not feat.pep[0] == 'M':  # start codon
					continue
				if feat in training_features:
					if feat.trans_splicing:     # checked all trans_splicing
						trans_splicing_records += [feat.rearrange_trans_splicing(record)]
						continue
					indices += [feat.index]
			features = [record.features[index] for index in indices]
			for i, feature in enumerate(features):
				start = feature.location.start - flank
				end = feature.location.end + flank
				start = 0 if start < 0 else start
				feat_record = record[start:end]		# sub-record
				feat_record.features = [feat for feat in feat_record.features \
								if (feat.location.start == feature.location.start - start and \
								    feat.location.end == feature.location.end - start)]
				if feature.strand == -1:
					feat_record = feat_record.reverse_complement()
				feat_record.name = feat_record.id = '{}_{}-{}'.format(record.id, start, end)
				if not feat_record.features[0].type == 'source':
					source = SeqFeature(
								FeatureLocation(0, len(feat_record.seq)), 
								type='source', strand=1)
					feat_record.features = [source] + feat_record.features
				SeqIO.write(feat_record, fout, 'genbank')
			for feat_record in trans_splicing_records:
				SeqIO.write(feat_record, fout, 'genbank')
			gene_number += len(indices) + len(trans_splicing_records)
		return gene_number

	def align_seqs(self, gene_names, seqfiles, outdir):
		alnfiles = []
		for gene_name, seqfile in zip(gene_names, seqfiles):
			outaln = '{}/{}.aln'.format(outdir, gene_name)
#			#cmd = 'mafft --auto {} > {}'.format(seqfile, outaln)
			cmd = 'mafft --auto {} | prepareAlign | mafft --auto - > {}'.format(seqfile, outaln)
			run_cmd(cmd)
			alnfiles += [outaln]
		return alnfiles
	def hmmbuild_alns(self, gene_names, alnfiles, outdir):
		hmmfiles = []
		for gene_name, alnfile in zip(gene_names, alnfiles):
			outhmm = '{}/{}.hmm'.format(outdir, gene_name)
			cmd = 'hmmbuild -n {} {} {}'.format(gene_name, outhmm, alnfile)
			run_cmd(cmd)
			hmmfiles += [outhmm]
		return hmmfiles
	def msa2prfl(self, gene_names, alnfiles, outdir):
		prflfiles = []
		for gene_name, alnfile in zip(gene_names, alnfiles):
			outprfl = '{}/{}.prfl'.format(outdir, gene_name)
			cmd = 'msa2prfl.pl {} --setname={} > {}'.format(alnfile, gene_name, outprfl)
			run_cmd(cmd)
			prflfiles += [outprfl]
		return prflfiles
	def write_orthologs(self, gene_name, group, d_seqs, outseq):
		with open(outseq, 'w') as fout:
			for gene in sorted(group):
				print >>fout, '>{}\n{}'.format(gene, d_seqs[gene])
	def filter_paralogs(self, genes, d_seqs, tmpfix):
		'''for paralogs, cluster for respective ones'''
		group = OrthoMCLGroupRecord(genes=genes)
		single_genes = []
		for sp, genes in group.spdict.items():
#			#gene = max(genes, key=lambda x: len(d_seqs[x]))
			if len(genes) > 1:
				tmpseq = '{}.paralogs'.format(tmpfix)
				clstseq = '{}.paralogs.clust'.format(tmpfix)
				with open(tmpseq, 'w') as f:
					for gene in sorted(genes):
						print >>f, '>{}\n{}'.format(gene, d_seqs[gene])
				genes = self.cluster_genes(self.cdhit, tmpseq, clstseq)
			single_genes += genes
		return single_genes
	def cluster_genes(self, cdhit, input, output):
		cmd = '{} -i {} -o {} -c 0.95 -aS 0.95 -G 0 -T 4 -d 0 -M 2000'.format(cdhit, input, output)
		run_cmd(cmd)
		return [rc.id for rc in SeqIO.parse(output, 'fasta')]
	def get_gene_name(self, genes, d_names, min_count=3):
		'''get gene name by max count'''
		gene_names = [d_names[gene] for gene in genes]
		return self.get_most_common(gene_names)
	def get_most_common(self, names):
		name_count = Counter(names)
		max_name, max_count = max(name_count.items(), key=lambda x: x[1])
		return max_name, name_count

def get_defualt_dbdir():
	rootdir = os.path.dirname(os.path.realpath(__file__))
	dbdir = '{}/../db'.format(rootdir)
	return os.path.realpath(dbdir)

def main():
	'''example:
-gbfiles *.gb[.gz] -organ mt -taxon rosids				  # use records of rosids and name db as mt-rosids
-gbfiles *.gb[.gz] -organ mt -taxon rosids -custom myself # use records of rosids but name db as mt-myself
-gbfiles *.gb[.gz] -organ mt -custom myself               # use all records in gbfiles and name db as mt-myself
'''
	args = makeArgparse()
	print >>sys.stderr, args.__dict__
	DB(**args.__dict__).makedb()

def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-gbfiles', type=str, nargs='+', required=True, default=None, 
						help="genbank files")
	parser.add_argument('-organ', type=str, choices=['mt', 'pt'], default='mt',
						help="mt (mitochondrion) or pt (plastid) to build [default=%(default)s]")
	parser.add_argument('-dbdir', type=str, default=get_defualt_dbdir(),
						help="directory to output [default=%(default)s]")
	parser.add_argument('-taxon', type=str, default=None,
						help="taxon to filter out, such as Embryophyta [default=%(default)s]")
	parser.add_argument('-custom', type=str, default=None,
						help="use custom mark instead taxon [default=%(default)s]")
	args = parser.parse_args()
	return args

if __name__ == '__main__':
	main()
