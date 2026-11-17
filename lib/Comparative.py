from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import object
from past.utils import old_div
import sys, os
import re
import glob
import itertools
from Bio import SeqIO
from .RunCmdsMP import run_cmd, logger
from .OrthoFinder import catAln
from .small_tools import mkdirs, parse_kargs, tr_numeric
from itertools import combinations

class PhyloPipeline(object):
	def __init__(self, indir, 
			outprefix=None,
			gnid = True,
			tmpdir='./tmp',
			types=['cds'], #, 'rna'], 
			root=None,
			max_gene_missing=50,
			min_shared=50,	# min_taxa_missing
			min_cov=50, # min coverage
			):
		self.indir = indir
		self.tmpdir = tmpdir
		self.types = types
		self.min_shared = min_shared
		self.max_gene_missing = max_gene_missing
		self.min_cov = min_cov
		self.outprefix = outprefix
		self.gnid = gnid
		self.root = root
		if outprefix is None:
			self.outprefix = '{}-{}'.format(
				os.path.basename(indir.strip('/')), '-'.join(types))
		
	def run(self):
		mkdirs(self.tmpdir)
		self.d_gnid = self.get_genome_id()

		# bin
		self.genes = self.bin_seqs()
		# align
		alnfiles = self.align_seqs()
		# concat
		concat_alignment = '{}.aln'.format(self.outprefix)
		with open(concat_alignment, 'w') as fout:
			catAln(alnfiles, fout)
		# trim
		trimed_alignment = concat_alignment + '.trimal'
		self.trimal(concat_alignment, trimed_alignment)
		# tree
		self.iqtree(trimed_alignment)

	def iqtree(self, alnfile, opts=''):
		if self.root is not None:
			opts += ' -o {}'.format(self.root)
#		cmd = 'iqtree -redo -s {} -bb 1000 {} -nt AUTO > /dev/null'.format(alnfile, opts)
#		cmd = 'iqtree2 -redo -s {} -B 1000 {} -T 20 > /dev/null'.format(alnfile, opts)
		cmd = 'iqtree -redo -s {} -B 1000 {} -T 20 > /dev/null'.format(alnfile, opts)
		run_cmd(cmd, log=True)

	def trimal(self, inaln, outaln):
		cmd = 'trimal -in {} -automated1 > {}'.format(inaln, outaln)
		run_cmd(cmd)

	def align_seqs(self):
		alnfiles = []
		for gene in self.genes:
			seqfile = self.get_seqfile(gene)
			alnfile = self.get_alnfile(gene)
			cmd = 'mafft --auto {seq} > {aln}'.format(seq=seqfile, aln=alnfile)
			run_cmd(cmd)
			alnfiles += [alnfile]
		return alnfiles

	def bin_seqs(self):
		fastas = []
		prefixs = []
		for type in self.types:
			fastafiles, prefies = self.get_fastas(type)
			prefixs += prefies
			fastas += fastafiles
		d_seqs = {}
	
		d_sps = {}
		genes = set([])
		for rc, prefix in self.get_seqs(fastas, prefixs):
			gene = rc.id
			genes.add(gene)
			rc.taxon = rc.id = prefix
			if self.gnid:
				suffix = '-{}'.format(self.d_gnid[prefix])
				rc.id += suffix
				if self.root == prefix:
					self.root = rc.id
			try: d_seqs[gene] += [rc]
			except KeyError: d_seqs[gene] = [rc]
			try: d_sps[prefix] += [gene]
			except KeyError: d_sps[prefix] = [gene]
		# remove taxa with too many missing genes
		sp_to_remove = set([])
		for sp, _genes in list(d_sps.items()):
			non_missing = old_div(1e2*len(_genes), len(genes))
			logger.info('{}\t{}'.format(sp, non_missing))
			if (100-non_missing) > self.max_gene_missing:
				sp_to_remove.add(sp)
				logger.info('taxon {} only contains {}% genes, removed.\
                            '.format(sp, non_missing))
				continue
		logger.info('to be remove: {}'.format(sp_to_remove))
		nspecies = len(set(prefixs)) - len(sp_to_remove) 
		# write
		gene_to_remove = set([])
		for gene, records in list(d_seqs.items()):
			#print >> sys.stderr, [record.taxon for record in records][:10]
			records = [record for record in records if not record.taxon in sp_to_remove]
			shared = old_div(1e2*len(records), nspecies)
			logger.info('{}\t{}'.format(gene, shared))
			if shared < self.min_shared:
				gene_to_remove.add(gene)
				logger.info('gene {} is only shared by {}% species, removed.\
							'.format(gene, shared))
				continue
			seqfile = self.get_seqfile(gene)
			with open(seqfile, 'w') as fout:
				for record in records:
					SeqIO.write(record, fout, 'fasta')
		logger.info('removed taxon: {}/{}; removed genes: {}/{}'.format(
			len(sp_to_remove), len(d_sps), len(gene_to_remove), len(d_seqs)))
		return list(d_seqs.keys())

	def get_alnfile(self, gene):
		return '{}/{}.aln'.format(self.tmpdir, gene)
	def get_seqfile(self, gene):
		return '{}/{}.fa'.format(self.tmpdir, gene)
	
	def get_prefix(self, filename, suffix):
		pattern = r'(\S+){}'.format(suffix)
		filename = os.path.basename(filename)
		return re.compile(pattern).search(filename).groups()[0]
	def get_fastas(self, type):
		suffix = '.{}.fasta'.format(type)
		fastafiles = sorted(glob.glob('{}/*{}'.format(self.indir, suffix)))
		prefixs = [self.get_prefix(fasta, suffix) for fasta in fastafiles]
		return fastafiles, prefixs
	def get_files(self, suffix):
		files = sorted(glob.glob('{}/*{}'.format(self.indir, suffix)))
		prefixs = [self.get_prefix(fl, suffix) for fl in files]
		return files, prefixs
	def get_genome_id(self):
		suffix = '.fsa'
		logger.info(f"DEBUG: looking for *{suffix} in {self.indir}")
		gbfiles = sorted(glob.glob('{}/*{}'.format(self.indir, suffix)))
		logger.info(f"DEBUG: found fsa list = {gbfiles}")
    
		prefixs = []
		for fasta in gbfiles:
			logger.info(f"DEBUG: run get_prefix on {fasta}")
			p = self.get_prefix(fasta, suffix)
			prefixs.append(p)
			logger.info(f"DEBUG: get_prefix return: {p}")
		logger.info(f"DEBUG: final prefixs list = {prefixs}")

		d = {}
		for prefix, gbfile in zip(prefixs, gbfiles):
			logger.info(f"DEBUG: start reading {gbfile}, prefix={prefix}")
			for rc in SeqIO.parse(gbfile, 'fasta'):
				logger.info(f"DEBUG: first seq id: {rc.id}")
				break
			d[prefix] = rc.id
		logger.info(f"DEBUG: get_genome_id done, d={d}")
		return d

	def get_seqs(self, fastas, prefixs):
		for fasta, prefix in zip(fastas, prefixs):
			for rc in SeqIO.parse(fasta, 'fasta'):
#				if re.compile(r'-\d+$').search(rc.id):
#					continue
				attr = self.get_attributes(rc.description)
				cov_val = float(attr['cov'])
				if hasattr(self, "min_cov") and cov_val < self.min_cov:
					continue
				yield rc, prefix

	def get_attributes(self, desc):
		desc = ' '.join(desc.split(' ')[1:])
		return dict(kv.split('=') for kv in desc.split(';') if kv)
class Summary(PhyloPipeline):
	def __init__(self, indir,):
		self.indir = indir
		self.output = '{}.summry'.format(self.indir.rstrip('/'))
	def summary(self):
		from .Gff import GffGenes
		genes = []
		d_info = {}
		for gff, prefix in zip(*self.get_files('.gff3')):
			d_genes = {}
			for rc in GffGenes(gff):
				gene = rc.gene
				try: name = gene.attributes['gene']
				except KeyError:
					print('[WARN] gene {} has no cov in {}; please check'.format(gene.id, gff), file=sys.stderr)
					continue
				#print(gff, gene, gene.attributes)
				try: cov = gene.attributes['cov']
				except KeyError:
					print('[WARN] gene {} has no cov in {}; please check'.format(name, gff), file=sys.stderr)
					continue
				try: d_genes[name] += [cov]
				except KeyError: d_genes[name] = [cov]
				genes += [(rc.rna_type, name)]
			d_info[prefix] = d_genes

		f = open(self.output, 'w')
		sps = sorted(d_info.keys())
		line = ['gene', 'type'] + sps
		print('\t'.join(line), file=f)
		for type, gene in sorted(set(genes)):
			covs = [d_info[sp].get(gene, ['0']) for sp in sps]
#			print covs
			covs = ['/'.join(sorted(cov, key=lambda y:-float(y))) for cov in covs]
			line = [gene, type] + covs
			print('\t'.join(line), file=f)
		f.close()
class KaKsPipeline(PhyloPipeline):
    def __init__(self, indir, tmpdir='/dev/shm/tmp/', ncpu=20, codontable=1):
        self.indir = indir
        self.tmpdir = tmpdir
        self.ncpu = ncpu
        self.types = ['cds', 'pep']
        self.output = '{}.kaks.xls'.format(self.indir.rstrip('/'))
        # 密码子表说明：
        # 1=标准密码子（陆生植物线粒体）
        # 11=叶绿体/质体密码子
        # 4=霉菌线粒体（植物不要用）
        self.codontable = codontable

    def parse_pal2nal_paml(self, path):
        records = {}
        with open(path, "r") as f:
            lines = [line.strip() for line in f if line.strip()]
        ptr = 1
        while ptr < len(lines):
            sp = lines[ptr]
            ptr += 1
            seq = ""
            while ptr < len(lines):
                line = lines[ptr]
                if not all(c in "ATCGNatcgn-" for c in line):
                    break
                seq += line
                ptr += 1
            records[sp] = seq
        return records

    def run(self):
        files = self.cat_seqs()
        self.run_pal2nal_kaks(files)

    def run_pal2nal_kaks(self, files):
        from itertools import combinations
        aln_dir = os.path.join(self.tmpdir, "gene_aln")
        mkdirs(aln_dir)
        cds_cat = files['cds']
        pep_cat = files['pep']
        gene_dict = files['gene_dict']
        kaks_all = []

        for gene, id_list in gene_dict.items():
            logger.info(f"[KAKS] processing gene: {gene}")
            gene_pep = os.path.join(aln_dir, f"{gene}.pep.fa")
            gene_cds = os.path.join(aln_dir, f"{gene}.cds.fa")
            gene_pep_aln = os.path.join(aln_dir, f"{gene}.pep.aln")
            gene_paml = os.path.join(aln_dir, f"{gene}.paml")

            # 纯Biopython提取该基因所有物种序列
            self.extract_target_seqs(pep_cat, gene_pep, id_list)
            self.extract_target_seqs(cds_cat, gene_cds, id_list)

            # 蛋白比对 mafft
            cmd = f"mafft --auto {gene_pep} > {gene_pep_aln}"
            run_cmd(cmd, log=True)

            # pal2nal 输出 paml，规避多物种axt bug
            cmd = (f"pal2nal.pl {gene_pep_aln} {gene_cds} "
                   f"-output paml -codontable {self.codontable} > {gene_paml}")
            run_cmd(cmd, log=True)

            if not (os.path.exists(gene_paml) and os.path.getsize(gene_paml) > 0):
                logger.warning(f"[KAKS] {gene} paml file empty or missing, skip")
                continue

            # 解析paml，拆两两axt
            sp_dict = self.parse_pal2nal_paml(gene_paml)
            sp_list = list(sp_dict.keys())
            if len(sp_list) < 2:
                logger.warning(f"[KAKS] {gene} only has {len(sp_list)} species, skip")
                continue

            # 遍历所有两两组合
            for sp1, sp2 in combinations(sp_list, 2):
                s1 = sp_dict[sp1]
                s2 = sp_dict[sp2]
                # 关键：清洗竖线，防止shell把 | 识别成管道
                sp1_safe = sp1.replace("|", "_")
                sp2_safe = sp2.replace("|", "_")
                pair_tag = f"{sp1_safe}_vs_{sp2_safe}"
                pair_axt = os.path.join(aln_dir, f"{gene}.{pair_tag}.axt")
                pair_kaks = os.path.join(aln_dir, f"{gene}.{pair_tag}.kaks")

                with open(pair_axt, "w") as fw:
                    fw.write(f">{sp1} {sp2}\n")
                    fw.write(f"{s1}\n")
                    fw.write(f"{s2}\n")

                # 跑KaKs
                cmd = f"KaKs_Calculator -i {pair_axt} -o {pair_kaks} -c {self.codontable}"
                run_cmd(cmd, log=True)

                if os.path.exists(pair_kaks) and os.path.getsize(pair_kaks) > 0:
                    kaks_all.append(pair_kaks)

        # ========== 这里就是修改点 ==========
        with open(self.output, 'w') as fout:
            header_wrote = False
            for fn in kaks_all:
                with open(fn) as f:
                    for line in f:
                        if line.strip().startswith("Sequence"):
                            if not header_wrote:
                                fout.write(line)
                                header_wrote = True
                        else:
                            fout.write(line)
        logger.info(f"[KAKS] Done. Merged result: {self.output}")


    def extract_target_seqs(self, in_fa, out_fa, target_ids):
        """纯Biopython，根据id列表提取序列 + 兜底去重"""
        records = []
        seen_id = set()
        for rc in SeqIO.parse(in_fa, "fasta"):
            if rc.id in target_ids and rc.id not in seen_id:
                records.append(rc)
                seen_id.add(rc.id)
        SeqIO.write(records, out_fa, 'fasta')

    def cat_seqs(self):
        d_genes = {}
        d_files = {}
        logger.info("[DEBUG cat_seqs] enter cat_seqs, types = %s", self.types)
        for typ in self.types:
            catfasta = '{}/{}.fasta'.format(self.tmpdir, typ)
            d_files[typ] = catfasta
            fastas, prefixs = self.get_fastas(typ)
            logger.info(f"[DEBUG] typ={typ}, fastas={fastas}, prefixs={prefixs}")
            fout = open(catfasta, 'w')
            cnt = 0
            seen_sp_gene = set()  # 新增：记录已经写过的 (物种prefix,gene)
            for rc, prefix in self.get_seqs(fastas, prefixs):
                attr = self.get_attributes(rc.description)
                gene = attr['gene']
                key = (prefix, gene)
                if key in seen_sp_gene:
                    # 同一个物种同一个基因，已经存过，跳过这条转录本
                    continue
                seen_sp_gene.add(key)
                # 保留 | 分隔物种和基因
                rc.id = '{}|{}'.format(prefix, gene)
                try:
                    d_genes[gene].append(rc.id)
                except KeyError:
                    d_genes[gene] = [rc.id]
                SeqIO.write(rc, fout, 'fasta')
                cnt += 1
            fout.close()
            logger.info(f"[DEBUG] typ={typ}, wrote {cnt} records")
        logger.info(f"[DEBUG cat_seqs] total unique gene collected: {len(d_genes)}")
        logger.info(f"[DEBUG cat_seqs] gene list: {list(d_genes.keys())}")

        d_files['gene_dict'] = d_genes
        return d_files

def main():
	subcmd = sys.argv[1]
	indir = sys.argv[2]
	kargs = parse_kargs(sys.argv)
	if subcmd == 'phylo':
		try:
			root = sys.argv[3]
		except IndexError:
			root = None
		PhyloPipeline(indir, root=root, **kargs).run()
	elif subcmd == 'kaks':
		# 从kargs取出codontable，取不到就默认4（植物线粒体）
		codontable = int(kargs.pop("codontable", 1))
		KaKsPipeline(indir, codontable=codontable, **kargs).run()
	elif subcmd == 'summary':
		Summary(indir).summary()
	else:
		raise ValueError("unknown subcmd: {}".format(subcmd))

if __name__ == '__main__':
	main()
