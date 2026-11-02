# Short Introduction

- **Tool purpose**: Automatic annotation pipeline for chloroplast (pt) / mitochondrial (mt) organelle genomes.
- **Applicable taxa**: Mainly plants; good performance for animal protein‑coding gene annotation.
- **Main outputs**:
  `cov_report_Ara.tsv`, `filtered.gff3`, `sequence.fasta.fsa`, `sequence.fasta.gff3`, `sequence.fasta.rna.fasta`, `filtered.cds.fasta`, `sequence.fasta.cds.fasta`, `sequence.fasta.gb`, `sequence.fasta.pep.fasta`, `stop_report.tsv`
# Quick install and start

```bash
# download
git config --global http.postBuffer 5242880000
git clone https://github.com/zhangrengang/OGAP
cd OGAP

# install
conda env create -f OGAP.yaml
conda activate OGAP

# start
cd test
python ../OGAP.py Arabidopsis_thaliana-mt.gb -mt -o mt_out
python ../OGAP.py Arabidopsis_thaliana-mt.fa -mt -o mt_out -sp Arabidopsis_thaliana
python ../OGAP.py Arabidopsis_thaliana-pt.gb -pt -o pt_out
python ../OGAP.py Arabidopsis_thaliana-mt.gb -mt -pt -o mt_out
```
---

# Installation
Core dependencies are managed by conda via `OGAP.yaml`.
Binaries including `blat`, `tbl2asn`, `asn2gb` are bundled in the `bin/` folder of this repository.

+ [python 3](https://www.python.org/) (tested: 3.12.13; minimum supported: 3.8)
    + [biopython](https://biopython.org/): >=1.85, tested on 1.87 (required if orthofinder downstream function is used)
    + [networkx](http://networkx.github.io/): >=3.6.1
+ [hmmsearch](http://hmmer.org/) (3.1 ~ 3.4, backward compatible): compatible with HMMER3/f database format
+ [exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) for coding genes annotation (2.2 ~ 2.4.0)
+ [augustus](http://bioinf.uni-greifswald.de/webaugustus/) for coding genes annotation (3.1 ~ 3.5)
> Note: Augustus training function is available, but the prfl generation module has cross-version compatibility issues caused by libstdc++. Please validate prfl files manually.
+ [tRNAscan-SE](http://trna.ucsc.edu/software/) for tRNA genes annotation
+ [blat](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/) for rRNA genes annotation (bundled in bin/)
+ [tbl2asn](https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/) output sqn file for submitting to GenBank (bundled in bin/)
+ [asn2gb](https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/asn2gb/) output genbank file (bundled in bin/)

+ without `-taxon` option to automatically get taxon from organism
    + [ete3](http://etetoolkit.org/) for taxonomy mapping from organism (optional)
+ `-trn_struct` option to plot tRNA secondary structure
    + [RNAplot in ViennaRNA](https://www.tbi.univie.ac.at/RNA/) (optional)
    + ps2pdf (optional)
+ `-draw_map` option to draw Genome Map
    + [OGDraw](https://chlorobox.mpimp-golm.mpg.de/OGDraw-Downloads.html) (optional, not in conda)
+ `-compare_map` option to draw Genome Map together with the raw genbank record
    + OGDraw (optional, not in conda)
    + latex (optional, not in conda)
+ `-repeat` option to annotate repeat region
    + [vmatch](http://www.vmatch.de/) (optional)
    + [trf](http://tandem.bu.edu/trf/trf.html) (optional)

- to compare multiple annotations (Comparative.py, optional)
    - phylogenetics
        - [mafft](https://mafft.cbrc.jp/alignment/software/)
        - [trimal](http://trimal.cgenomics.org/)
        - [iqtree](http://www.iqtree.org/)
    - KaKs
        - [ParaAT](http://bigd.big.ac.cn/tools/paraat)
        - [KaKs_Calculator](https://bigd.big.ac.cn/tools/kaks)

* to build custom database (ogap_build_db.py, optional)
    * [OrthoFinder 3.1.5](https://github.com/davidemms/OrthoFinder): NOT bundled with OGAP, install separately if homology grouping is needed
    * [cd-hit](http://cd-hit.org/): required for clustering

---

## Get OGAP
```
git clone https://github.com/zhangrengang/OGAP
```

## Quick Start
#### Mitochondrion genome in GenBank format
```bash
cd OGAP/test
python ../OGAP.py Arabidopsis_thaliana-mt.gb -mt -o mt_out
```
By default, the organism name will be extracted from the GenBank `ORGANISM` field, and the annotation database is automatically selected via taxonomy mapping.


#### Mitochondrion genome in FASTA format
```bash
python ../OGAP.py Arabidopsis_thaliana-mt.fa -mt -o mt_out -sp Arabidopsis_thaliana
```
FASTA files lack embedded organism information; the `-sp` argument is required for automatic taxonomy mapping.

#### Mitochondrion genome with manually specified taxon (-taxon)
```bash
python ../OGAP.py Arabidopsis_thaliana-mt.fa -mt -o mt_out -sp Arabidopsis_thaliana -taxon rosids
```

#### Multiple taxon databases are supported
```bash
python ../OGAP.py Arabidopsis_thaliana-mt.fa -mt -o mt_out -sp Arabidopsis_thaliana -taxon rosids malvids
```
#### Plastid / chloroplast genome (use `-pt` mode instead of `-mt`)

```bash
python ../OGAP.py Arabidopsis_thaliana-pt.gb -pt -o pt_out
```
#### Annotate sequences containing both mitochondrial and chloroplast-derived fragments
```bash
python ../OGAP.py Arabidopsis_thaliana-mixed.gb -mt -pt -o mt_out
```

> ⚠️ Warning: `-mt -pt` combined mode searches both mitochondrial and plastid gene models.
> This is only recommended for mixed organelle contigs. Do NOT use this for pure mitochondrial or pure plastid sequences, as it increases runtime and false positives.


### Pipeline for multiple genomes (optional downstream workflow)
```bash
for sp in Arabidopsis_thaliana Vitis_vinifera Oryza_sativa Salix_suchowensis Citrus_sinensis
do
    python ../OGAP.py genbank/$sp.gb -mt -prefix $sp -outdir re_anno &> $sp.log
done

# Comparative analysis (optional, requires extra comparative dependencies)
python ../lib/Comparative.py summary re_anno/
python ../lib/Comparative.py phylo re_anno/
python ../lib/Comparative.py kaks re_anno/
```

#### Core parameter description

| Parameter | Description |
|---|---|
| `-pt` / `-mt` | Specify organelle type; alias for `--ogtype pt/mt` |
| `-sp` | Species name, matches AUGUSTUS species model |
| `-taxon` | Taxonomic group, used for internal pipeline filtering and scoring |
| `--no-rrna-alias` | Disable rRNA alias normalization, revert to original annotation baseline |
| `-min_rrn_score=0.6` | HMM score threshold for rRNA; 0.6 commonly‑used for mitochondria, adjust per taxon |
| `-trn-opts='-O -X 5'` | Arguments passed to tRNAscan‑SE; commonly‑used settings for mitochondrial mode |
| `-o` | Output directory path |

- **Note:** 
The & runs jobs in parallel. Adjust the loop carefully on shared computing clusters to avoid overloading resources.
Comparative analysis functions are optional and require additional tools: mafft, trimal, iqtree, ParaAT, KaKs_Calculator.
- **General note for all test examples:**
The Arabidopsis sequence files are not bundled in the repository. These examples show syntax only.
---

### Build custom organelle database
- **ogap_build_db.py, optional**
  
OGAP allows users to build custom mitochondrial or plastid databases for specialized taxonomic groups.
This step is **NOT required for normal annotation**; you only need it when the default built-in databases are not suitable for your target species.

Required dependency: `cd-hit`
Optional dependency: OrthoFinder 3.1.5 (only for homology grouping during database construction, not bundled)

##### Basic usage (build a mitochondrial database)
```bash
# pt
python3 ./OGAP_master/ogap_build_db.py \
  -organ pt \
  -gbfiles ./*.gb \
  -custom <taxon_name> \
  --dbroot-dir ./db/pt \
  --orthofinder-exe ./orthofinder-3.1.5/orthofinder.py \
  --overwrite

# mt
python3 ./OGAP_master/ogap_build_db.py \
  -organ mt \
  -gbfiles ./*.gb \
  -custom <taxon_name> \
  --dbroot-dir ./db/mt \
  --orthofinder-exe ./orthofinder-3.1.5/orthofinder.py \
  --overwrite
```

- `-organ pt/mt`: build plastid or mitochondrial database
- `-gbfiles ./*.gb`: input a set of GenBank reference files
- `-custom <taxon_name>`: name for your custom taxon (used later with OGAP `-taxon`)
- `--dbroot-dir`: output directory for the built database
- `--orthofinder-exe`: - `--orthofinder-exe`: safe fallback option; manually supply the full path to orthofinder.py, works regardless of whether OrthoFinder is in PATH
- `--overwrite`: overwrite existing database files

