# Short Introduction

- **Tool purpose**: Automatic annotation pipeline for chloroplast (pt) / mitochondrial (mt) organelle genomes.
- **Applicable taxa**: Mainly plants; good performance for animal protein‑coding gene annotation.

# Quick install and start
> **Recommendation:**
 Use a dedicated conda environment for OGAP.
This section covers the standard workflow for **organelle genome annotation (core function).**
We use the lighter `OGAP.yaml` environment here, which does not include OrthoFinder.
This setup avoids dependency version conflicts for tools such as AUGUSTUS, HMMER.
OrthoFinder is only required for custom database building (`makedb.py`), and OrthoFinder itself also recommends isolated environments for stable execution.

```bash
# download
git config --global http.postBuffer 5242880000
git clone https://github.com/zhangrengang/OGAP
cd OGAP

# install
conda env create -f ./config/OGAP.yaml
conda activate OGAP

# start
cd test
# Input: GenBank  | Mode: -mt (mitochondrion)
python ../OGAP.py Arabidopsis_thaliana-mt.gb -mt -o mt_out
# Input: FASTA    | Mode: -mt (mitochondrion, -sp REQUIRED)
python ../OGAP.py Arabidopsis_thaliana-mt.fa -mt -o mt_out -sp Arabidopsis_thaliana
# Input: GenBank  | Mode: -pt (plastid)
python ../OGAP.py Arabidopsis_thaliana-pt.gb -pt -o pt_out
# Input: GenBank  | Mode: -mt -pt (combined search)
python ../OGAP.py Arabidopsis_thaliana-mt.gb -mt -pt -o mt_out
```
## Input requirements
- GenBank format (`.gb` / `.gbff`): Complete GenBank flatfile with `ORGANISM` field; preferred input type.
- FASTA format (.fa / .fasta): Only raw nucleotide sequence. `-sp <species_name>` is required. The name will first be used for taxonomic lookup. If lookup fails, it serves only as a label placeholder.
- Mixed contig assemblies are supported. The `-mt -pt` combined mode is designed for this scenario, to detect plastid fragments integrated into mitochondrial sequences.
---

# Installation
Core dependencies are managed by conda via `OGAP.yaml`.
Binaries including `blat`, `tbl2asn`, `asn2gb` are bundled in the `bin/` folder of this repository.

+ [python 3](https://www.python.org/) (tested: 3.12.13; minimum supported: 3.8)
    + [biopython](https://biopython.org/): >=1.85, tested on 1.87 (required if orthofinder downstream function is used)
    + [networkx](http://networkx.github.io/): >=3.6.1
    + [numpy](https://numpy.org/) >=2.3.2
    + [scipy](https://scipy.org/) >=1.16
    + [rich](https://rich.readthedocs.io/) >=14.1.0
    + [scikit-learn](https://scikit-learn.org/) >=1.7.1
+ [hmmsearch](http://hmmer.org/) (3.1 ~ 3.4, backward compatible): compatible with HMMER3/f database format
+ [exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) for coding genes annotation (2.2 ~ 2.4.0)
+ [augustus](http://bioinf.uni-greifswald.de/webaugustus/) for coding genes annotation (3.1 ~ 3.5)
> Note: Augustus training function is available, but the prfl generation module has cross-version compatibility issues caused by libstdc++. Please validate prfl files manually.
+ [tRNAscan-SE](http://trna.ucsc.edu/software/) for tRNA genes annotation
+ [blat](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/) for rRNA genes annotation (bundled in bin/)
+ [tbl2asn](https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/) output sqn file for submitting to GenBank (bundled in bin/)
+ [asn2gb](https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/asn2gb/) output genbank file (bundled in bin/)
+ [MAFFT](https://mafft.cbrc.jp/alignment/software/) multiple sequence alignment
+ [Trimal](http://trimal.cgenomics.org/) alignment trimming
+ [IQ-TREE](http://www.iqtree.org/) phylogenetic tree inference
+ [FAMSA](https://github.com/GarrettGunn/FAMSA) multiple sequence alignment
+ [FastTree](http://www.microbesonline.org/fasttree/) >=2.1.11, approximate phylogenetic tree
+ [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) local sequence alignment
+ [pal2nal](http://www.bork.embl.de/pal2nal/) convert protein alignment to codon alignment
+ [KaKs_Calculator2](https://bigd.big.ac.cn/tools/kaks) calculate Ka/Ks values

Python packages installed via pip:
+ [future](https://python-future.org/)
+ [ete3](http://etetoolkit.org/) >=3.1.3, <4.0 (for taxonomy mapping)
+ [lazy-property](https://pypi.org/project/lazy-property/)
+ [xopen](https://pypi.org/project/xopen/)

+ without `-taxon` option to automatically get taxon from organism
    + ete3 (optional, already listed above)
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
        - MAFFT
        - Trimal
        - IQ-TREE
    - Ka/Ks calculation
        - pal2nal
        - KaKs_Calculator2

- to build custom database (makedb.py, optional)
    - [OrthoFinder 3.1.5](https://github.com/davidemms/OrthoFinder): NOT bundled with OGAP, install separately if homology grouping is needed
    - [cd-hit](http://cd-hit.org/): required for sequence clustering
---

## Get OGAP
Clone the repository if you prefer manual environment setup.
```bash
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
Use `-taxon` to explicitly select reference database(s) from the `./db` folder, overriding the automatic taxon-based database selection.
The name can be either an OGAP built-in database or a custom database built by `makedb.py`.

#### Multiple taxon databases are supported
```bash
python ../OGAP.py Arabidopsis_thaliana-mt.fa -mt -o mt_out -sp Arabidopsis_thaliana -taxon rosids malvids
```
You can supply **multiple database names** to `-taxon`.
OGAP will query all these databases together for annotation.
These names correspond to database entries under your `./db` directory, either built-in OGAP databases or custom databases created with `makedb.py`.

#### Plastid / chloroplast genome (use `-pt` mode instead of `-mt`)
```bash
# pt: plastid (chloroplast) genome annotation
python ../OGAP.py Arabidopsis_thaliana-pt.gb -pt -o pt_out

# mt: mitochondrial genome annotation
python ../OGAP.py Arabidopsis_thaliana-mt.gb -mt -o mt_out
```


#### Annotate sequences containing both mitochondrial and chloroplast-derived fragments
```bash
python ../OGAP.py Arabidopsis_thaliana-mixed.gb -mt -pt -o mt_out
```
Combined `-mt -pt` mode searches against both mitochondrial and plastid gene models.
It can be used to identify plastid genes integrated into mitochondrial contigs, a common feature in plant mitochondrial genomes.
This mode also helps assess the confidence of mitochondrial gene annotations.

#### Multiple-genome pipeline (optional downstream workflow)
`-prefix` assigns a species-specific prefix for output files, so annotation files from different species will not overwrite each other.
```bash
# Batch annotation for multiple organelle genomes
for sp in Arabidopsis_thaliana Vitis_vinifera Oryza_sativa
do
    python ../OGAP.py ../test/genbank/$sp.gb \
    -mt -prefix $sp \
    -o re_anno \
    &> $sp.log
done
```
>**Note:**
>- This loop runs OGAP annotation for each species separately. All results are written into the shared re_anno/ folder.
>- The **`` & ``** submits jobs in parallel. On shared clusters, control the number of concurrent jobs to avoid resource overload.
>- This step only generates individual annotation files for each species, and does not perform comparative analysis.

#### Comparative analysis (optional)
Run these three subcommands on the folder containing all species annotation results (`re_anno/`) to get orthogroup summary, multiple-sequence alignments, species tree and Ka/Ks estimates.
```bash
python ../compare.py summary re_anno/
python ../compare.py phylo re_anno/
python ../compare.py kaks re_anno/
```
>**Note:**
>- These subcommands belong to compare.py, separate from the main OGAP.py annotation parameters.
>- Extra dependencies (mafft, trimal, iqtree, pal2nal, KaKs_Calculator) are required for this module.
>- pal2nal is only required for the kaks step.

##### **General note for all test examples:**
The Arabidopsis sequence files are not bundled in the repository. These examples show syntax only.

### Core parameter description
| Parameter | Description |
| ---- | ---- |
| `-pt` / `-mt` | Enable chloroplast (`-pt`) or mitochondrial (`-mt`) annotation mode. Combine `-pt -mt` for mixed mode. |
| `-sp` | Species label. This label is used for taxonomic mapping (**ete3**) and downstream comparative analysis. For FASTA input mode, a valid scientific name is required for taxon recognition; if unrecognized, it only serves as a custom identifier. |
| `-taxon` | Name(s) of reference database(s) under `./db`; can be built-in OGAP databases or custom databases built by `makedb.py`. Multiple database names separated by space are supported. |
| `-prefix` | Prefix for output files when running batch annotation. |
| `-o` / `--outdir` |Output directory path; all results will be saved here.|
| `-min_rrn_score=0.6` | HMM score threshold for rRNA identification; 0.6 is a typical cutoff for animal mitochondrial data, adjust according to your taxon. Usually not modified for plant organelle annotation. |
| `-trn-opts='-O -X 5'` | Extra arguments passed to tRNAscan‑SE; this preset is intended for animal mitochondrial tRNA prediction, and generally not needed for plants. |

#### Main outputs:
Output files share the same prefix `<name>`, which is derived from the basename of your input FASTA or GenBank file.

`<name>.fsa`, `<name>.gff3`, `<name>.rna.fasta`, `<name>.cds.fasta`, `<name>.gb`, `<name>.pep.fasta`, `<name>.tbl`, `<name>.sqn`

- `.tbl` / `.sqn`: NCBI submission-ready feature table and sqn files
- `.fsa`: cleaned input genomic sequence
- `.gff3`: final annotation in GFF3 format
- `.rna.fasta`: rRNA and tRNA nucleotide sequences
- `.cds.fasta`: CDS (coding DNA sequences)
- `.pep.fasta`: translated protein sequences
- `.gb`: GenBank format annotation file

> **Note:** For batch multiple-genome annotation, all species results are placed in the shared output directory. Each species uses its own `<name>` (controlled by `-prefix`
) for file prefixes to avoid overwriting.

#### Comparative analysis outputs (../compare.py)
These files use the user-defined project prefix `<project>` (usually the name you assign to the comparison output directory), and this prefix is independent of the single-genome `<name>` prefix.

`<project>.summry`, `<project>.kaks.xls`, `<project>-cds.aln`, `<project>-cds.aln.trimal`, `<project>-cds.aln.trimal.iqtree`, `<project>-cds.aln.trimal.treefile`

- `<project>.summry`: Summary table of gene presence/absence and copy number across all species (stored in the top-level project directory)
- `<project>.kaks.xls`: Ka/Ks calculation results for orthologous gene pairs (stored in the top-level project directory)
- `<project>-cds.aln`: Raw codon alignment of orthogroups
- `<project>-cds.aln.trimal`: Trimmed codon alignment after Trimal filtering
- `<project>-cds.aln.trimal.iqtree`: IQ-TREE log file
- `<project>-cds.aln.trimal.treefile`: Final species tree (Newick format)
- `kaks_origin.log`: Log for Ka/Ks calculation pipeline
- `phylo_origin.log`: Log for phylogenetic reconstruction pipeline

> Note: Species-level annotation files (`<species>.gb`, `<species>.gff3`, `<species>.cds.fasta`, etc.) are also generated inside the comparison directory. These are individual annotation outputs for each taxon, not comparative results.
---


### Build custom organelle database
- **makedb.py(optional)**

OGAP allows users to build custom mitochondrial or plastid databases for specialized taxonomic groups.
This step is **NOT required for normal annotation**; you only need it when the default built-in databases are not suitable for your target species.

Required dependency: `cd-hit`
Optional dependency: OrthoFinder 3.1.5 (used for homology grouping during database construction; not bundled with conda environment)

##### **Environment preparation**
To run `makedb.py` with OrthoFinder homology grouping, use the all-in-one conda environment and download the precompiled OrthoFinder release:
```bash
# Environment installation
conda deactivate
conda env create -f ./config/orthofinder_build.yaml
conda activate ogap-build

# OrthoFinder download & setup
wget https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.5/orthofinder-linux-intel-3.1.5.tar.gz
tar -zxvf orthofinder-linux-intel-3.1.5.tar.gz
```
>**Note:** We use the precompiled Linux Intel build orthofinder-linux-intel-3.1.5.tar.gz. The smaller orthofinder-3.1.5.tar.gz is source code and not recommended, as extra compilation is required.

#### Basic usage
Three usage modes:

**Mode 1: Build database from all input GenBank files (default path, recommended)**
All input GenBank records are used directly, without taxon filtering.
`-taxon` defines your custom database name, which will be invoked later in OGAP annotation with `-taxon`.
`--dbroot-dir` is omitted here; the database will be saved automatically to `OGAP_master/db/{mt,pt}/<custom_db_name>`.
```bash
# pt
python3 ./OGAP_master/makedb.py \
  -organ pt \
  -gbfiles ./*.gb \
  -taxon <custom_db_name> \
  --orthofinder-exe ./orthofinder-linux-intel-3.1.5/orthofinder.py \
  --overwrite

# mt
python3 ./OGAP_master/makedb.py \
  -organ mt \
  -gbfiles ./*.gb \
  -taxon <custom_db_name> \
  --orthofinder-exe ./orthofinder-linux-intel-3.1.5/orthofinder.py \
  --overwrite
```
- `-taxon`: name of your custom database, used for calling this database in OGAP annotation
- All provided .gb records are included; no taxon filtering is performed

**Mode 2: Build database with taxon filtering (default path, recommended)**
Filter input GenBank records and retain only sequences belonging to the target taxon, then build the database.
The database name is the taxon name you supplied to --filter-taxon.
`--dbroot-dir` is omitted here; the database will be saved automatically to `OGAP_master/db/{mt,pt}/<taxon_name>`.
```bash
# pt
python3 ./OGAP_master/makedb.py \
  -organ pt \
  -gbfiles ./*.gb \
  --filter-taxon <taxon_to_filter> \
  --orthofinder-exe ./orthofinder-linux-intel-3.1.5/orthofinder.py \
  --overwrite

# mt
python3 ./OGAP_master/makedb.py \
  -organ mt \
  -gbfiles ./*.gb \
  --filter-taxon <taxon_to_filter> \
  --orthofinder-exe ./orthofinder-linux-intel-3.1.5/orthofinder.py \
  --overwrite
```
- `--filter-taxon`: filter input GenBank records, keep only sequences assigned to this taxon
- For batch building databases for many taxa, use `--filter-taxon-list` (one taxon name per line in the text file)
**Important:** `--filter-taxon` / `--filter-taxon-list` cannot be used together with -taxon.

**Mode 3: Advanced — custom output root directory**
Add `--dbroot-dir` to specify your own output root folder, instead of the default `OGAP_master/db/`.
This flag works for both unfiltered (`-taxon`) and filtered (`--filter-taxon`) workflows.
```bash
# Example: unfiltered custom database with custom root
python3 ./OGAP_master/makedb.py \
  -organ mt \
  -gbfiles ./*.gb \
  -taxon <custom_db_name> \
  --dbroot-dir /your/custom/path \
  --orthofinder-exe ./orthofinder-linux-intel-3.1.5/orthofinder.py \
  --overwrite
```
`--dbroot-dir`: custom root directory for database output (for advanced users)
>**Note:** OrthoFinder is required for homology grouping and classification profile generation. If you omit `--orthofinder-exe`, this homology grouping step will be skipped.
>After construction, OGAP can automatically discover databases under the default `OGAP_master/db/pt` or `OGAP_master/db/mtwhen `using `-taxon`.


#### **Parameter overview**
- `-organ pt/mt`: build plastid or mitochondrial database
- `-gbfiles ./*.gb`: input a set of GenBank flat files (GBFF). Both `.gb` and `.gbff` are supported. Files are parsed as GenBank format regardless of filename suffix, provided the content is valid GBFF.
- `-taxon`: unique name for your custom database, used to call this database in OGAP annotation. Mutually exclusive with `--filter-taxon` / `--filter-taxon-list`.
- `--filter-taxon`: keep only GenBank records matching this taxon name. Mutually exclusive with `-taxon`.
- `--filter-taxon-list`: path to a text file, one taxon name per line for batch filtering and database construction. Mutually exclusive with `-taxon`.
- `--dbroot-dir`: custom root directory to store the built database. If omitted, defaults to `OGAP_master/db/`.
- `--orthofinder-exe`: manually provide the full path to `orthofinder.py`. This is a safe option and does not require OrthoFinder to be in your PATH
- `--overwrite`: overwrite existing database files if they already exist


#### Outputs of `makedb` database construction
The database directory is named `<db_name>`. All generated files and folders are listed below:
<db_name>/
- `agst.tgz`  : Archive of AGST training files
-  `cds.info`  : Metadata table for CDS sequences
- `genome.info` : Basic information of source genomes
- `hmms.tgz` : Compressed archive of gene HMM profiles
- `msas.tgz` : Compressed archive of multiple sequence alignments
- `prfl.tgz`  : Archive of auxiliary profile files
- `protein.info`  : Metadata table for protein sequences
- `rna.info`  : Metadata table for rRNA and tRNA genes
- `seqs.tgz` :  Compressed archive of raw CDS and protein sequences
- `species.treefile` : Species tree inferred from orthogroups
- `taxonomy.info` : Plain-text taxonomy summary
- `taxonomy.json.gz` : Compressed JSON taxonomy mapping file
- `version ` : Database version marker
