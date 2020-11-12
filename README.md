### Installation ###
Dependencies:
+ [python 2.7](https://www.python.org/)  
    + [biopython](https://biopython.org/): quickly install by `pip install biopython<=1.76`  
    + [networkx](http://networkx.github.io/): quickly install by `pip install networkx<2.0`  
    + [lazy_property](https://github.com/jackmaney/lazy-property): quickly install by `pip install lazy-property`  
+ [hmmsearch 3.1x or 3.2x](http://hmmer.org/): compatible with HMMER3/f database format  
+ [exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) for coding genes annotation
+ [augustus](http://bioinf.uni-greifswald.de/webaugustus/) for coding genes annotation
+ [tRNAscan-SE](http://trna.ucsc.edu/software/) for tRNA genes annotation
+ [blat](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/) for rRNA genes annotation
+ [tbl2asn](https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/) output sqn file for submitting to GenBank
+ [asn2gb](https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/asn2gb/) output genbank file
+ without `-taxon` option to automatically get taxon from organism
    + [ete3](http://etetoolkit.org/) for taxonomy mapping from organism
+ `-trn_struct` option to plot tRNA secondary structure
    + [RNAplot in ViennaRNA](https://www.tbi.univie.ac.at/RNA/)  
    + ps2pdf  
+ `-draw_map` option to draw Genome Map
    + [OGDraw](https://chlorobox.mpimp-golm.mpg.de/OGDraw-Downloads.html) 
+ `-compare_map` option to draw Genome Map together with the raw genbank record
    + OGDraw  
	+ latex  
+ `-repeat` option to annotate repeat region
    + [vmatch](http://www.vmatch.de/)  
    + [trf](http://tandem.bu.edu/trf/trf.html)  

- to compare multiple annotations (Comparative.py)  
    - phylogenetics
        - [mafft](https://mafft.cbrc.jp/alignment/software/)  
        - [trimal](http://trimal.cgenomics.org/)  
        - [iqtree](http://www.iqtree.org/)  
    - KaKs
        - [ParaAT](http://bigd.big.ac.cn/tools/paraat)  
        - [KaKs_Calculator](https://bigd.big.ac.cn/tools/kaks)  

* to build custom database (Database.py)  
    * [OrthoFinder](https://github.com/davidemms/OrthoFinder)  
    * [cd-hit](http://cd-hit.org/)  

- OGAP:
```
git clone https://github.com/zhangrengang/OGAP
```

### Quick Start

##### mitochondrion genome in genbank format
```
cd OGAP/test
python ../OGAP.py Arabidopsis_thaliana-mt.gb -mt -o mt_out
```
By default, organism name will be extract from the genbank file (ORGANISM) and database will be selected by taxonomy mapping from organism, automatically.
##### mitochondrion genome in fasta format
```
python ../OGAP.py Arabidopsis_thaliana-mt.fa -mt -o mt_out -sp Arabidopsis_thaliana
```
By default, database will be selected by taxonomy mapping from organism (`-sp`), automatically.
##### mitochondrion genome with database specified (`-taxon`)
```
python ../OGAP.py Arabidopsis_thaliana-mt.fa -mt -o mt_out -sp Arabidopsis_thaliana -taxon rosids
```
##### multiple database are supported
```
python ../OGAP.py Arabidopsis_thaliana-mt.fa -mt -o mt_out -sp Arabidopsis_thaliana -taxon rosids malvids
```

##### plastid/chloroplast genome is similar but change to `-pt` mode
```
python ../OGAP.py Arabidopsis_thaliana-pt.gb -pt -o pt_out
```

### Pipeline for multiple genomes (an example)
```
for sp in Arabidopsis_thaliana Vitis_vinifera Oryza_sativa Salix_suchowensis Citrus_sinensis
do
	python ../OGAP.py genbank/$sp.gb -mt -prefix $sp -outdir re_anno &> $sp.log
done

python ../lib/Comparative.py phylo re_anno/
python ../lib/Comparative.py kaks re_anno/
```
