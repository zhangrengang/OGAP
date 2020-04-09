OGAP: Organelle Genome Annotation Pipeline

### Installation ###
Dependencies:
+ [python 2.7](https://www.python.org/)  
    + [biopython](https://biopython.org/): quickly install by `pip install biopython`  
    + [networkx](http://networkx.github.io/)  
    + [lazy_property](https://github.com/jackmaney/lazy-property)  
+ [hmmsearch 3.1x or 3.2x](http://hmmer.org/): compatible with HMMER3/f database format  
+ [exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
+ [augustus](http://bioinf.uni-greifswald.de/webaugustus/)
+ [tRNAscan-SE](http://trna.ucsc.edu/software/)
+ [blat](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/)
+ [tbl2asn](https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/)
+ [asn2gb](https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/asn2gb/)
+ without `-taxon` option  
    + [ete3](http://etetoolkit.org/) for taxonomy mapping from organism
+ `-trn_struct` option  
    + [RNAplot in ViennaRNA](https://www.tbi.univie.ac.at/RNA/)  
    + ps2pdf  
+ `-drawgenemap` option  
    + [drawgenemap](https://chlorobox.mpimp-golm.mpg.de/OGDraw-Downloads.html)  
+ `-repeat` option  
    + [vmatch](http://www.vmatch.de/)  
    + [trf](http://tandem.bu.edu/trf/trf.html)  

- to compare multiple genomes (Comparative.py)  
    - phylogenetics
        - [mafft](https://mafft.cbrc.jp/alignment/software/)  
        - [trimal](http://trimal.cgenomics.org/)  
        - [iqtree](http://www.iqtree.org/)  
    - KaKs
        - [ParaAT](http://bigd.big.ac.cn/tools/paraat)  
        - [KaKs_Calculator](https://bigd.big.ac.cn/tools/kaks)  

*   to make database (Database.py)  
    *   [OrthoFinder](https://github.com/davidemms/OrthoFinder)  
    *   [cd-hit](http://cd-hit.org/)  

