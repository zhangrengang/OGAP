tp=mt
tx=rosids
[ $tp = mt ] && gb=~/database/RefSeq/mitochondrion/mitochondrion.*.genomic.gbff.gz
[ $tp = pt ] && gb=~/database/RefSeq/plastid/plastid.*.genomic.gbff.gz
python /share/home/nature/users/zhangrenang/ogap/lib/makedb.py -gbfiles $gb -organ $tp -taxon $tx &> build.$tp-$tx.log
