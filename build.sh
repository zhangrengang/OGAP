tp=mt
tx=rosids
[ $1 ] && tp=$1
[ $2 ] && tx=$2
[ $tp = mt ] && gb=~/database/RefSeq/mitochondrion/mitochondrion.*.genomic.gbff.gz
[ $tp = pt ] && gb=~/database/RefSeq/plastid/plastid.*.genomic.gbff.gz
python /share/home/nature/users/zhangrenang/ogap/makedb.py -gbfiles $gb -organ $tp -taxon $tx 2>&1 | tee build.$tp-$tx.log
