OG=$1
f1=$2
f2=$3

cd "AUTesting/"$OG

iqtree -s $f1 -z $f2 -m LG+G -zb 10000 -n 0 -au