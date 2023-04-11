# this script requires two others:
EIGENSTRAT_TO_VCF=~/bin/eigenstrat2vcf_RT.py # this is a custom script converting eigenstrat to vcf formatted data
COMPUTED=~/bin/dating-master/bin/computed # this is the computed algorithm


echo "============================================"
echo "COMPUTED"
echo "============================================"

NCHROM=$1
EIGENPREFIX=$2

echo $EIGENPREFIX

cut -f3 -d" " $EIGENPREFIX.ind >$EIGENPREFIX.pop
gunzip -f $EIGENPREFIX.*gz
cp $EIGENPREFIX.snp $EIGENPREFIX.tmp
awk '{print $1, $2, $3, $4, "A", "T" }' $EIGENPREFIX.tmp >$EIGENPREFIX.snp
python $EIGENSTRAT_TO_VCF -r $EIGENPREFIX -p $EIGENPREFIX.pop --anc0 >$EIGENPREFIX.vcf
echo "#IID SEX POP" >$EIGENPREFIX.ind2
cat $EIGENPREFIX.ind >>$EIGENPREFIX.ind2
# ascertainment
plink2 --vcf $EIGENPREFIX.vcf --freq \
--autosome-num $NCHROM \
--pheno $EIGENPREFIX.ind2 \
--loop-cats POP \
--out $EIGENPREFIX
paste <(cut -f5 $EIGENPREFIX.YRI.afreq) <(cut -f5 $EIGENPREFIX.CEU.afreq) <(cut -f5 $EIGENPREFIX.Vindija.afreq) | tail -n +2 >$EIGENPREFIX.afreq
paste $EIGENPREFIX.afreq <(grep -vP "^#" $EIGENPREFIX.vcf) >$EIGENPREFIX.vcf2
awk '$3 > 0 { print $0 }' $EIGENPREFIX.vcf2 | awk '$2 > 0 { print $0 }' | awk '$2 < 0.1 { print $0 }' >$EIGENPREFIX.vcf3
cat <(grep -P "^#" $EIGENPREFIX.vcf) <(cut -f4- $EIGENPREFIX.vcf3) >$EIGENPREFIX.vcf4
#
python $EIGENSTRAT_TO_VCF -v $EIGENPREFIX.vcf4 -o $EIGENPREFIX.asc
awk '{print $1, $2, $4*1e-8, $4, "A", "T" }' $EIGENPREFIX.asc.snp >$EIGENPREFIX.asc.gpos.snp
cat $EIGENPREFIX.asc.geno | cut -c 51-100 >$EIGENPREFIX.asc.s.geno
head -n100 $EIGENPREFIX.asc.ind | tail -n50 >$EIGENPREFIX.asc.s.ind

# ComputeD
echo "genotypename: $EIGENPREFIX.asc.s.geno" >$EIGENPREFIX.dpar
echo "snpname: $EIGENPREFIX.asc.gpos.snp" >>$EIGENPREFIX.dpar
echo "indivname: $EIGENPREFIX.asc.s.ind" >>$EIGENPREFIX.dpar
echo "output: $EIGENPREFIX.computed.txt" >>$EIGENPREFIX.dpar
echo "binsize: 0.00001" >>$EIGENPREFIX.dpar
echo "maxdis: 0.01" >>$EIGENPREFIX.dpar

$COMPUTED -p $EIGENPREFIX.dpar

rm $EIGENPREFIX.vcf* $EIGENPREFIX.tmp $EIGENPREFIX.pop $EIGENPREFIX.*afreq $EIGENPREFIX.ind2 $EIGENPREFIX.log
rm $EIGENPREFIX.*geno $EIGENPREFIX.*geno1 $EIGENPREFIX.*snp $EIGENPREFIX.*ind

#___
