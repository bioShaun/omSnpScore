#!/bin/bash

bedtools="/data/software/bedtools/bedtools-2.27.1/bin/bedtools"

varScore_csv=$1    # mutant_0.3_0.7_wild_0.3_0.7.snp_num.window.w20.s5.var.score.csv
qtlseqr_ed_csv=$2  # qtlseqr.10000000.RIL.ed.0.3.5.csv
snp_freq_csv=$3
out_prefix=$4

# Data prepare

sed 's/,/\t/g' $varScore_csv | grep "^chr"  > ${out_prefix}/varscore.bed
$bedtools map -o max -c 4 -a /data/scripts/circos/1m.window -b ${out_prefix}/varscore.bed |awk '$4!="."'  > ${out_prefix}/varscore.circos.bed

cut -f1,2 -d"," $snp_freq_csv |sed 's/,/\t/' |awk -v OFS="\t" '{print $1,$2-1,$2,1}'|grep "chr" > ${out_prefix}/snp_density.bed
$bedtools map -o sum -c 4 -a /data/scripts/circos/1m.window -b ${out_prefix}/snp_density.bed |awk '$4!="."'  > ${out_prefix}/snp_density.circos.bed

cut -f1,2,28 -d"," $qtlseqr_ed_csv |sed 's/,/\t/g'|grep "^chr"|awk -v OFS="\t" '{print $1,$2,$2+1,$3}' > ${out_prefix}/ed.bed
$bedtools map -o max -c 4 -a /data/scripts/circos/1m.window -b ${out_prefix}/ed.bed |awk '$4!="."'  > ${out_prefix}/ed.circos.bed

cut -f1,2,19 -d"," $qtlseqr_ed_csv |sed 's/,/\t/g'|grep "^chr"|awk -v OFS="\t" '{print $1,$2,$2+1,$3}' > ${out_prefix}/Gprime.bed
$bedtools map -o max -c 4 -a /data/scripts/circos/1m.window -b ${out_prefix}/Gprime.bed |awk '$4!="."'  > ${out_prefix}/Gprime.circos.bed

cut -f1,2,15 -d"," $qtlseqr_ed_csv |sed 's/,/\t/g'|grep "^chr"|awk -v OFS="\t" '{print $1,$2,$2+1,$3}' > ${out_prefix}/qtlseq.bed
$bedtools map -o mean -c 4 -a /data/scripts/circos/10m.window -b ${out_prefix}/qtlseq.bed |awk '$4!="."'  > ${out_prefix}/qtlseq.circos.bed

#cut -f1,2,25 -d"," $qtlseqr_ed_csv |sed 's/,/\t/g'|grep "^chr"|awk -v OFS="\t" '{print $1,$2,$2+1,$3}' > ${out_prefix}/CI_95.bed
#$bedtools map -o mean -c 4 -a /data/scripts/circos/10m.window -b ${out_prefix}/CI_95.bed |awk '$4!="."'  > ${out_prefix}/CI_95.circos.bed

#awk -v OFS="\t" '{print $1,$2,$3,-1-$4}' ${out_prefix}/CI_95.circos.bed > ${out_prefix}/CI_95.neg.circos.bed
#awk -v OFS="\t" '{print $1,$2,$3,1+$4}' ${out_prefix}/CI_95.circos.bed > ${out_prefix}/CI_95.pos.circos.bed

#cut -f1,2,26 -d"," $qtlseqr_ed_csv |sed 's/,/\t/g'|grep "^chr"|awk -v OFS="\t" '{print $1,$2,$2+1,$3}' > ${out_prefix}/CI_99.bed

#$bedtools map -o mean -c 4 -a /data/scripts/circos/10m.window -b ${out_prefix}/CI_99.bed |awk '$4!="."'  > ${out_prefix}/CI_99.circos.bed
#awk -v OFS="\t" '{print $1,$2,$3,-1-$4}' ${out_prefix}/CI_99.circos.bed > ${out_prefix}/CI_99.neg.circos.bed
#awk -v OFS="\t" '{print $1,$2,$3,1+$4}' ${out_prefix}/CI_99.circos.bed > ${out_prefix}/CI_99.pos.circos.bed


# run circos
/home/zxchen/miniconda3/bin/circos -conf ${out_prefix}/circos.conf
