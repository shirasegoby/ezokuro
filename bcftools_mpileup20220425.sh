#!/bin/sh
READDIR=/Users/shotarohirase/Desktop/エゾクロ解析re/
REF=/Users/shotarohirase/Desktop/エゾクロ解析re/reference/Haliotis.fa

BAMS=""
for bam in `ls ./ezokuro_realigned_bam/*.bam.realigned.bam.uniq.bam`
do
  BAMS=${BAMS}" ${bam} "
done

/Users/shotarohirase/Desktop/Genomic_tools/bcftools-1.6/bcftools mpileup -Ou --annotate FORMAT/DP -f $REF $BAMS | /Users/shotarohirase/Desktop/Genomic_tools/bcftools-1.6/bcftools call --threads 10 -mv -Oz -o ezokuro.vcf.gz


#/Users/shotarohirase/Desktop/Genomic_tools/vcftools/src/cpp/vcftools --gzvcf family_F.vcf.gz --minDP 20 --min-meanDP 20 --max-missing-count 10 --out family_F_DP20 --recode