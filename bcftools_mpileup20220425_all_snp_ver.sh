#!/bin/sh
READDIR=/Users/shotarohirase/Desktop/エゾクロ解析re/
REF=/Users/shotarohirase/Desktop/エゾクロ解析re/reference/Haliotis.fa

BAMS=""
for bam in `ls ./ezokuro_realigned_bam/*.bam.realigned.bam.uniq.bam`
do
  BAMS=${BAMS}" ${bam} "
done

/Users/shotarohirase/Desktop/Genomic_tools/bcftools-1.6/bcftools mpileup -Ou --annotate FORMAT/DP -f $REF $BAMS | /Users/shotarohirase/Desktop/Genomic_tools/bcftools-1.6/bcftools call --threads 10 -m -Oz -o ezokuro_all_snp.vcf.gz

