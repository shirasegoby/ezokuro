#!/bin/sh
REF=Haliotis.fa

BAMS=""
for bam in `ls ./ezokuro_realigned_bam/*.bam.realigned.bam.uniq.bam`
do
  BAMS=${BAMS}" ${bam} "
done

bcftools mpileup -Ou --annotate FORMAT/DP -f $REF $BAMS | bcftools call --threads 10 -m -Oz -o ezokuro_all_snp.vcf.gz

