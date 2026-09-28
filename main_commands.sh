#!/bin/bash

# Population genomic analysis pipeline
#
# Required software:
#   bcftools
#   vcftools
#   PLINK
#   ADMIXTURE
#   dDocent filtering scripts
#   genomics_general
#   easySFS
#   bgzip/tabix
#
# Additional R/Perl scripts are required for some analyses.


########################################
# 1. SNP calling
########################################

# Run bcftools_mpileup.sh before starting this script.


########################################
# 2. SNP filtering
########################################

vcftools \
    --gzvcf ezokuro.vcf.gz \
    --max-alleles 2 \
    --minDP 10 \
    --min-meanDP 10 \
    --minQ 20 \
    --remove-indels \
    --out ezokuro_filtered \
    --recode


########################################
# 3. dDocent filtering
########################################

sh dDocent/scripts/pop_missing_filter.sh \
    ezokuro_filtered.recode.vcf \
    popmap.txt \
    0.2 \
    2 \
    ezokuro_filtered.recode.dDocent


########################################
# 4. Convert VCF to PLINK format
########################################

vcftools \
    --vcf ezokuro_filtered.recode.dDocent.recode.rename2.vcf \
    --plink \
    --out ezokuro_filtered.recode.dDocent.recode.rename2

sort ezokuro_filtered.recode.dDocent.recode.rename2.ped \
    > ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped


########################################
# 5. SNP filtering for ADMIXTURE
########################################

vcftools \
    --vcf ezokuro_filtered.recode.dDocent.recode.rename2.vcf \
    --thin 1000 \
    --plink \
    --out ezokuro_filtered.recode.dDocent.recode.rename2.thin1000


########################################
# 6. Convert PLINK genotypes to PED12
########################################

plink \
    --file ezokuro_filtered.recode.dDocent.recode.rename2.thin1000 \
    --recode12 \
    --out ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12

sort ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12.ped \
    > ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12.sorted.ped


########################################
# 7. ADMIXTURE analysis
########################################

for K in {1..9}
do
    admixture \
        --cv \
        ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12.sorted.ped \
        "$K" \
        | tee "Admixture_log${K}.out"
done

grep -h CV ./*log*.out | sed 's/:/\t/g' > CV_error

# Plot ADMIXTURE results using:
# Pophelper_ADMIXTURE_plot2.R


########################################
# 8. OutFLANK analysis
########################################

mkdir -p outflank

cp ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped outflank/
cp ezokuro_filtered.recode.dDocent.recode.rename2.map outflank/

cd outflank || exit 1


# Modify MAP file

sed -i.bak 's/:/\t/g' \
    ezokuro_filtered.recode.dDocent.recode.rename2.map

rm -f ezokuro_filtered.recode.dDocent.recode.rename2.map.bak

SNP_number=$(
    wc -l < ezokuro_filtered.recode.dDocent.recode.rename2.map
)

seq -s $'\nSNP_' 0 "$SNP_number" > tmp.txt

SNP_number2=$((SNP_number + 1))

sed '1d' tmp.txt > tmp2.txt
sed "${SNP_number2}d" tmp2.txt > tmp3.txt

cut -f 4,5 \
    ezokuro_filtered.recode.dDocent.recode.rename2.map \
    > tmp4.txt

cut -f 2 \
    ezokuro_filtered.recode.dDocent.recode.rename2.map \
    | paste - tmp3.txt tmp4.txt \
    | sed 's/\t/ /g' \
    > ezokuro_filtered.recode.dDocent.recode.rename2.sorted.map

rm -f tmp*.txt


# Convert to additive genotype format for OutFLANK

plink \
    --file ezokuro_filtered.recode.dDocent.recode.rename2.sorted \
    --recodeA \
    --out ezokuro_filtered.recode.dDocent.recode.rename2.sorted.OUTflank2 \
    --allow-extra-chr

# Run:
# outflank_ezo_kuro.R

cd ..


########################################
# 9. pcadapt analysis
########################################

mkdir -p pcadapt

cp \
    outflank/ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped \
    pcadapt/

cp \
    outflank/ezokuro_filtered.recode.dDocent.recode.rename2.sorted.map \
    pcadapt/

# Run:
# PCAdapt_ezokuo.R


########################################
# 10. ADMIXTURE using common outlier SNPs
########################################

cd pcadapt || exit 1

cut -d " " -f 1,4 \
    pcadapt_PC1_outflank_common_outlier.map \
    > pcadapt_PC1_outflank_common_outlier.map.list

vcftools \
    --vcf ../ezokuro_filtered.recode.dDocent.recode.rename2.vcf \
    --positions pcadapt_PC1_outflank_common_outlier.map.list \
    --out pcadapt_PC1_outflank_common_outlier.map.list \
    --recode

vcftools \
    --vcf ../ezokuro_filtered.recode.dDocent.recode.rename2.vcf \
    --positions pcadapt_PC1_outflank_common_outlier.map.list \
    --out pcadapt_PC1_outflank_common_outlier.map.list \
    --plink


# Convert to additive genotype format

plink \
    --file pcadapt_PC1_outflank_common_outlier.map.list \
    --recodeA \
    --out pcadapt_PC1_outflank_common_outlier.map.list \
    --allow-extra-chr


########################################
# 11. ADMIXTURE using outlier SNPs
########################################

mkdir -p ADMIXTURE

cp \
    pcadapt_PC1_outflank_common_outlier.map.list.ped \
    ADMIXTURE/pcadapt_PC1_outflank_common_outlier.map.list.ped

sort \
    ADMIXTURE/pcadapt_PC1_outflank_common_outlier.map.list.ped \
    > ADMIXTURE/pcadapt_PC1_outflank_common_outlier.map.list.sort.ped

cp \
    pcadapt_PC1_outflank_common_outlier.map.list.map \
    ADMIXTURE/pcadapt_PC1_outflank_common_outlier.map.list.sort.map

cd ADMIXTURE || exit 1

plink \
    --file pcadapt_PC1_outflank_common_outlier.map.list.sort \
    --recode12 \
    --out pcadapt_PC1_outflank_common_outlier.map.list.sort.ped12


for K in {1..9}
do
    admixture \
        --cv \
        pcadapt_PC1_outflank_common_outlier.map.list.sort.ped12.ped \
        "$K" \
        | tee "Admixture_log${K}.out"
done

grep -h CV ./*log*.out | sed 's/:/\t/g' > CV_error

cd ../..


########################################
# 12. IBD and IBE analyses
########################################

# Run:
# pairtial_mantel_final.R


########################################
# 13. SNP calling using all SNPs
########################################

# Run bcftools_mpileup_all_snp.sh before this step.

vcftools \
    --gzvcf ezokuro_all_snp.vcf.gz \
    --max-alleles 2 \
    --minDP 10 \
    --min-meanDP 10 \
    --minQ 20 \
    --remove-indels \
    --out ezokuro_all_filtered \
    --recode


########################################
# 14. Population genomic statistics
########################################

# NOTE:
# The original workflow uses:
# ezokuro_all_filtered2.recode.vcf
#
# Confirm how this file is generated before running the following steps.

bgzip ezokuro_all_filtered2.recode.vcf

tabix ezokuro_all_filtered2.recode.vcf.gz

python2 VCF_processing/parseVCF.py \
    -i ezokuro_all_filtered2.recode.vcf.gz \
    --skipIndels \
    | bgzip \
    > ezokuro_all_filtered2.recode.geno.gz


python2 popgenWindows.py \
    -w 20000 \
    -s 20000 \
    -m 100 \
    -g ezokuro_all_filtered2.recode.geno.gz \
    -o ezokuro_all_filtered2.recode.geno.popGenwindow.gz \
    -f phased \
    -T 5 \
    -p PO01_08Tu \
    -p PO02_KZ \
    -p PO03_AmH \
    -p PO04_KTT \
    -p PO05_KB \
    -p PO06_KD \
    -p PO07_Kgg \
    -p PO08_Fs \
    -p PO09_Ib \
    -p PO10_CB \
    -p PO11_KS \
    -p PO12_SO \
    -p PO13_L \
    -p PO14_I \
    -p PO15_KM \
    -p PO16_Sb \
    -p SJ01_RIN \
    -p SJ02_06RM \
    -p SJ03_MTN \
    -p SJ04_AmNt \
    -p SJ05_AKC \
    -p SJ06_ATS \
    -p SJ07_Tb \
    -p SJ08_YM \
    -p SJ09_AW \
    -p SJ10_Ng2 \
    -p SJ11_IF \
    -p SJ12_N \
    -p SJ13_SN \
    -p SJ14_Sg \
    -p SJ15_Ns \
    --popsFile popmap2.txt


########################################
# 15. Gene enrichment of outlier SNPs
########################################

# Test whether putative outlier SNP loci are significantly
# enriched within annotated gene regions.

perl outlier_SNP_overlap_gene_permu.pl


########################################
# 16. fastsimcoal2: preparation of SFS
########################################

python easySFS.py \
    -i ezokuro_all_filtered2.recode.vcf \
    -a \
    -p popmap_easySFS \
    -f \
    -o ezo_kuroSFS \
    --prefix ezo_kuroSFS \
    --preview


python easySFS.py \
    -i ezokuro_all_filtered2.recode.vcf \
    -a \
    -p popmap_easySFS \
    -f \
    -o ezo_kuroSFS \
    --prefix ezo_kuroSFS \
    --proj=40,40



