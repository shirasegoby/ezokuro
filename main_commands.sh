#!/bin/sh

#SNP calling

#Use bcftools_mpileup.sh

#SNP filtering
vcftools --gzvcf ezokuro.vcf.gz --max-alleles 2  --minDP 10 --min-meanDP 10 --minQ 20 --remove-indels    --out  ezokuro_filtered --recode

#dDocent filtering
sh './dDocent/scripts/pop_missing_filter.sh' ezokuro_filtered.recode.vcf ./popmap.txt 0.2 2 ezokuro_filtered.recode.dDocent

# convert vcf into ped
vcftools --vcf ezokuro_filtered.recode.dDocent.recode.rename2.vcf --plink  --out  ezokuro_filtered.recode.dDocent.recode.rename2
sort ezokuro_filtered.recode.dDocent.recode.rename2.ped > ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped

# SNP filtering for ADMIXTURE
vcftools --vcf ezokuro_filtered.recode.dDocent.recode.rename2.vcf  --thin 1000 --plink  --out  ezokuro_filtered.recode.dDocent.recode.rename2.thin1000 

# convert into ped12
./plink_mac/plink  --file ezokuro_filtered.recode.dDocent.recode.rename2.thin1000  --recode12 --out ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12
sort ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12.ped > ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12.sorted.ped

# run ADMIXTURE
for K in 1 2 3 4 5 6 7 8 9;
do ./admixture_macosx-1.3.0/admixture --cv ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12.sorted.ped $K | tee Admixture_log${K}.out; 
done
grep -h CV ./*log*.out |gsed "s/:/\t/g" - > CV_error

# Make ADMIXTURE plot
# run Pophelper_ADMIXTURE_plot2.R

# Outflank analysis
cp ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped ./outflank/
cp ezokuro_filtered.recode.dDocent.recode.rename2.map ./outflank/
cd ./outflank/

# modify map file
gsed -i -e "s/:/\t/g" ezokuro_filtered.recode.dDocent.recode.rename2.map
SNP_number=`wc -l ezokuro_filtered.recode.dDocent.recode.rename2.map | awk '{print $1}'`
seq -s "\nSNP_" 0 $SNP_number > tmp.txt
SNP_number2=$((SNP_number+1))
sed '1d' tmp.txt > tmp2.txt 
sed "${SNP_number2}d" tmp2.txt  >tmp3.txt
cut -f 4,5 ezokuro_filtered.recode.dDocent.recode.rename2.map > tmp4.txt
cut -f 2 ezokuro_filtered.recode.dDocent.recode.rename2.map |paste - tmp3.txt tmp4.txt |cat - | gsed -e "s/\t/ /g" - > ezokuro_filtered.recode.dDocent.recode.rename2.sorted.map
rm tmp*.txt


#Outflank
#codeA file make
./plink_mac/plink --file ezokuro_filtered.recode.dDocent.recode.rename2.sorted   --recodeA --out  ezokuro_filtered.recode.dDocent.recode.rename2.sorted.OUTflank2 --allow-extra-chr

# run outflank_ezo_kuro.R

# pcadapt analysis
cp ./outflank/ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped ./pcadapt/
cp ./outflank/ezokuro_filtered.recode.dDocent.recode.rename2.sorted.map ./pcadapt/

# run PCAdapt_ezokuo.R

#ADMIXTURE using common outlier SNPs
cut -d " " -f 1,4 pcadapt_PC1_outflank_common_outlier.map > pcadapt_PC1_outflank_common_outlier.map.list
vcftools --vcf ../ezokuro_filtered.recode.dDocent.recode.rename2.vcf --positions /Users/shotarohirase/Desktop/エゾクロ解析re/common_outlier/pcadapt_PC1_outflank_common_outlier.map.list --out /Users/shotarohirase/Desktop/エゾクロ解析re/common_outlier/pcadapt_PC1_outflank_common_outlier.map.list --recode
vcftools --vcf ../ezokuro_filtered.recode.dDocent.recode.rename2.vcf --positions /Users/shotarohirase/Desktop/エゾクロ解析re/common_outlier/pcadapt_PC1_outflank_common_outlier.map.list --out /Users/shotarohirase/Desktop/エゾクロ解析re/common_outlier/pcadapt_PC1_outflank_common_outlier.map.list --plink

#convert to codeA
/Users/shotarohirase/Desktop/Genomic_tools/plink_mac/plink --file pcadapt_PC1_outflank_common_outlier.map.list   --recodeA --out  pcadapt_PC1_outflank_common_outlier.map.list --allow-extra-chr

mkdir ADMIXTURE
cd ADMIXTURE
cp ../pcadapt_PC1_outflank_common_outlier.map.list.ped ./pcadapt_PC1_outflank_common_outlier.map.list.ped
sort ./pcadapt_PC1_outflank_common_outlier.map.list.ped > ./pcadapt_PC1_outflank_common_outlier.map.list.sort.ped
cp ../pcadapt_PC1_outflank_common_outlier.map.list.map  ./pcadapt_PC1_outflank_common_outlier.map.list.sort.map 
./plink_mac/plink  --file ./pcadapt_PC1_outflank_common_outlier.map.list.sort  --recode12 --out ./pcadapt_PC1_outflank_common_outlier.map.list.sort.ped12

# run admixture
for K in 1 2 3 4 5 6 7 8 9;
do ./admixture_macosx-1.3.0/admixture --cv ./pcadapt_PC1_outflank_common_outlier.map.list.sort.ped12.ped $K | tee Admixture_log${K}.out; 
done
grep -h CV ./*log*.out |gsed "s/:/\t/g" - > CV_error

# IBD and IBE analyses
# use pairtial_mantel_final.R

# SNP calling (all SNPs)
#use bcftools_mpileup_all_snp.sh

vcftools --gzvcf ./ezokuro_all_snp.vcf.gz --max-alleles 2  --minDP 10 --min-meanDP 10 --minQ 20 --remove-indels    --out ezokuro_all_filtered --recode

bgzip ezokuro_all_filtered2.recode.vcf
tabix ezokuro_all_filtered2.recode.vcf.gz
python2 ./VCF_processing/parseVCF.py  -i ezokuro_all_filtered2.recode.vcf.gz --skipIndels | bgzip > ezokuro_all_filtered2.recode.geno.gz
python2 popgenWindows.py  -w 20000 -s 20000 -m 100 -g   ./ezokuro_all_filtered2.recode.geno.gz  -o  ezokuro_all_filtered2.recode.geno.popGenwindow.gz -f phased -T 5 -p PO01_08Tu  -p PO02_KZ  -p PO03_AmH  -p PO04_KTT  -p PO05_KB  -p PO06_KD  -p PO07_Kgg  -p PO08_Fs  -p PO09_Ib  -p PO10_CB  -p PO11_KS  -p PO12_SO  -p PO13_L  -p PO14_I  -p PO15_KM  -p PO16_Sb  -p SJ01_RIN  -p SJ02_06RM  -p SJ03_MTN  -p SJ04_AmNt  -p SJ05_AKC  -p SJ06_ATS  -p SJ07_Tb  -p SJ08_YM  -p SJ09_AW  -p SJ10_Ng2  -p SJ11_IF  -p SJ12_N  -p SJ13_SN  -p SJ14_Sg  -p SJ15_Ns --popsFile ./popmap2.txt

#whether these 448 putative SNP loci were significantly enriched within annotated gene regions
perl outlier_SNP_overlap_gene_permu.pl

#fastsimcoal2
python ./easySFS.py -i ezokuro_all_filtered2.recode.vcf -a -p popmap_easySFS -f -o ezo_kuroSFS --prefix ezo_kuroSFS --preview
python ./easySFS.py -i ezokuro_all_filtered2.recode.vcf -a -p popmap_easySFS -f -o ezo_kuroSFS --prefix ezo_kuroSFS --proj=40,40

# Run fastsimcoal2 using script in fastsimcoal2 folder




