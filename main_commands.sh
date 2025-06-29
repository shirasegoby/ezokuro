
#!/bin/sh

# SNP calling
BAMS=""
for bam in `ls ./ezokuro_realigned_bam/*.bam.realigned.bam.uniq.bam`
do
  BAMS=${BAMS}" ${bam} "
done

bcftools mpileup -Ou --annotate FORMAT/DP -f $REF $BAMS | bcftools call --threads 10 -mv -Oz -o ezokuro.vcf.gz


#filtering SNP
vcftools --gzvcf ezokuro.vcf.gz  --max-alleles 2  --minDP 10 --min-meanDP 10 --minQ 20 --remove-indels --out  ezokuro_filtered --recode

sh './dDocent/scripts/pop_missing_filter.sh' ezokuro_filtered.recode.vcf ./popmap.txt 0.2 2 ezokuro_filtered.recode.dDocent


#convert to genepop
java -Xmx2048m -Xms1024m -jar /Users/shotarohirase/Desktop/Genomic_tools/PGDSpider_2.1.1.5/PGDSpider2-cli.jar \
-inputfile ezokuro_filtered.recode.dDocent.recode.rename2.vcf -outputfile ezokuro_filtered.recode.dDocent.recode.rename.gen \
-outputformat GENEPOP -spid vcf2genepop.spid

#sort genepop
grep -v -e "SNP" -e "Pop" ezokuro_filtered.recode.dDocent.recode.rename.gen |sort -  > ezokuro_filtered.recode.dDocent.recode.rename2.gen 

grep -e "SNP" -e "Pop" ezokuro_filtered.recode.dDocent.recode.rename.gen > ezokuro_filtered.recode.dDocent.recode.rename2.gen.header

cat ezokuro_filtered.recode.dDocent.recode.rename2.gen.header ezokuro_filtered.recode.dDocent.recode.rename2.gen  > ezokuro_filtered.recode.dDocent.recode.rename3.gen

rm ezokuro_filtered.recode.dDocent.recode.rename2.gen ezokuro_filtered.recode.dDocent.recode.rename2.gen.header

gsed -i '1s/^/#all_pop\n/' ezokuro_filtered.recode.dDocent.recode.rename3.gen

# mofify genepop file by manual

# make PO and SJ genepop files
grep -e "#" -e "SNP" -e "Pop" -e "PO" ezokuro_filtered.recode.dDocent.recode.rename3.gen > ezokuro_filtered.recode.dDocent.recode.rename.PO.gen
grep -e "#" -e "SNP" -e "Pop" -e "SJ" ezokuro_filtered.recode.dDocent.recode.rename3.gen > ezokuro_filtered.recode.dDocent.recode.rename.SJ.gen


# convert to ped file
vcftools --vcf ezokuro_filtered.recode.dDocent.recode.rename2.vcf --plink  --out  ezokuro_filtered.recode.dDocent.recode.rename2

sort ezokuro_filtered.recode.dDocent.recode.rename2.ped > ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped

cp ezokuro_filtered.recode.dDocent.recode.rename2.map ezokuro_filtered.recode.dDocent.recode.rename2.sorted.map

/Users/shotarohirase/Desktop/Genomic_tools/plink_mac/plink  --file ezokuro_filtered.recode.dDocent.recode.rename2.sorted --recode12 --out ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped12


# ADMIXTURE without filtering

mkdir ./ADMIXTURE2_remove_Ng2_04_nofilter

cp ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped12.ped ./ADMIXTURE2_remove_Ng2_04_nofilter
cd ./ADMIXTURE2_remove_Ng2_04_nofilter

for K in 1 2 3 4 5 6 7 8 9;

do admixture --cv ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped12.ped $K | tee Admixture_log${K}.out; 

done

grep -h CV ./*log*.out |gsed "s/:/\t/g" - > CV_error

# make barplot by Pophelper_ADMIXTURE_plot2.R

# ADMIXTURE with filtering

# filtering
vcftools --vcf ezokuro_filtered.recode.dDocent.recode.rename2.vcf  --thin 1000 --plink  --out  ezokuro_filtered.recode.dDocent.recode.rename2.thin1000 

#convert to ped12
/Users/shotarohirase/Desktop/Genomic_tools/plink_mac/plink  --file ezokuro_filtered.recode.dDocent.recode.rename2.thin1000  --recode12 --out ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12

sort ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12.ped > ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12.sorted.ped

cp ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12.sorted.ped ./ADMIXTURE2_remove_Ng2_04

cd ./ADMIXTURE2_remove_Ng2_04

for K in 1 2 3 4 5 6 7 8 9;

do /Users/shotarohirase/Desktop/Genomic_tools/ADMIXTURE/admixture_macosx-1.3.0/admixture --cv ezokuro_filtered.recode.dDocent.recode.rename2.thin1000.ped12.sorted.ped $K | tee Admixture_log${K}.out; 

done

grep -h CV ./*log*.out |gsed "s/:/\t/g" - > CV_error

# make barplot by Pophelper_ADMIXTURE_plot2.Rでグラフをつくる。


# calculate genetic diversity with genepop file and genodive
cd ../

cp ezokuro_filtered.recode.dDocent.recode.rename3.gen ./Genetic_diversity/

# genome scan with outflank

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

#convert to codeA
/Users/shotarohirase/Desktop/Genomic_tools/plink_mac/plink --file ezokuro_filtered.recode.dDocent.recode.rename2.sorted   --recodeA --out  ezokuro_filtered.recode.dDocent.recode.rename2.sorted.OUTflank2 --allow-extra-chr

#running outflank_ezo_kuro.R

cd ../


# genome scan with pcadapt

cp ./outflank/ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped ./pcadapt/
cp ./outflank/ezokuro_filtered.recode.dDocent.recode.rename2.sorted.map ./pcadapt/

# running PCAdapt_ezokuo.R


#common outlier

cd ./common_outlier

#outflankとPCAdaptのcommon outlierを抽出する。
#southern_allele_freq*.plのSNP数（i）を変更する。
#common_outlier.Rで、common outlierと頻度分布の図をつくる。
#common outlierのpedファイルをつくる。0,1,2になっているやつ。

cut -d " " -f 1,4 pcadapt_PC1_outflank_common_outlier.map > pcadapt_PC1_outflank_common_outlier.map.list

vcftools --vcf ../ezokuro_filtered.recode.dDocent.recode.rename2.vcf --positions /Users/shotarohirase/Desktop/エゾクロ解析re/common_outlier/pcadapt_PC1_outflank_common_outlier.map.list --out /Users/shotarohirase/Desktop/エゾクロ解析re/common_outlier/pcadapt_PC1_outflank_common_outlier.map.list --recode
vcftools --vcf ../ezokuro_filtered.recode.dDocent.recode.rename2.vcf --positions /Users/shotarohirase/Desktop/エゾクロ解析re/common_outlier/pcadapt_PC1_outflank_common_outlier.map.list --out /Users/shotarohirase/Desktop/エゾクロ解析re/common_outlier/pcadapt_PC1_outflank_common_outlier.map.list --plink

#convert to codeA
/Users/shotarohirase/Desktop/Genomic_tools/plink_mac/plink --file pcadapt_PC1_outflank_common_outlier.map.list   --recodeA --out  pcadapt_PC1_outflank_common_outlier.map.list --allow-extra-chr


#common outlierでADMIXTUREをする
mkdir ADMIXTURE
cd ADMIXTURE

cp ../pcadapt_PC1_outflank_common_outlier.map.list.ped ./pcadapt_PC1_outflank_common_outlier.map.list.ped

sort ./pcadapt_PC1_outflank_common_outlier.map.list.ped > ./pcadapt_PC1_outflank_common_outlier.map.list.sort.ped

cp ../pcadapt_PC1_outflank_common_outlier.map.list.map  ./pcadapt_PC1_outflank_common_outlier.map.list.sort.map 

/Users/shotarohirase/Desktop/Genomic_tools/plink_mac/plink  --file ./pcadapt_PC1_outflank_common_outlier.map.list.sort  --recode12 --out ./pcadapt_PC1_outflank_common_outlier.map.list.sort.ped12

for K in 1 2 3 4 5 6 7 8 9;

do /Users/shotarohirase/Desktop/Genomic_tools/ADMIXTURE/admixture_macosx-1.3.0/admixture --cv ./pcadapt_PC1_outflank_common_outlier.map.list.sort.ped12.ped $K | tee Admixture_log${K}.out; 

done

grep -h CV ./*log*.out |gsed "s/:/\t/g" - > CV_error

#IDの付いた2.Qを作る。
cut -f 1 pcadapt_PC1_outflank_common_outlier.map.list.sort.ped |paste - pcadapt_PC1_outflank_common_outlier.map.list.sort.ped12.2.Q.plus_indID 



################################################################################################
# POとSJのそれぞれのgenepopファイルを使って、IBDとIBDを行う。

#pairtial_mantel.Rで実行する



##########################################################
# piを算出するために、すべてのsnpをコールする

bcftools_mpileup20220425_all_snp_ver.sh

vcftools --gzvcf /Users/shotarohirase/Desktop/エゾクロ解析re/ezokuro_all_snp.vcf.gz --max-alleles 2  --minDP 10 --min-meanDP 10 --minQ 20 --remove-indels    --out ezokuro_all_filtered --recode


#dDocentのスクリプトでfilteringする。
#時間かかりすぎるので、やらない。
#sh '/Users/shotarohirase/Desktop/Genomic_tools/dDocent/scripts/pop_missing_filter.sh' ezokuro_all_filtered.recode.vcf ./popmap.txt 0.2 2 ezokuro_all_filtered.recode.dDocent

#個体名をrenameする
cp ezokuro_all_filtered.recode.vcf ./rename/
sh rename.sh
cp ./rename/ezokuro_all_filtered.recode.vcf ./

gsed -e "s/YM08/SJ08_YM08/g" ezokuro_all_filtered.recode.vcf > ezokuro_all_filtered2.recode.vcf

bgzip ezokuro_all_filtered2.recode.vcf

tabix ezokuro_all_filtered2.recode.vcf.gz

python2 /Users/shotarohirase/Desktop/Genomic_tools/genomics_general_Martin_Simon/VCF_processing/parseVCF.py  -i ezokuro_all_filtered2.recode.vcf.gz --skipIndels | bgzip > ezokuro_all_filtered2.recode.geno.gz

python2 /Users/shotarohirase/Desktop/Genomic_tools/genomics_general_Martin_Simon/popgenWindows.py  -w 20000 -s 20000 -m 100 -g   ./ezokuro_all_filtered2.recode.geno.gz  -o  ezokuro_all_filtered2.recode.geno.popGenwindow.gz -f phased -T 5 -p PO01_08Tu  -p PO02_KZ  -p PO03_AmH  -p PO04_KTT  -p PO05_KB  -p PO06_KD  -p PO07_Kgg  -p PO08_Fs  -p PO09_Ib  -p PO10_CB  -p PO11_KS  -p PO12_SO  -p PO13_L  -p PO14_I  -p PO15_KM  -p PO16_Sb  -p SJ01_RIN  -p SJ02_06RM  -p SJ03_MTN  -p SJ04_AmNt  -p SJ05_AKC  -p SJ06_ATS  -p SJ07_Tb  -p SJ08_YM  -p SJ09_AW  -p SJ10_Ng2  -p SJ11_IF  -p SJ12_N  -p SJ13_SN  -p SJ14_Sg  -p SJ15_Ns --popsFile ./popmap2.txt



# Red abalone研究との関連

grep -i -w -f <(echo -e "contig147173\ncontig148731\ncontig83864\ncontig89165\ncontig24507\ncontig90217\ncontig104796\ncontig46443\ncontig83915\ncontig88626\ncontig14237\ncontig26566\ncontig115444\ncontig83774\ncontig83577\ncontig25921\ncontig90598\ncontig102343\ncontig88857\nContig88612\ncontig83945\ncontig84742\ncontig24566\ncontig82770\ncontig88923\ncontig24491\ncontig147114\ncontig27448\ncontig16119\ncontig13864\ncontig137187\ncontig137221\ncontig83712\ncontig73178\ncontig102168\ncontig24945\ncontig88834")  /Users/shotarohirase/Desktop/エゾクロ解析re/outlier_gene_overlap_analysis/Red_abalone_RNAseq_supplement_doi_10_5061_dryad_85p80__v20121031/supplementary_materials/suppl_3.txt   > /Users/shotarohirase/Desktop/エゾクロ解析re/outlier_gene_overlap_analysis/Red_abalone_RNAseq_supplement_doi_10_5061_dryad_85p80__v20121031/supplementary_materials/suppl_3_candidate_contig2.txt 

/Users/shotarohirase/Desktop/Genomic_tools/ncbi-blast-2.15.0+/bin/makeblastdb -in /Users/shotarohirase/Desktop/エゾクロ解析re/outlier_gene_overlap_analysis/Red_abalone_RNAseq_supplement_doi_10_5061_dryad_85p80__v20121031/supplementary_materials/suppl_3_candidate_contig.fa  -dbtype nucl
/Users/shotarohirase/Desktop/Genomic_tools/ncbi-blast-2.15.0+/bin/tblastn -db /Users/shotarohirase/Desktop/エゾクロ解析re/outlier_gene_overlap_analysis/Red_abalone_RNAseq_supplement_doi_10_5061_dryad_85p80__v20121031/supplementary_materials/suppl_3_candidate_contig.fa -query /Users/shotarohirase/Desktop/エゾクロ解析re/outlier_gene_overlap_analysis/candidate_gene_uniq.bed.pep -evalue 0.0001 -outfmt 6 > /Users/shotarohirase/Desktop/エゾクロ解析re/outlier_gene_overlap_analysis/candidate_gene_uniq.bed.pep.red_abalone_contig_tblasn0.0001




























