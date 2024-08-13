#!/usr/bin/perl -w


open(OUT, "> permu_overlap_number_annotated_bed.txt") or die("error :$!");


system ("intersectBed -a pcadapt_PC1_outflank_common_outlier.map.list.map.bed -b ./transcripts.fasta.transdecoder.genome.gtf.bed > temp2.bed");
system ("sort temp2.bed |uniq - > temp3.bed");
open(FILE, "temp3.bed") or die "Can't open temp2.bed: $!";
@array  = <FILE>; 
$lines = @array;

print OUT "oberved\t$lines\n";

system ("rm temp3.bed");
system ("rm temp2.bed");
system ("rm temp.bed");


for ($i=1; $i<=1000; $i++){

print "$i\t";

my @array="";

#ランダムに448SNPを抽出し、geneとオーバーラップする数を調べる。
system ("shuf -n 448 ./ezokuro_filtered.recode.dDocent.recode.rename2.map.bed  > temp.bed");

system ("intersectBed -a temp.bed -b ./transcripts.fasta.transdecoder.genome.gtf.bed > temp2.bed");

system ("sort temp2.bed |uniq - > temp3.bed");


open(FILE, "temp3.bed") or die "Can't open temp2.bed: $!";
@array  = <FILE>; 
$lines = @array;

print OUT "$lines\n";

system ("rm temp3.bed");
system ("rm temp2.bed");
system ("rm temp.bed");

}


system ("sort permu_overlap_number_annotated_bed.txt > permu_overlap_number_annotated_bed_sorted.txt");