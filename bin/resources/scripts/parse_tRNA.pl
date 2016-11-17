#!/usr/bin/perl -w

# Generate tRNA annotation files
# New genome (example rn4):
#  tRNA files need to be generated from the genomic tRNA database (http://gtrnadb.ucsc.edu/)
#  bedfile: mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, chromStart, chromEnd, name, score, strand from rn4.tRNAs"  > rn4_tRNA.bed
#  bed12: (http://gtrnadb.ucsc.edu/Rnorv/Rnorv-by-locus-txt.html) copy from text version of the tRNAs by locus page of the genomic tRNA database (rn4_tRNA.txt - tab delim)
#  then run scripts/parse_tRNA.pl rn4_tRNA.bed rn4_tRNA.txt > rn4_tRNA12.bed
#  mature tRNA file:
#  fastaFromBed -s -split -fi ../genome/rn4.fa -bed rn4_tRNA12.bed -s rn4_mature_tRNA_LIB.fa
#  J.Baran

$bedfile=shift;
$exon=shift;

open (EX,$exon) or die;
while(<EX>) {
   chomp();
   my($ch,$num,$st,$end,$AA,$AC,$leader,$trailer,$exst,$exed,$s1,$s2,$s3,$n) =split(/\t/,$_);
   $key = "$ch.tRNA$num-$AA$AC";
   $HexonSt{$key}=$exst;
   $HexonEd{$key}=$exed;
}
close (EX);
open (BED, $bedfile) or die;
while(<BED>) {
   chomp();
   my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$_);
   if ($strand eq "+") {
   $es=$HexonSt{$name} ;
   $ee=$HexonEd{$name} ;
   }
   else {
   $ee=$HexonSt{$name} ;
   $es=$HexonEd{$name} ;
   }
   $offset = $es -$start;
   $len = $ee-$es;
   if ($es==0) {
      $N=1;
      $blocksz=$end-$start;
      $blockst=$start;
   }
   else {
      $N =2;
      $A=$es-$start-1;
      $B=$end-$ee;
      $blocksz=$A .",".$B; 
      $C=$ee-$start;
      $blockst="0,". $C;
   }
# print ">>$name\t$start\t$es\t$ee\t$offset\t$len\n";
   print "$_";
   print "\t$start\t$end\t0,0,0\t$N\t$blocksz\t$blockst\n";

}
close(BED);
