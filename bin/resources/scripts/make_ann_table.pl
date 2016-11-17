#!/usr/bin/perl -w
# perl make_ann_table.pl rno.gff3 rno_mature.fa rno_hairpin.fa rno_table.txt
# works on miRbase r19 files. Generate annotation tables from gff file and mature and hairpin fasta files.
# J.Baran
my %gffHash=();

$gff=shift;

open (IN, $gff)or die;
while ($line = <IN>) {
   chomp($line);
   unless ($line =~/\#/) {
      my ($chr,$dum,$type,$start,$stop,$d2,$strand,$d3,$info) = split(/\t/,$line);
      unless ($type eq "miRNA") {

	 my($acc,$id,$name) = split(/;/,$info);
	 my($code,$mir) = split(/=/,$name);
	 $mir = lc($mir);
	 @parts = split(/-/,$mir);
	 $base=$parts[0]."-".$parts[1]."-".$parts[2];


	 $field = join(":",$chr,$start,$stop,$strand,$mir);
	 push(@{$gffHash{$base}} , $field );
      }
   }
}
close(IN);



%mirseqs =();
# mouseMirnas
$mirs = shift;
open (IN, $mirs);
while ($line = <IN>) {
   chomp($line);
   if ($line =~ />/) {
      @arr = split(" ", $line);
      $mir = substr($arr[0], 1);
      $mir = lc($mir);
      $id = $arr[1];
      @parts = split(/-/,$mir);
      $base=$parts[0]."-".$parts[1]."-".$parts[2];

   }
   $matureseq = <IN>;
   chomp($matureseq);
   $matureseq =~ s/U/T/g;
   $end="";
   if ($parts[-1] =~/p/) {
      $end = $parts[-1];
      pop(@parts);
   }
   $name = join("-",@parts);
   $nparts = @parts;
   foreach $entry (@{$gffHash{$base}}) {
      @e_parts = split(/:/,$entry);
      if(($nparts ==3) or ($name eq $e_parts[-1])){
	 $kname = $e_parts[-1];
	 $kname = join("-",$e_parts[-1] , $end) unless ($end eq "");
	 $mirseqs{$kname} {M} = $matureseq;
	 $mirseqs{$kname} {C} = $e_parts[0];
	 $mirseqs{$kname} {S} = $e_parts[1];
	 $mirseqs{$kname} {E} = $e_parts[2];
	 $mirseqs{$kname} {D} = $e_parts[3];
      }
   }
}
close (IN);

$hps=shift;
open (IN, $hps)or die;
$line = <IN>;
while (defined ($line)) {
   chomp($line);
   if ($line =~ />/) {
      @arr = split(" ", $line);
      $mir = substr($arr[0], 1);
      $mir=lc($mir);
      $id = $arr[1];
      @parts = split(/-/,$mir);
      $base=$parts[0]."-".$parts[1]."-".$parts[2];
   }
   $line = <IN>; chomp($line);
   $hp = "";
   until (!defined($line) ||($line =~/>/)) {
      $hp= $hp.$line;
      $line = <IN>;
      chomp($line) if defined($line);
   }
   $hp =~ s/U/T/g;

   $name = join("-",@parts);
   $nparts = @parts;
   foreach $entry (@{$gffHash{$base}}) {
      @e_parts = split(/:/,$entry);
      if(($nparts ==3) or ($name eq $e_parts[-1])){
	 $kname = join("-",$e_parts[-1] , "5p");
	 if (exists ($mirseqs{$kname})){
	    $mirseqs{$kname}{H} = $hp;
	 }
	 $kname = join("-",$e_parts[-1] , "3p");
	 if (exists ($mirseqs{$kname})){
	    $mirseqs{$kname}{H} = $hp;
	 }
	 $kname = $e_parts[-1];
	 if (exists ($mirseqs{$kname})){
	    $mirseqs{$kname}{H} = $hp;
	 }
      }
   }
}
close (IN);


$file =shift;
open (OUT, ">$file");

foreach $k(keys(%mirseqs)){
  print OUT "$k\t$mirseqs{$k}{C}\t$mirseqs{$k}{S}\t$mirseqs{$k}{E}\t$mirseqs{$k}{D}\t$mirseqs{$k}{M}\t$mirseqs{$k}{H}\t\n";
}


