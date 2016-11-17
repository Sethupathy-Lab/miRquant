#!/usr/bin/perl -w
# perl gen_table_mirLoc.pl mmu_table.txt > mmu_tableL.bed
# Generates coordinates for mature miRNAs from precursor coordiante table
# J. Baran

%mirseqs =();
$mirs = shift;
open (IN, $mirs);
while ($line = <IN>) {
   chomp($line);
   my ($mName,$mchr,$mSt,$mEd,$mStr,$mSeq,$mHp) = split(/\t/,$line);
   $ind = index($mHp, $mSeq);                                                                                                                                                        
   if ($mStr eq "+") {
      $locA = $mSt + $ind;
      $locB = $mSt + $ind + length($mSeq) -1;
   }
   else {
      $locB = $mEd - $ind  if ($mStr eq "-");
      $locA = $mEd - $ind -(length($mSeq)-1);
   }
   $mName2="$mName:$mSeq";
   $out = join ("\t",$mchr,$locA,$locB,$mName2,1,$mStr);

   print "$out\n";
   
   $baseName = $mName;
   $baseName =~ s/-\dp//ge; 
   push (@{$mirseqs{$baseName}},$mSeq);
   $mirHp{$baseName} = $mHp;

}

for $key (keys(%mirHp)) {
   if (scalar(@{$mirseqs{$key}})==2) {
      $hP = $mirHp{$key};
      $seq1 = ${$mirseqs{$key}}[0];
      $seq2 = ${$mirseqs{$key}}[1];

      $ind1 = index($hP,$seq1);
      $ind2 = index($hP,$seq2);

      if ($ind1 < $ind2) {
	 $sep = $ind2 - ($ind1+length($seq1) -1);
	 $l = length($seq1) - 1;
	 print "$key\t$sep\t$ind2\t$ind1\t$l\n";
      }
      else {
	 $sep = $ind1 - ($ind2+length($seq2) -1);
	 $l = length($seq2) - 1;
	 print "$key\t$sep\t$ind1\t$ind2\t$l\n";
      }


   }
}

