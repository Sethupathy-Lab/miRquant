#! /nas02/home/r/o/roachjm/bin/perl -w

# mirBaseFa2mmuOther.pl mature.fa rno
# mirBaseFa2mmuOther.pl hairpin.fa rno
# Pulls out annotations for a specific species (rno) from input fasta files.
# J.Baran


my $filename = $ARGV[0];
my $species = $ARGV[1];

if ($species eq 'mmu') {
   $spName='Mus musculus';
}
elsif ($species eq 'rno') {
   $spName='Rattus norvegicus';
}
else { # hsa
   $spName='Homo sapiens';
}
my @parts= split(/\//,$filename);
print "@parts";
print "\n$parts[$#parts]\n";
my $name = $parts[$#parts];
$parts[$#parts] = join("_",$species , $name);
my $mmufile= join('/', @parts);
$parts[$#parts] = join("_","oth" ,$name);
my $othfile= join('/', @parts);
print "SPECIES: $mmufile\n";
# print "OTH: $othfile\n";
open(FA, $filename)
   or die "Couldn't open $filename for reading!";
open(SPECIES,">$mmufile");
# open(OTH,">$othfile");
my $type =0;
while (<FA>) {
   chomp;
   my $currentLine = $_;

   if (/>/) { #Header Line
      $index = index ($currentLine,$spName);
      if ($index != -1) {
	 $type = 1;
      }
      else{
	 $type = 0;
      }

   }
   print SPECIES "$currentLine\n", if ($type==1);
# print OTH "$currentLine\n", if ($type==0);

        

}
close(FA);
close(SPECIES);
# close(OTH);

#######################################################
