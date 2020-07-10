#!/usr/bin/perl -w
use strict;


#fasta::REFPtemp.kept.afterN.fasta:ARV_6-273-R.3 -              ARV_6-273-R.3  -              N    664    24.55
use strict;

open IO, "$ARGV[0]";# output of "infoseq name.fasta"
my %elements;
while (<IO>) {
    chomp $_;
    my @row=split(/\s+/, $_);
    my $seq=$row[2];
    my $len=$row[5];
    my $cg=$row[6];
   # print "$seq\t$len\t$cg\n";
   if ($len=~m/\d+/){
    $elements{$len}{$cg}{$seq}=1;
   }
}

print"Name\tLength\t%GC\n";
foreach my $size (sort {$b <=> $a} keys %elements){
    foreach my $tri (sort {$a <=> $b} keys %{$elements{$size}}){
        foreach my $ctg (sort {$a cmp $b} keys %{$elements{$size}{$tri}}){
            print "$ctg\t$size\t$tri\n";
        }
    }
    
}