use strict;

#fasta::orf_genomic_all.fasta:YAL001C -              YAL001C        -              N    3573   37.08                      TFC3 SGDID:S000000001, Chr I from 151166-147594, Genome Release 64-2-1, reverse complement, Verified ORF, "Subunit of RNA polymerase III transcription initiation factor complex; part of the TauB domain of TFIIIC that binds DNA at the BoxB promoter sites of tRNA and similar g


my %lens;
my %cgs;

open IS, $ARGV[0]; #infoseqfile 


while(<IS>){
    my @line=split(/\s+/, $_);
    $lens{$line[2]}=$line[5];
    $cgs{$line[2]}=$line[6];
    #print "|$line[2]|$lens{$line[2]}|$cgs{$line[2]}|\n";

}

close IS;
#YAR010C_NumOfGenes_27	YJL033W	57.14	168	59	13	757	916	1813	1975	1.9	32.5
open MVIII, $ARGV[1];

while(<MVIII>){
        chomp $_;
        my @line=split(/\s+/, $_);
        print "$_\t$lens{$line[0]}\t$lens{$line[1]}\t$cgs{$line[0]}\t$cgs{$line[1]}\n";
        
}

close MVIII;