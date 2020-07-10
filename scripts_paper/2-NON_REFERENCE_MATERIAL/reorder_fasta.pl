#!/usr/bin/perl -w
use strict;
use POSIX;
#use lib '/Users/mdechiara/Downloads/BioPerl-1.6.924';
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;

my $fasta_library = $ARGV[1]; #assembly
my $database = Bio::DB::Fasta->new("$fasta_library") or die "Failed to creat Fasta DP object on fasta library\n";

open IO, "$ARGV[0]"or die "NOT OPENED THE INFOSEQ";# output of "order_infoseq on infoseq name.fasta"
my %elements;
while (<IO>) {
    chomp $_;
    my @row=split(/\s+/, $_);
    if ($row[1]=~m/\d+/){
        my $query=$row[0];
        #print "|$query|\n";
        my $getseq = $database->get_Seq_by_id($query);
        my $sequio;
        my $seq_obj = Bio::Seq->new(
			      -seq => $getseq->seq($query),
			  #    -length      => $stop-$Locuses{$stop}{'start'}+1,
			      -display_id  => $query );
        my $io = Bio::SeqIO->new(-format => "fasta", -file => ">>$ARGV[1].ordered.fasta" );  
        $io->write_seq($seq_obj);
    }
    
}
