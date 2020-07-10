#!/usr/bin/perl

use strict;
use lib '/Users/mdechiara/Downloads/BioPerl-1.6.924';

use Bio::Seq;
use Bio::SeqIO;

my %accepted;
open IO, "$ARGV[0]"; #list of orf to be taken 
while (<IO>) {
    chomp $_;
    my @data=split(/\t/, $_);
    $accepted{$data[1]}{'len'}=$data[2];
    $accepted{$data[1]}{'clus'}=$data[0];
    
}

my $seqi = Bio::SeqIO->new('-format' =>  "fasta",
			    -file => "$ARGV[1]"
				      );


unless (-e "$ARGV[1].collapsed.fromgraph.fasta"){
system "touch $ARGV[1].collapsed.fromgraph.fasta";
}


while(my $seq = $seqi->next_seq ) {
    my $len=$seq->length();
    my $id = $seq->id();
    my $nuc = $seq->seq();
    chomp $id;
    my $newid=$id;
    if ($accepted{$id}{'len'}>1) {
        $newid="$id"."_NumOfGenes_"."$accepted{$id}{'len'}";
    }
    if($accepted{$id}{'len'}>0) {
        if ($accepted{$id}{'len'}>1) {
            $newid="$id"."_NumOfGenes_"."$accepted{$id}{'len'}";
        }
        my $seq_obj = Bio::Seq->new(
			      -seq => $nuc,
			      -display_id  => "$newid" );
	my $io = Bio::SeqIO->new(-format => "fasta", -file => ">>$ARGV[1].collapsed.fromgraph.fasta" );  
        $io->write_seq($seq_obj);
    }

    
}
