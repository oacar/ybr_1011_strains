use strict;
use warnings;
#use lib '/Users/mdechiara/Downloads/BioPerl-1.6.924';

use Bio::SeqIO;
my $coordfile=$ARGV[0];
my $refseq=$ARGV[1];#my $queryseq=$ARGV[2];
my $out="$refseq.kept.fasta";

#$out=~s/.fa/.kept.fasta/;

my %refhash;

my $seqio = Bio::SeqIO->new('-format' =>  "fasta",
			    -file => "$refseq"
				      );

while(my $seq = $seqio->next_seq ) {
     
    my $len=$seq->length();
    my $id = $seq->id();#    print "$coordfile\t$id\t$len\n";

    $refhash{$id}=$len;
}


my %CTGrefgen;
my $c=0;
open COORDFILE, "$coordfile";

while (<COORDFILE>){
        chomp $_;
        my @line=split(/\t/, $_);
        my $start=$line[1];#	print "$start\t";
        my $end=$line[2];
        my $ctg=$line[0];
	$ctg =~s/-*//;
	my $len=$end-$start;
        $CTGrefgen{$ctg}{$c}{'a'}=$start;
        $CTGrefgen{$ctg}{$c}{'b'}=$end;
        $c++;
}
close COORDFILE;



my $seqio2 = Bio::SeqIO->new('-format' =>  "fasta",
                            -file => "$refseq"
                                      );

my $io = Bio::SeqIO->new(-format => "fasta", -file => ">$out" );

while(my $seq = $seqio2->next_seq ) {
    my $len=$seq->length();
    my $id = $seq->id(); #   print "$id\t$len\n";
    my $seqs =$seq->seq();
#    unless (exists $CTGrefgen{$id}) {
   #     my $seq_obj = Bio::Seq->new(
#			      -seq => $seqs,
#			       -display_id  => "$id|R.C" );
#	my $io = Bio::SeqIO->new(-format => "fasta", -file => ">>$out" );  
#	$io->write_seq($seq_obj);
 #   }
    
    
    foreach my $reg (sort {$a <=> $b} keys %{$CTGrefgen{$id}}){
	my $substring = $seq->trunc($CTGrefgen{$id}{$reg}{'a'},$CTGrefgen{$id}{$reg}{'b'});

	my $seqO=$substring->seq();

	my $seq_obj = Bio::Seq->new(
			      -seq => $seqO,
			       #-display_id  => "$id-R.$reg|RegL.$refhash{$id}-FR.$CTGrefgen{$id}{$reg}{'a'}-$CTGrefgen{$id}{$reg}{'b'}" );
				-display_id  => "$id-R.$reg");
	#my $io = Bio::SeqIO->new(-format => "fasta", -file => ">>$out" );  
	$io->write_seq($seq_obj);
    }
}











