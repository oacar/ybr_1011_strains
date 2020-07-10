use strict;
#use lib '/Users/mdechiara/Downloads/BioPerl-1.6.924';
use Bio::SeqIO;

my %idlen;


my $seqio = Bio::SeqIO->new('-format' =>  "fasta",
			    -file => "$ARGV[0]"    #fasta to be check
				      );

my $c=0;
my %order;
my %order2;
while(my $seq = $seqio->next_seq ) {

    my $len1=$seq->length();
    my $id1 = $seq->id();
    $idlen{$id1}=$len1;
    $order{$id1}=$c;
    $order2{$c}=$id1;
    $c++;

}

my %keeped;
my %ends;
my $thr=100;
my $d=0;
open T, "$ARGV[1]"; #rem2.merge
while (<T>) {
    chomp $_;
    $d++;
    my @row=split(/\t/, $_);
    my $start= $row[1];
    my $stop=$row[2];
    my $sequ=$row[0];
    
    my $sequstart=1;
    
    if (exists $ends{$sequ}){
        $sequstart=1+$ends{$sequ};
    }
    my $totlen=$idlen{$sequ};
    my $upperlim=$totlen-$thr;
    if  ($start < $thr){
        $start = 1;
    }
    if ($stop>$upperlim) {
        $stop =$totlen;
    }
    
    print "$sequ\t$sequstart\t$start\t---$stop---$totlen\n";
    $ends{$sequ}=$stop;    
}

foreach my $reg (sort {$a cmp $b } keys %ends){
    if ($ends{$reg}<$idlen{$reg}) {
	print "$reg\t$ends{$reg}\t$idlen{$reg}\t---\n";
    }
}

foreach my $keep (sort {$a cmp $b } keys %idlen){
    unless (exists $ends{$keep}) {
	print "$keep\t1\t$idlen{$keep}\t---\n";
    }
    
}
