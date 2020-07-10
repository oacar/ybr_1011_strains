use strict;
use warnings;
#use lib '/Users/mdechiara/Downloads/BioPerl-1.6.924';
use Bio::SeqIO;

unless (-e $ARGV[0]){die "No external regions\n";}

my $seqio = Bio::SeqIO->new('-format' =>  "fasta",
                            -file => "$ARGV[0]");

my $out="$ARGV[0]";
$out=~s/.fasta/.afterN.fasta/;

while(my $seqi = $seqio->next_seq){
    my %count;
#    print "A\t";
    my $len=$seqi->length();
    my $id = $seqi->id(); 
    my $seq0=$seqi->seq();
    

    foreach my $str (split //, $seq0) {
        $count{$str}++;
    }
    my $percN=0;
    if (exists $count{'N'}){
    $percN=$count{'N'}/$len;
    }
 
#    print "$id - $len - $percN\n";    
    
    if ($percN<$ARGV[1]){
    my $seq_obj = Bio::Seq->new(
    -seq => $seq0 , 
    -display_id  => $id);

    my $io = Bio::SeqIO->new(-format => "fasta", -file => ">>$out" );
    $io->write_seq($seq_obj);
    }
}            

                                                                                                                               
