#!/usr/bin/perl -w
use strict;
use POSIX;
#use lib '/Users/mdechiara/Downloads/BioPerl-1.6.924';
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;

my $blastpath="/ccc/work/cont007/fg0006/fg0006/Programs/ncbi-blast-2.2.30+/bin";
my $bedpath="/ccc/cont007/home/fg/fg/products/bedtools-2.24.0/bin";


my $Tid=95; #threshold for identity percentage
my $Tlen=200; #threshold for alignment length
my $subcheck=500;#size to slide the sequence
my $ref="sgd"; #blastdb of the reference

my $fasta_library = $ARGV[0]; #assembly
my $firstphase = "$blastpath/blastn -query $fasta_library -db $ref -outfmt '6 std qlen slen' -perc_identity 90 -reward 1 -penalty -5 -gapopen 5 -gapextend 5 -no_greedy > $ARGV[1]";
print "$firstphase\n";
system($firstphase);

my @m8 = qx "cat $ARGV[1] |wc -l";
$m8[0]=~s/\s//g;
print "$m8[0]\n";
#}



my $database = Bio::DB::Fasta->new("$fasta_library") or die "Failed to creat Fasta DP object on fasta library\n";

#AGN_2-6949      chr3    99.39   3462    19      1       360     3819    307747  304286  0.0      6658   3819    316620

my %eliminate;
my $d=0;

my %idlen;




#----------
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
#--------------------------

open(my $log, ">", "$ARGV[0].log") or die $!;
my $nol=0;
open T, "$ARGV[1]"; # blastn -outfmt 6 table
while (<T>) {
    $nol++;
    print $log "Line  $nol out of $m8[0]\n";
    chomp $_;
    my @row=split(/\s+/, $_);
    my $query= $row[0];
    my $hit= $row[1];
    my $id=$row[2];
    my $len=$row[3];
    my $pos1A=$row[6];
    my $pos1B=$row[7];
    my $pos2A=$row[8];
    my $pos2B=$row[9];
    my $qlen=$row[12];
    my $slen=$row[13];
    if (($id>$Tid)&&($qlen<=100)){
            $eliminate{$query}{$d}{'A'}=1;
            $eliminate{$query}{$d}{'B'}=$qlen;
            $d++;
    }
    elsif (($id>$Tid)&&($len<=$Tlen)&&($qlen-$len<75)){
            $eliminate{$query}{$d}{'A'}=1;
            $eliminate{$query}{$d}{'B'}=$qlen;
            $d++;
    }
    elsif (($id>$Tid)&&($len>$Tlen)){
        if ($len<1500) {
            $eliminate{$query}{$d}{'A'}=$pos1A;
            $eliminate{$query}{$d}{'B'}=$pos1B;
            $d++;
        }
        else{
            my $last=ceil($len/$subcheck);
            my @array = (1..$last);
            print "$query\n";
            my $slice = $database->get_Seq_by_id($query);
         #   print "id", $slice->id(), "\n";
            my $i;
            my $e;
            my $in=$pos1A;
            my $end=$pos1A+$subcheck;            
            while ($end - $pos1B <$subcheck ) {
                my $so;
                unless ($end > $pos1B){
                    $so = $slice->trunc($in,$end);
                    $i=$in;
                    $e=$end;
                    $in=$in+$subcheck;
                    $end=$end+$subcheck;
                }
                else {
                    $end=$pos1B;
                    $i=$in;
                    $e=$end;
                    $so = $slice->trunc($in,$end);
                    $end=$end+(2*$subcheck);
                }
                my $seq=$so->seq();
                my $cmd1="echo \">$query\n$seq\" | $blastpath/blastn  -db $ref -outfmt '6 std' -perc_identity 90 -reward 1 -penalty -5 -gapopen 5 -gapextend 5 -no_greedy 2>>err.log";
                #print "1 $cmd1\n";
                #system($cmd1);
                #print "2 $cmd1\n";
                #print "IN - END -- $i - $e\n";
                my $lensub=$e-$i;
                my @out = qx "$cmd1";
                foreach my $j (@out){
                    chomp $j;
                    #print "OUT $j\n";
                    my @innerrow=split(/\s+/, $j);
                    if (($innerrow[2]>$Tid)&&($innerrow[3]>0.75*$lensub)){
                        $eliminate{$query}{$d}{'A'}=$i-1+$innerrow[6];
                        $eliminate{$query}{$d}{'B'}=$i-1+$innerrow[7];
                        $d++;
                    }
                }
            }
        }
    }
}

close T;


open(BED, ">", "$ARGV[1].temp.rem.bed" ) or die "not opening temp";
foreach my $ord (sort {$a <=> $b} keys %order2){
    my $to_be_elim=$order2{$ord};
        foreach my $k (sort {$a <=> $b} keys %{$eliminate{$to_be_elim}}){
            print BED "$to_be_elim\t$eliminate{$to_be_elim}{$k}{'A'}\t$eliminate{$to_be_elim}{$k}{'B'}\t$idlen{$to_be_elim}\n";
        }
}
close BED;
system ("cat $ARGV[1].temp.rem.bed | $bedpath/sortBed -i - |$bedpath/mergeBed -i - -d 100 > $ARGV[1].rem.bed 2>>err.log");
unlink ("$ARGV[1].temp.rem.bed");

system("perl invert_bed2.pl $ARGV[0] $ARGV[1].rem.bed > $ARGV[1].temp.keep.bed");
system ("cat $ARGV[1].temp.keep.bed | $bedpath/sortBed -i - |$bedpath/mergeBed -i - -d 100 > $ARGV[1].keep.bed 2>>err.log");
#unlink ("$ARGV[1].temp.rem.bed");
system ("perl noshorterthen200.pl $ARGV[1].keep.bed >$ARGV[1].kept.bed 2>>err.log");
system("perl bedRegions2.pl $ARGV[1].kept.bed $ARGV[0]");



