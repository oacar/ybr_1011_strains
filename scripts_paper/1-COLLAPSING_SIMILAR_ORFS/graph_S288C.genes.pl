#!/usr/bin/perl
use strict;
use Graph::Undirected;
my $g = Graph::Undirected->new;

#
#0 query
#1 hit
#2 id
#3 aln len
#4 5 mis gap 
#6 q start
#7 q stop
#8 h start
#9 h stop
#10 evalue
#11 points
#12 len query
#13 len hit
#14 cg hit
#15 cg query
#16 ratio aln len/query len
#17 ratio aln len/hit len
#18 max (16, 17)


#YAR010C	YPR137C-B	100	1323	0	0	1	1323	1	1323	0	1834.9	1323	5269	41.04	38.49	1	0.251091289	1

my %lens;

#parameter#

my $min_identity=95;
my $min_overl=0;
my $min_maxoverlperc=0.75;

my $param= "id.$min_identity.overl.$min_overl.overlper.$min_maxoverlperc.YRQ";

open T, ">graph.$ARGV[0].$param.dat";

open LIST, $ARGV[0];


while(<LIST>){
#    print $_;
    chomp $_;
    my @data=split(/\t/, $_);
    $lens{$data[0]}=$data[12];
    $lens{$data[1]}=$data[13];
    #print "$data[2] - $data[3] - $data[18]\n";
    my $ratio1=$data[3]/$data[12];
    my $ratio2=$data[3]/$data[13];
    my $ratio_max=max($ratio1,$ratio2);
    if (($data[2]>$min_identity)&&($data[3]>$min_overl)&&($ratio_max>$min_maxoverlperc)) {
        $g->add_edge($data[1],$data[0]);
    }
        
}


my %clusters;
my %b = $g->betweenness;

#foreach my $ver (sort {$a cmp $b} keys %b){
#    print "$ver $b{$ver}\n";
#}

my @c=$g->connected_components();

for (my $i=0;$i<scalar @c;$i++){
    my $num= scalar @{ $c[$i] };
print T "Cluster $i $num\n";
    $clusters{$i}{'num'}=$num;
    #$clusters{$i}{'longest'}=$c[$i][0];
    #$clusters{$i}{'lenofthelon'}=$lens{$c[$i][0]};
    $clusters{$i}{'betname'}=$c[$i][0];
    $clusters{$i}{'betmax'}=$b{$c[$i][0]};
    #$clusters{$i}{'shortest'}=$c[$i][0];
    #$clusters{$i}{'lenoftheshor'}=$lens{$c[$i][0]};
    
  for (my $j=0;$j<scalar @{ $c[$i] };$j++){
       print T "$c[$i][$j]\t$lens{$c[$i][$j]}\n";
       #if ($lens{$c[$i][$j]} > $clusters{$i}{'lenofthelon'}) {
        #$clusters{$i}{'longest'}=$c[$i][$j];
        #$clusters{$i}{'lenofthelon'} = $lens{$c[$i][$j]};
       #}
       
       my @num=split(/-/,$c[$i][$j]);
       my $code=$num[0];
       my @betnum=split(/-/,$clusters{$i}{'betname'});
       my $betcode=$betnum[0];
       #print "$code\t$betcode\n";
       
       if (($c[$i][$j]=~m/^[Y|R|Q]/g)||(($c[$i][$j]!~m/^[Y|R|Q]/g)&&($clusters{$i}{'betname'}!~m/^[Y|R|Q]/g))||(($c[$i][$j]=~m/^[Y|R|Q]/g)&&($clusters{$i}{'betname'}=~m/^[Y|R|Q]/g))){
       
      # if (($betcode>62)&&($code<61)) {
        #        $clusters{$i}{'betname'}=$c[$i][$j];
       #         $clusters{$i}{'betmax'} = $lens{$c[$i][$j]};
      # }
        if (($c[$i][$j]=~m/^[Y|R|Q]/)&&($clusters{$i}{'betname'}!~m/^[Y|R|Q]/)) {
                $clusters{$i}{'betname'}=$c[$i][$j];
                $clusters{$i}{'betmax'} = $lens{$c[$i][$j]};
       }
       
        if ($b{$c[$i][$j]} > $clusters{$i}{'betmax'}) {
        $clusters{$i}{'betname'}=$c[$i][$j];
        $clusters{$i}{'betmax'} = $lens{$c[$i][$j]};
       }
       if ($b{$c[$i][$j]} == $clusters{$i}{'betmax'}) {
        if ($lens{$c[$i][$j]} < $clusters{$i}{'betname'}){
            $clusters{$i}{'betname'}=$c[$i][$j];
            $clusters{$i}{'betmax'} = $lens{$c[$i][$j]};
        }
       }
       }
  }
}
close T;

open SU, ">graph_summary.$ARGV[0].$param.dat";




foreach my $clus_id (sort {$a <=> $b} keys %clusters){
    print SU "$clus_id\t$clusters{$clus_id}{'betname'}\t$clusters{$clus_id}{'num'}\n";
}





sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
}
