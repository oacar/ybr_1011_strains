use strict;


open T, "$ARGV[0]";
while (<T>) {
    chomp $_;

    my @row=split(/\t/, $_);
    my $start= $row[1];
    my $stop=$row[2];
    
    if ($stop-$start>200){
	print "$_\n";
    }
}
