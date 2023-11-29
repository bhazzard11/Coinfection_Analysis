#!/usr/bin/perl -w                                                                                                                                                                                    

use strict;
use warnings;

##Takes list of unique snps and makes new list of positions for mpileup analysis##   

my $number_args = $#ARGV + 1;  ##Checking command line is getting what it expects

if ($number_args != 3) {  
    print "Wrong entry. Please enter Sample1_unique_snps.txt Sample1_unique_snps.txt Output.txt\n";
    exit;  
}  

open(my $IN1, "<".$ARGV[1]); #path to unique snp list for sample 1
open(my $IN2, "<".$ARGV[2]); #path to unique snp list for sample 2
open(my $OUT, ">".$ARGV[3]); #path to output

my %sites; #make list of sites where there is a snp and which sample it came from
while (my $NIH = <$IN1>){ #list for genotype 1
    chomp $NIH;
    my @fields = split /\t/, $NIH;
    my $position = $fields[0];
    $sites{$position} = 1;
}

while (my $CH = <$IN2>){ #list for genotype 2
    chomp $CH;
    my @fields = split /\t/, $CH;
    my $position = $fields[0];
    $sites{$position} =1;
}

my $count = 0;
for my $pos (sort keys %sites){
    $count++;
    my @fields = split /_/, $pos;
    my $num = $fields[3];
    my $chr = join("_", $fields[0], $fields[1], $fields[2]);
    print $OUT "$chr\t$num\n";
}

print "Number of positions = $count\n";
close($IN1);
close($IN2);
close($output);
