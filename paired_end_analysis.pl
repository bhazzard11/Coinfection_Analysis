#!/usr/bin/perl -w                                                                                                                                                      \
                                                                                                                                                                         

use strict;
use warnings;

##Takes paired end bam file and determines number of insertions, deletions, and inversions based on samflag.
#Outputs text file with number of each type and location.

my $number_args = $#ARGV + 1;  

if ($number_args != 2) {  
    print "Wrong entry. Please enter paired end BAM file and output name.txt\n";
    exit;  
}  

open(my $BAM, "-|", "samtools", "view",$ARGV[0]); #Opens BAM file
open(my $output, ">".$ARGV[1]); #Opens output

my %del_flags = ("145", 1, "97", 2, "81", 3, "161", 4); #Samflag list for deletion
my %inv_flags = ("177",1, "113",2, "65",3, "129",4); #Samflag list for inversions
my %ins_flags = ("73", 1, "133", 2, "69", 3, "137", 4); #Samflag list for insertions
my %deletions;
my %duplication;
my %inversion;
my %insertion;
while (<$BAM>){ #Reads bam file and counts all the samflags
    chomp;
    my @fields = split /\t/;
    my $sam_flag = $fields[1];
    my $chr = $fields[2];
    my $position = $fields[3];
    my $window = int($position/1000);
    my $read_pair = $fields[6];
    my $read_distance = $fields[8];
    my $location = join("_", $chr, $window);
    if ($del_flags{$sam_flag}){
        if ($read_distance >=1000){
            if (defined $deletions{$location}){
                $deletions{$location}++;
            } else {
                $deletions{$location} = 1;
            }
		}
    } elsif ($inv_flags{$sam_flag}){
        if (defined $inversion{$location}){
                $inversion{$location}++;
        } else {
            $inversion{$location} = 1;
        }
    } elsif ($ins_flags{$sam_flag}){
		if (defined $insertion{$location}){
            $insertion{$location}++;
        } else {
            $insertion{$location} = 1;
        }
    }
}
close($BAM);

foreach my $del (sort keys %deletions){ #Prints deletions
    if ($deletions{$del} > 10){
        print $output "$del\tDeletion\t$deletions{$del}\n";
    }
}
foreach my $ins (sort keys %insertion){ #Prints insertions
    if ($insertion{$ins} > 10){
        print $output "$ins\tInsertion\t$insertion{$ins}\n";
    }
}

foreach my $inv (sort keys %inversion){ #Prints inversions
    if ($inversion{$inv} > 10){
        print $output "$inv\tInversion\t$inversion{$inv}\n";
    }
}
close($output);
exit;
