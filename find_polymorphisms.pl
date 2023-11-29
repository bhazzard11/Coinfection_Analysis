#!/usr/bin/perl -w                                                                                                                                                                                   \
                                                                                                                                                                                                      

use strict;
use warnings;

###reads mpileup output to determine ratio of polymorphesisms, used to determine if sample is polyclonal###                                                                                           
#Takes input of mpileup file and output name: perl make_site_list.pl path/to/sample.mpileup path/to/OUTNAME                                                                                                              
my $number_args = $#ARGV + 1;  

if ($number_args != 2) {  
    print "Wrong entry. Please enter paired mpileup and output name.txt\n";
    exit;  
}  

open(my $IN2, "<".$ARGV[0]); #mpileup
open (my $output, ">", "$ARGV[1]_polyclonal_ratio.txt"); #path to output for the ratio of alleles at that position.                                                               
open (my $snps, ">", "$ARGV[1]_snp_list.txt");   #path to output for the list of positions for each snp compared to the reference.                                                                               
open (my $masked, ">", "$ARGV[1]_mask_list.txt");    #path to output for list of positions that are undeterminable.                                                                           


my %alt_list = ("A", 1, "T", 2, "G", 3, "C", 4, "a", 5, "t", 6, "g", 7, "c", 8); #possible snp list
my %mask;
my $mask_count = 0;
while (my $read = <$mpileup>){ #read mpileup
    chomp $read;
    my @fields = split /\t/, $read; #split line on tab
    my $len = scalar(@fields);
#    print "$len\n";                                                                                                                                                                                  
    my $chrm = $fields[0]; #chromosone
    my $pos = $fields[1]; #position
    my $norm = $fields[2]; #base in reference genome
    my $num = $fields[3]; #number of reads at that position
    my $type = $fields[4]; #mpileup designation for each read at that position
#    my $len = $fields[2];                                                                                                                                                                            
    my $total = 0;
    my $alt = 0;
    my $snp;
    my @list = ($chrm, $pos);
    my $name = join("_", @list);
    if ($num >= 50){ #make sure at least 50 reads at position
		 if ($len == 6){ 
            for my $c (split //, $type){
                $total++;
                if ($alt_list{$c}){
                    $snp = $c;
                    $alt++;
                } else {}
            }
        }
        if ($alt >= 1){
            my $ratio = $alt/$total;
            print $output "$name\t$ratio\n"; #print ratio of allele1:allele2 for that position across all cells
            if ($ratio >= 0.3){
                if ($ratio <= 0.8){
                    print $masked "$name\n";
                } elsif ($ratio >= .9){
                    print $snps "$name\t$norm\t$snp\t$num\n"; #print snp location, the reference allele, the snp allele, and the number of reads at that position
                }
            }
        } else {
            my $ratio = 0;
            print $output "$name\t$ratio\n"; #print ratio of allele1:allele2 for that position across all cells
        }
    }
}

close($mpileup);
close($output);
close($snps);
