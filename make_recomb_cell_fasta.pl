#!/usr/bin/perl -w                                                                                                                                                 \                                 \
                                                                                                                                                                                                      

use strict;
use warnings;

###Uses top most expressed positions to determine number of recombination events in sample##                                                                                                         \
# Outputs fasta of alleles at each position with at least 200 cells for each cell#

my $number_args = $#ARGV + 1;  ##Checking command line is getting what it expects

if ($number_args != 6) {  
    print "Wrong entry. Please enter Sample1_unique_snps.txt Sample2_unique_snps.txt cell_list.txt sample_mpileup position_counts.txt Output.fasta";
    exit;  
}  

open(my $NIH_snps, "<".$ARGV[0]); #path to unique snp list for sample 1
open(my $Chesson_snps, "<".$ARGV[1]); #path to unique snp list for sample 2
open(my $cell_list, "<".$ARGV[2]);#list of cells with cell name as first column (can come from count table)
open(my $sample, "<".$ARGV[3]);#sample mpileup
open(my $position_list, "<".$ARGV[4]);#sample counts per snp positions
open(my $OUT, ">".$ARGV[5]); #path to output counts per position intermediate file
                                                                                                                                                                                                     
my %all_pos;
my %NIH_list;
while (my $NIH = <$NIH_snps>){
    chomp $NIH;
    my @fields = split /\t/, $NIH;
    $NIH_list{$fields[0]} = $fields[2];
    $all_pos{$fields[0]} = 1;
}
close($NIH_snps);

my %Ch_list;
while (my $Ch = <$Chesson_snps>){
    chomp $Ch;
    my @fields = split /\t/, $Ch;
    $Ch_list{$fields[0]} = $fields[2];
    $all_pos{$fields[0]} = 1;
}
close($Chesson_snps);

my %approved_cells;
my $rows = 0;
while (my $list = <$cell_list>){
    chomp $list;
    my @fields = split /\t/, $list;
    if ($rows >= 1){
        my $barcode = $fields[0];
        $approved_cells{$barcode} = 1;
    } else {
        $rows++;
    }
}
close($cell_list);

my %top_pos;
while (my $top = <$position_list>){
    chomp $top;
    my @fields = split /\t/, $top;
    if ($fields[1] >= 200){
        $top_pos{$fields[0]} = 1;
    }
}
close($position_list);

my %extra_ch = ("\$", 1, "^", 2, "+", 3, "-", 4, "]", 5, "\"", 6);
my %position_map;
my %cell_map;
my %dup_check;
my %location_count;
my %pos_list;
my %cell_strain;
while (my $line = <$sample>){
    chomp $line;
    my @fields = split /\t/, $line;
    my $location = join("_", $fields[0], $fields[1]);
    my @reads = split //, $fields[4];
    my @name = split /,/, $fields[6];
    my $name_num = scalar(@name);
    my $base = $fields[2];
    if (defined $top_pos{$location}){
        my %major_allele;
        my %major_strain;
        if (defined $NIH_list{$location}){
            my $snp = $NIH_list{$location};
            my $strain = "NIH";
            my $counter = 0;
            my $skip = 0;
            for my $read (@reads){
                if ($skip >= 1){
                    $skip = $skip - 1;
                } else {
                    if ($read eq $snp){
                        my @read_name = split /_/, $name[$counter];
                        my $barcode = $read_name[2];
                        my $UMI = $read_name[3];
                        if (defined $dup_check{$barcode}{$UMI}){
                        } else {
                            $dup_check{$barcode}{$UMI} = 1;
                            $major_allele{$barcode}{$read}++;
                            $major_strain{$barcode}{$strain}++;
                        }
 					} elsif ($read eq "\."){
                        my @read_name = split /_/, $name[$counter];
                        my $barcode = $read_name[2];
                        my $UMI = $read_name[3];
                        my $other = "CH";
                        if (defined $dup_check{$barcode}{$UMI}){
                        } else {
                            $dup_check{$barcode}{$UMI} = 1;
                            $major_allele{$barcode}{$base}++;
                            $major_strain{$barcode}{$other}++;
                        }
                    } elsif (defined $extra_ch{$read}){
                        $counter = $counter- 1;
                    } elsif ($read =~ /\d+/){
                        $counter = $counter- 1;
                        $skip = $read;
                    }
                    $counter++;
                }
            }
        } elsif (defined $Ch_list{$location}){
         	my $snp = $Ch_list{$location};
            my $strain = "CH";
            my $counter = 0;
            my $skip = 0;
            for my $read (@reads){
                if ($skip >= 1){
                    $skip = $skip - 1;
                } else {
                    if ($read eq $snp){
                        my @read_name = split /_/, $name[$counter];
                        my $barcode = $read_name[2];
                        my $UMI = $read_name[3];
                        if (defined $dup_check{$barcode}{$UMI}){
                        } else {
                            $dup_check{$barcode}{$UMI} = 1;
                            $major_allele{$barcode}{$read}++;
                            $major_strain{$barcode}{$strain}++;
                        }
                    } elsif ($read eq "\."){
                        my @read_name = split /_/, $name[$counter];
                        my $barcode = $read_name[2];
                        my $UMI = $read_name[3];
                        my $other = "NIH";
                        if (defined $dup_check{$barcode}{$UMI}){
                        } else {
                            $dup_check{$barcode}{$UMI} = 1;
                            $major_allele{$barcode}{$base}++;
                            $major_strain{$barcode}{$other}++;
                        }
                    } elsif (defined $extra_ch{$read}){
                        $counter = $counter- 1;
                    } elsif($read =~ /\d+/){
                        $counter = $counter- 1;
                        $skip = $read;
                    }
                    $counter++;
                }
            }
        }
		foreach my $cell_al (sort keys %major_allele){
            if (defined $approved_cells{$cell_al}){
                my $highest = (sort {$major_allele{$cell_al}{$a} <=> $major_allele{$cell_al}{$b}} keys %{$major_allele{$cell_al}})[0];
                $cell_map{$cell_al}{$location} = $highest;
                my $highest_strain = (sort {$major_strain{$cell_al}{$a} <=> $major_strain{$cell_al}{$b}} keys %{$major_strain{$cell_al}})[0];
                if ($highest_strain eq "NIH"){
                    $cell_strain{$cell_al}{$location} = "NIH";
                    my $name = "NIH";
                    $location_count{$location}{$name}++;
                } elsif ($highest_strain eq "CH"){
                    $cell_strain{$cell_al}{$location} = "CH";
                    my $name = "CH";
                    $location_count{$location}{$name}++;
                }

            }
        }
    }
}
close($sample);

my %recomb_cells;
foreach my $cell (sort keys %cell_strain){
    my $NIH_count = 0;
    my $total_count = 0;
    foreach my $nuc (sort keys %{$cell_strain{$cell}}){
        if ($cell_strain{$cell}{$nuc} eq "NIH"){
            $NIH_count++;
            $total_count++;
        } else {
            $total_count++;
        }
    }
    my $ratio = $NIH_count/$total_count;
    if ($ratio <= 0.9){
        $recomb_cells{$cell} = 1;
    }
}

foreach my $barcode (sort keys %cell_map){
    my @seq;
    my $N_count = 0;
    my $total_pos = 0;
    foreach my $allele (sort keys %top_pos){
        if (defined $cell_map{$barcode}{$allele}){
            my $letter = uc $cell_map{$barcode}{$allele};
            push @seq, $letter;
        } else {
            push @seq, "N";
            $N_count++;
        }
        $total_pos++;
    }
    my $math = $N_count/$total_pos;
    if ($math <= 0.60){
        print $output ">$barcode\n@seq\n";
    }
}
close($output);


