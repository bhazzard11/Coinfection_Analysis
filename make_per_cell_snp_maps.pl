#!/usr/bin/perl -w                                                                                                                                                        \                                

use strict;
use warnings;

###Creates snp map for each cell, giving either N, C, or NA for each position with at least 200 cells##       
                                                                                                                                          
my $number_args = $#ARGV + 1;  ##Checking command line is getting what it expects

if ($number_args != 5) {  
    print "Wrong entry. Please enter Sample1_unique_snps.txt Sample2_unique_snps.txt cell_list.txt sample_mpileup Output.txt";
    exit;  
}  

open(my $NIH_snps, "<".$ARGV[0]); #path to unique snp list for sample 1
open(my $Chesson_snps, "<".$ARGV[1]); #path to unique snp list for sample 2
open(my $cell_list, "<".$ARGV[2]);#list of cells with cell name as first column (can come from count table)
open(my $sample, "<".$ARGV[3]);#sample mpileup
open(my $OUT, ">".$ARGV[4]); #path to output snp map for each cell


my %NIH_list; #Genotype 1 snp list
while (my $NIH = <$NIH_snps>){
    chomp $NIH;
    my @fields = split /\t/, $NIH;
    $NIH_list{$fields[0]} = $fields[2];
}
close($NIH_snps);

my %Ch_list; #Genotype 2 snp list
while (my $Ch = <$Chesson_snps>){
    chomp $Ch;
    my @fields = split /\t/, $Ch;
    $Ch_list{$fields[0]} = $fields[2];
}
close($Chesson_snps);

my %approved_cells; #list of cells
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

my %extra_ch = ("\$", 1, "^", 2, "+", 3, "-", 4, "]", 5, "\"", 6); #weird mpileup designations that are not needed
my %position_map;
my %cell_map;
my %dup_check;
my %location_count;
my %pos_list;
while (my $line = <$sample>){
    chomp $line;
    my @fields = split /\t/, $line;
    my $location = join("_", $fields[0], $fields[1]);
    my @reads = split //, $fields[4];
    my @name = split /,/, $fields[6];
    my $name_num = scalar(@name);
    if (scalar(@reads) >= 2){
        my %major_allele; #determine major allele at each location for each cell
        if (defined $NIH_list{$location}){ #if location matches one in genotype 1 list check if snp matches if not give other genotype
            my $strain = "NIH";
            my $snp = $NIH_list{$location};
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
                            $major_allele{$barcode}{$strain}++;
                        }
                    } elsif ($read eq "\."){
                        my @read_name = split /_/, $name[$counter];
                        my $barcode = $read_name[2];
                        my $UMI = $read_name[3];
                        my $other = "CH";
                        if (defined $dup_check{$barcode}{$UMI}){
                        } else {
                            $dup_check{$barcode}{$UMI} = 1;
                            $major_allele{$barcode}{$other}++;
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
		} elsif (defined $Ch_list{$location}){ #if location matches one in genotype 2 list check if snp matches if not give other genotype
            my $strain = "CH";
            my $snp = $Ch_list{$location};
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
                            $major_allele{$barcode}{$strain}++;
                        }
                    } elsif ($read eq "\."){
                        my @read_name = split /_/, $name[$counter];
                        my $barcode = $read_name[2];
                        my $UMI = $read_name[3];
                        my $other = "NIH";
                        if (defined $dup_check{$barcode}{$UMI}){
                        } else {
                            $dup_check{$barcode}{$UMI} = 1;
                            $major_allele{$barcode}{$other}++;
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
        foreach my $cell_al (sort keys %major_allele){ #Determines major allele at given location and gives either N or C designation for the genotype at that location.
            if (defined $approved_cells{$cell_al}){
                my $strain1;
                my $strain2;
                if (defined $major_allele{$cell_al}{"NIH"}){
                    $strain1 = $major_allele{$cell_al}{"NIH"};
                } else {
                    $strain1 = 0;
                }
                if (defined $major_allele{$cell_al}{"CH"}){
                    $strain2 = $major_allele{$cell_al}{"CH"};
                } else {
                    $strain2 = 0;
                }
                my $strain;
                if ($strain1 > $strain2){
                    $strain = "N"; #N for strain 1
                } else {
                    $strain = "C"; #C for strain 2
                }
                $position_map{$cell_al}{$location} = $strain; #defines cell, location with the strain of major allele
                $pos_list{$location}++; #counts number of cells with reads at that location
 			}
        }
    }
}
close($sample);

my @snp_list; #creates list of snps for output based on the number reads at each location only keep locations with at least 10 cells being defined at that location
foreach my $loc (sort keys %pos_list){
    if ($pos_list{$loc} >= 10){
        push @snp_list, $loc;
    }
}

print $output join("\t","Cell",@snp_list), "\n";
foreach my $cell_name (sort keys %position_map){ #Outputs for each cell the allele at each position, fills in NA if cell does not have a defined allele at that location
    print $output "$cell_name\t";
    my @data;
    foreach my $allele (sort @snp_list){
        if (defined $position_map{$cell_name}{$allele}){
            push @data, ($position_map{$cell_name}{$allele});
        } else {
            push @data, "NA";
        }
    }
    print $output join("\t", @data), "\n";
}
close($output);