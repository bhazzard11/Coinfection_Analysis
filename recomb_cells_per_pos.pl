#!/usr/bin/perl -w                                                                                                                                                                                        \
                                                                                                                                                                                                           
use strict;
use warnings;

###Counts number of cells have a read at a given position##      
# Creates an intermediate file that give the number of cells with reads at a given snp position. 
                                                                                                                                           
my $number_args = $#ARGV + 1;  ##Checking command line is getting what it expects

if ($number_args != 5) {  
    print "Wrong entry. Please enter Sample1_unique_snps.txt Sample2_unique_snps.txt cell_list.txt sample_mpileup Output.txt";
    exit;  
}  

open(my $NIH_snps, "<".$ARGV[0]); #path to unique snp list for sample 1
open(my $Chesson_snps, "<".$ARGV[1]); #path to unique snp list for sample 2
open(my $cell_list, "<".$ARGV[2]);#list of cells with cell name as first column (can come from count table)
open(my $sample, "<".$ARGV[3]);#sample mpileup
open(my $OUT, ">".$ARGV[4]); #path to output counts per position intermediate file


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

my %approved_cells; #List of cell names
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

my %extra_ch = ("\$", 1, "^", 2, "+", 3, "-", 4, "]", 5, "\"", 6);#weird mpileup designations that are not needed
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
        my %major_allele;
        if (defined $NIH_list{$location}){
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
        } elsif (defined $Ch_list{$location}){
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
        foreach my $cell_al (sort keys %major_allele){
            if (defined $approved_cells{$cell_al}){
                $pos_list{$location}++;
            }
        }
    }
}
close($sample);

print $output "Position\tCells\n";
my $over_count;
foreach my $allele (sort keys %pos_list){
    print $output "$allele\t$pos_list{$allele}\n";
}



