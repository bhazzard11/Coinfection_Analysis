#!/usr/bin/perl -w                                                                                                                                                                                   \
                                                                                                                                                                                                      
use strict;
use warnings;

###Count NIH vs Chesson snps for each cell in a single cell rna-seq sample###  
# Outputs for each cell from a single cell rna-seq experiment the number of alleles that originate from two different genotypes.                                                                                                                                                           
my $number_args = $#ARGV + 1;  ##Checking command line is getting what it expects

if ($number_args != 6) {  
    print "Wrong entry. Please enter Sample1_unique_snps.txt Sample2_unique_snps.txt cell_list.txt sample_mpileup Output.txt CELLPREFIX";
    exit;  
}  

open(my $NIH_snps, "<".$ARGV[0]); #path to unique snp list for sample 1
open(my $Chesson_snps, "<".$ARGV[1]); #path to unique snp list for sample 2
open(my $cell_list, "<".$ARGV[2]);#list of cells with cell name as first column (can come from count table)
open(my $sample, "<".$ARGV[3]);#sample mpileup
open(my $OUT, ">".$ARGV[4]); #path to output
my $samplename = $ARGV[5]; #cell prefix

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

my %extra_ch = ("\$", 1, "^", 2, "+", 3, "-", 4, "]", 5, "\"", 6); #possible weird symbols for mpileup
my %CH_cells;
my %NIH_cells;
my %dup_check;
while (my $line = <$sample>){ #look through mpileup 
    chomp $line;
    my @fields = split /\t/, $line;
    my $location = join("_", $fields[0], $fields[1]); #location of read chromosome_position
    my @reads = split //, $fields[4]; #mpileup code
    my @name = split /,/, $fields[6]; #read name
    my $name_num = scalar(@name); #go through all cell names/codes
    if (scalar(@reads) >= 2){
        my %major_allele;
		 if (defined $NIH_list{$location}){ #If snp location is in Genotype 1 list check snp
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
    	} elsif (defined $Ch_list{$location}){ #If snp location is in Genotype 2 list check snp
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
                    } elsif ($read eq "."){
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
        foreach my $cell_al (sort keys %major_allele){ #If more then one read per cell at a location determine major allele at that location for each cell
            if (defined $approved_cells{$cell_al}){
                my $read_count = 0;
                foreach my $allele (sort keys %{$major_allele{$cell_al}}){
                    $read_count = $read_count + $major_allele{$cell_al}{$allele};
                }
                my $highest = (sort {$major_allele{$cell_al}{$a} <=> $major_allele{$cell_al}{$b}} keys %{$major_allele{$cell_al}})[0];                                                                                                                                                           
                    if ($highest eq "NIH"){
                        $NIH_cells{$cell_al}++;
                    } elsif ($highest eq "CH"){
                        $CH_cells{$cell_al}++;
                    }
#               }                                                                                                                                                                                     
            }
        }
    }
}
close($sample);

my $cell_count;
print $output "Cell\tChesson\tNIH\n"; #Outputs number of reads each cell has for each strain
foreach my $cell (sort keys %approved_cells){
    $cell_count++;
    print $output "$samplename$cell\t";
    if (defined $CH_cells{$cell}){
        print $output "$CH_cells{$cell}\t";
    } else {
        print $output "0\t";
    }
    if (defined $NIH_cells{$cell}){
        print $output "$NIH_cells{$cell}\n";
    } else {
        print $output "0\n";
    }
}
close($output);

print "Cell Count: $cell_count\n"; #Should match the number of cells in the count table


