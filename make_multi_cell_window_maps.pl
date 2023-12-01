#!/usr/bin/perl -w                                                                                                                                               \                                         

use strict;
use warnings;

###Breaks up genome into windows, size determined by user. Compresses allele assignment into dominant allele in the window.##
# Outputs file with allele assignment for each window for each cell where at least 50% of windows have a allele assignment (1 = allele 1, 2 =  allele 2, 3 = no allele assignment)#
# Outputs file with counts per cell of how many alleles make up window#
# Outputs file with allele pattern across all windows with missing values being "U", outputs as fasta format that can be used in MEGA to generate phylogenetic trees#                                                                                                                                       

my $number_args = $#ARGV + 1;  ##Checking command line is getting what it expects

if ($number_args != 6) {  
    print "Wrong entry. Please enter Sample1_unique_snps.txt Sample2_unique_snps.txt cell_list.txt sample_mpileup OUTPREFIX WindowSize";
    exit;  
}  

open(my $NIH_snps, "<".$ARGV[0]); #path to unique snp list for sample 1
open(my $Chesson_snps, "<".$ARGV[1]); #path to unique snp list for sample 2
open(my $cell_list, "<".$ARGV[2]);#list of cells with cell name as first column (can come from count table)
open(my $sample, "<".$ARGV[3]);#sample mpileup
my $samplename = $ARGV[4]; #cell prefix

open (my $output, ">", "$samplename_window_snp_assignments.txt");
open (my $output2, ">", "$samplename_window_snp_counts.txt");
open (my $output3, ">", "$samplename_window_snp_assignments.fasta");

#my $window_size = 200000; #default window size 200000 (what was used in Hazzard et al. 2024)
my $window_size = $ARGV[5]; 

my %NIH_list;
while (my $NIH = <$NIH_snps>){
    chomp $NIH;
    my @fields = split /\t/, $NIH;
    $NIH_list{$fields[0]} = $fields[2];
}
close($NIH_snps);

my %Ch_list;
while (my $Ch = <$Chesson_snps>){
    chomp $Ch;
    my @fields = split /\t/, $Ch;
    $Ch_list{$fields[0]} = $fields[2];
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

my %extra_ch = ("\$", 1, "^", 2, "+", 3, "-", 4, "]", 5, "\"", 6);
my %position_map;
my %cell_map;
my %dup_check;
my %location_count;
my %win_map;
my %win_list;
my %genotype_homo;
my %genotype_het;
my %snp_win;
while (my $line = <$sample>){
    chomp $line;
    my @fields = split /\t/, $line;
    my $window = int($fields[1]/$window_size); #200,000
    my $location = join("_", $fields[0], $fields[1]);
    my $location_win = join("_", $fields[0],$window);
    $snp_win{$location_win}{$location} = 1;
    #$win_list{$location_win}++;                                                                                                                                                                           
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
                if ($location_win =~ "Transfer"){
                } else {
                $win_list{$location_win}++;                                                                                                                                         
                my $highest = (sort {$major_allele{$cell_al}{$a} <=> $major_allele{$cell_al}{$b}} keys %{$major_allele{$cell_al}})[0];
                if ($highest eq "NIH"){
                    $cell_map{$cell_al}{$location} = "NIH";
                    my $name = "NIH";
                    $win_map{$cell_al}{$location_win}{$name}++;
                    $location_count{$location_win}{$name}++;
                    $position_map{$fields[0]}{$window}{$name}++;
                } elsif ($highest eq "CH"){
                    $cell_map{$cell_al}{$location} = "CH";
                    my $name = "CH";
                    $win_map{$cell_al}{$location_win}{$name}++;
                    $location_count{$location_win}{$name}++;
                    $position_map{$fields[0]}{$window}{$name}++;
				}
            }
        }
    }
}
close($sample);

my %recomb_cells;
foreach my $cell (sort keys %cell_map){
    my $NIH_count = 0;
    my $total_count = 0;
    foreach my $nuc (sort keys %{$cell_map{$cell}}){
        if ($cell_map{$cell}{$nuc} eq "NIH"){
            $NIH_count++;
            $total_count++;
        } else {
            $total_count++;
        }
    }
    my $ratio = $NIH_count/$total_count;
    if ($ratio <= 0.9){
        if ($total_count >= 200){
            $recomb_cells{$cell} = 1;
        }
    }
}

my %recomb_pats;
my %het_win;
my %cell_win_count;
my %snp_alt_freq;
my %snp_tot;
my %cell_order;
my @allwindows = sort keys %win_list;
foreach my $cell (sort keys %win_map){
    if (defined $recomb_cells{$cell}){
#       print $output "\n$cell";                                                                                                                                                                               
        my @pattern;
        my $cell_win_num = scalar keys %{$win_map{$cell}};
        foreach my $window (sort keys %win_list){
#           if ($win_list{$window} >= 100){                                                                                                                                                                    
            if ($window =~ "Transfer"){
            } else {
				if (defined $win_map{$cell}{$window}){
                    $cell_win_count{$cell}++;
                    my $cell_win_snp_count;
                    my $NIH_count = 0;
                    my $CH_count = 0;
                    foreach my $allele (sort keys %{$win_map{$cell}{$window}}){
                        if ($allele eq "NIH"){
                            $NIH_count = $win_map{$cell}{$window}{$allele};
                        } elsif ($allele eq "CH"){
                            $CH_count = $win_map{$cell}{$window}{$allele};
                        }
                        $cell_win_snp_count = $NIH_count + $CH_count;
                    }
                    if ($NIH_count > $CH_count){
                        #push @pattern, "N";                                                                                                                                                                   
                        my $snp_sum = $NIH_count + $CH_count;
                        if ($snp_sum >= 3){
                            $snp_alt_freq{$cell}{$window} = $CH_count;
                            $snp_tot{$cell}{$window} = $snp_sum;
                            my $snp_per = ($CH_count/$snp_sum)*100;
                            if ($snp_per <= 33){
                                push @pattern, "A";
                                my $color = 1;
                                $recomb_pats{$cell}{$window} = $color;
                                $het_win{$window}++;
                            } else {
                                my $color = 3;
                                push @pattern, "N";
                                $recomb_pats{$cell}{$window} = $color;
                            }
						} else {
                            push @pattern, "N";
                            my $color = 3;
                            $recomb_pats{$cell}{$window} = $color;
                        }
                    } elsif ($CH_count > $NIH_count){
                        #push @pattern, "C";                                                                                                                                                                   
                        my $snp_sum = $NIH_count + $CH_count;
                        $snp_alt_freq{$cell}{$window} = $NIH_count;
                        $snp_tot{$cell}{$window} = $snp_sum;
                        my $snp_per = ($NIH_count/$snp_sum)*100;
                        if ($snp_sum >= 3){
                            if ($snp_per <= 33){
                                push @pattern, "C";
                                my $color = 2;
                                $recomb_pats{$cell}{$window} = $color;
                                $het_win{$window}++;
                            } else {
                                push @pattern, "N";
                                my $color = 3;
                                $recomb_pats{$cell}{$window} = $color;
                            }
						} else {
                            push @pattern, "N";
                            my $color = 3;
                            $recomb_pats{$cell}{$window} = $color;
                        }
                    } else {
                        push @pattern, "N";
                    }
                } else {
                    push @pattern, "N";
                }
            }
        }
        my $strain = join(",", @pattern);
        print $output3 ">$cell\n@pattern\n";
        $cell_order{$strain} = $cell;
    }
}
close($output3);

print $output "cell_window\tsnps\tminer_al_freq\n";
foreach my $cell_win (sort keys %snp_alt_freq){
    foreach my $win_dow (sort keys %{$snp_alt_freq{$cell_win}}){
        my $cell_name = join("_", $cell_win, $win_dow);
        print $output "$cell_name\t$snp_tot{$cell_win}{$win_dow}\t$snp_alt_freq{$cell_win}{$win_dow}\n";
    }
}
close($output);

my @het_win_list = sort keys %het_win;
print $output2 "Cell\t", join("\t",@het_win_list);
my $cell_count = 0;
foreach my $strain_sort (sort keys %cell_order){
    my $recomb_cell = $cell_order{$strain_sort};
    $cell_count++;
    print $output2 "\n$recomb_cell";
    foreach my $recomb_win (sort keys %het_win){
        if (defined $recomb_pats{$recomb_cell}{$recomb_win}){                
        print $output2 "\t$recomb_pats{$recomb_cell}{$recomb_win}";
        } else {
            print $output2 "\t3";
        }
    }                                                                                                                                                                                                    
}
close($output2);



