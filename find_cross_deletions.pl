#!/usr/bin/perl -w                                                                                                                                 \                    \
                                                                                                                                                                         
use strict;
use warnings;

###Finds drops in genome coverage between two species/subspecies based on comparison to same reference.
#Step 1: Create sorted bam file for each sample.
#Step 2: Run samtools mpileup for each sample (samtools mpileup -f path/to/Genome.fasta path/to/sorted.bam > path/to/output.mpileup)
#Step 3: Run this script either genome guided or by base pair chunks.

###Takes Command line input of either a genome gtf or the KB distance for window size, and two mpileup files to compare.
#Outputs tab delimited file with gene or window name, if using genome guided will output length of the gene if not outputs user defined size of window, number of reads at the position for each input file.


my $number_args = $#ARGV + 1;  ##Checking command line is getting what it expects

if ($number_args != 4) {  
    print "Wrong entry. Please enter NumberKB or GenomeAnotation.gff Sample1.mpileup Sample2.mpileup Output.txt\n";
    exit;  
}  

open(my $IN1, "<".$ARGV[1]); #Sample1
open(my $IN2, "<".$ARGV[2]); #Sample2
open(my $OUT, ">".$ARGV[3]); #Output file

if ($ARGV[0] =~ /^\d+$/){#If argument 0 is a number use this analysis (number is the size you want window to be in bps)
	my $windowKB = $ARGV[0]; #Set window size
	my %IN1_genes;
	my %gene_list;
	while (<$IN1>){ #Reads sample1 file and stores postion counts for sample1
    	chomp;
    	my @fields = split /\t/;
    	my $chr = $fields[0];
    	my $location = $fields[1];
    	my $reads = $fields[3];
    	my $window = int($location/$windowKB);
    	my $position = join ("_", $chr, $window);    
        if (defined $IN1_genes{$position}){
            $IN1_genes{$position} = ($IN1_genes{$position} + $reads);
        } else {
            $IN1_genes{$position} = $reads;
            $gene_list{$position} = 1;
        }                                                                                                                                                                  
	}
	close($IN1);
	
	my %IN2_genes;
	while (<$IN2>){ #Reads sample2 file and stores postion counts for sample2
    	chomp;
    	my @fields = split /\t/;
    	my $chr = $fields[0];
    	my $location = $fields[1];
    	my $reads = $fields[3];
    	my $window = int($location/$windowKB);
    	my $position = join ("_", $chr, $window);                                                                                                                                                                                                                         
        if (defined $IN2_genes{$position}){
            $IN2_genes{$position} = ($IN2_genes{$position} + $reads);
        } else {
           	$IN2_genes{$position} = $reads;
        	$gene_list{$position} = 1;
    	}                                                                                                                                                                 
	}
	close($IN2);
	
	foreach my $gene (sort keys %gene_list){ #Goes through list of positions and prints to file
    	my $IN1_value = 0;
    	my $IN2_value = 0;
    	print $OUT "$gene\t$windowKB\t";
    	if (defined $IN1_genes{$gene}){
        	my $IN1_value = $IN1_genes{$gene};
        	print $OUT "$IN1_value\t";
        	if (defined $IN2_genes{$gene}){
            	my $IN2_value = $IN2_genes{$gene};
            	print $OUT "$IN2_value\n";
        	} else {
            	print $OUT "$IN2_value\n";
        	}
    	} elsif (defined $IN2_genes{$gene}){
        	my $IN2_value = $IN2_genes{$gene};
        	print $OUT "$IN1_value\t$IN2_value\n";
    	} else {
        	print $OUT "$IN1_value\t$IN2_value\n";
    	}
	}
	close($OUT);
	exit;
} else { #If a number isn't given will run as genome guided
	open(my $gtf, "<".$ARGV[0]); #Opening genome file
	my %gene_list;
	my %gene_list_start;
	while (<$gtf>){     #Read genome file and store each gene and all the posible read positions for the gene                                                                                                                                                    
  	 	chomp;                                                                                                                                                              
   		my @fields = split /\t/;                                                                                                                                            
   	 	my $type = $fields[2];                                                                                                                                              
   		if ($type eq "mRNA"){                                                                                                                                               
        	my $chromosome = $fields[0];                                                                                                                                    
        	my $info = $fields[8];                                                                                                                                          
        	my $start = $fields[3];                                                                                                                                         
        	my $end = $fields[4];                                                                                                                                           
      		my $len = $end - $start;                                                                                                                                         
        	my $strand = $fields[6];                                                                                                                                        
        	my @collapse_info = split /;/, $info;                                                                                                                           
        	my $gene_id = $collapse_info[1];                                                                                                                                
        	my $gene_name = substr($gene_id, 7);                                                                                                                            
       		for my $num ($start..$end){                                                                                                                                      
            	my $window_name = join("_", $chromosome, $num);                                                                                                             
           		if (defined $gene_list_start{$window_name}){                                                                                                                 
            	} else {                                                                                                                                                    
            		$gene_list_start{$window_name} = $gene_name;                                                                                                            
                	$gene_list{$gene_name} = $len;                                                                                                                          
           		}                                                                                                                                                            
    	   	}                                                                                                                                                             
	    } else {}                                                                                                                                                           
	}                                                                                                                                                                       
	close($gtf);     
	             
	my %IN1_genes;
	while (<$IN1>){#Reads sample1 file and stores postion counts for sample1
    	chomp;
    	my @fields = split /\t/;
    	my $chr = $fields[0];
    	my $location = $fields[1];
    	my $reads = $fields[3];
    	my $position = $location;                                                                                                                  
   	 	if (defined $gene_list_start{$position}){                                                                                                                           
       		my $gene = $gene_list_start{$position};                                                                                                                          
        	if (defined $IN1_genes{$gene}){
            	$IN1_genes{$gene} = ($IN1_genes{$gene} + $reads);
        	} else {
            	$IN1_genes{$gene} = $reads;
            	$gene_list{$gene} = $gene_list{$gene};
        	}
    	}                                                                                                                                                                   
	}
	close($IN1);
	
	my %IN2_genes;
	while (<$IN2>){#Reads sample2 file and stores postion counts for sample2
    	chomp;
    	my @fields = split /\t/;
    	my $chr = $fields[0];
    	my $location = $fields[1];
    	my $reads = $fields[3];
    	my $position = $location;                                                                                                            
   	 	if (defined $gene_list_start{$position}){                                                                                                                           
       		my $gene = $gene_list_start{$position};                                                                                                                          
        	if (defined $IN2_genes{$gene}){
            	$IN2_genes{$gene} = ($IN2_genes{$gene} + $reads);
        	} else {
            	$IN2_genes{$gene} = $reads;
            	$gene_list{$gene} = $gene_list{$gene};
        	}
    	}                                                                                                                                                                   
	}
	close($IN2);
	
	foreach my $gene (sort keys %gene_list){ #Goes through list of genes and prints to file
    	my $IN1_value = 0;
    	my $IN2_value = 0;
    	print $OUT "$gene\t$gene_list{$gene}\t";
    	if (defined $IN1_genes{$gene}){
        	my $IN1_value = $IN1_genes{$gene};
        	print $OUT "$IN1_value\t";
        	if (defined $IN2_genes{$gene}){
            	my $IN2_value = $IN2_genes{$gene};
            	print $OUT "$IN2_value\n";
        	} else {
            	print $OUT "$IN2_value\n";
        	}
    	} elsif (defined $IN2_genes{$gene}){
        	my $IN2_value = $IN2_genes{$gene};
        	print $OUT "$IN1_value\t$IN2_value\n";
    	} else {
        	print $OUT "$IN1_value\t$IN2_value\n";
    	}
	}
	close($OUT);
}
exit;


	
