# Coinfection_Analysis
Scripts used to analyze experimental coinfection of two strains of distantly related strains of the same species.

For Plasmodium:


For Cryptosporidium:
Overview of what I did and what each code does:
1.	Run mpileup on each sorted bam file (samtools mpileup -f path/to/Genome.fasta path/to/sorted.bam > path/to/output.mpileup)
2.	Ran find_cross_deletions.pl (perl find_cross_deletions.pl 1000 mapped_Cp_isra_RG.mpileup mapped_Cp_tom_RG.mpileup snp_comparison_Cp.txt). The output is a txt file with the: name of the window by chromosome_window, size of the window (I just ran at 1000bp but this can be changed), Input1 number of reads, Input2 number of reads. I didn’t do the math in what’s attached but you can do this yourself, all your looking for is a large discrepancy in the number of reads at a given location between the samples, excel can do this.
3.	Ran paired_end_analysis (perl paired_end_analysis.pl mapped_Cp_tom_RG.bam mapped_Cp_tom_RG_paired_end_analysis.txt). First input is the bam file (it doesn’t need to be sorted) then output. File output is location, what type it is (Insertion, Inversion, or Deletion), and number of reads that support it. 
Note: this is entirely based on the samflag
4.	Check regions of interest in IGV (I didn’t do this for Cp but it will be the same as what we did last week for the Ct samples).
![image](https://github.com/bhazzard11/Coinfection_Analysis/assets/107216382/7d5bc88f-acd2-4aae-b893-f946b257c459)
