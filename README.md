# Coinfection_Analysis
Scripts used to analyze experimental coinfection of two distantly related strains of the same species.

For Single Cell Analyis:

Generate list of snps for each genome from mapped WGS data: 
1. Run mpileup on each sorted bam file (samtools mpileup -f path/to/Genome.fasta path/to/sorted.bam > path/to/output.mpileup)
2. Run find_polymorphisms.pl to find the location of all differences between your sample and the reference (Default is >50 reads per position, with at least 90% difference             from reference). Outputs a list of snps (Out_snp_list.txt), list of areas with increased variability (OUT_mask_list.txt), and ratio of polyclonality at a given location            (OUT_polyclonal_ratio.txt).
3. Run make_site_list.pl to generate list of unique snps between the two strains.
4. Find deletions between two strains: Run find_cross_deletions.pl (perl find_cross_deletions.pl 1000 path/to/sample1.mpileup path/to/sample2.mpileup snp_comparison_OUT.txt). The output is a txt file with the: name of the window by chromosome_window, size of the window you inputed, Input1 number of reads, Input2 number of reads.
5. For each single cell sample: Run mpileup on each sorted bam file using the site list generated from WGS data (samtools mpileup -l site_list.txt -f path/to/Genome.fasta --output-extra QNAME    path/to/sorted.bam > path/to/output.mpileup)
6. Run genotype_individual_cells.pl, takes input of unique snps of each genotype, a list of cells for each sample (needs to match read name in mpilup), and mpileup, and     outputs .txt file with cell name, the number of reads at each location that corrispond to each genotype (cell must have at least 5 reads at each location).
7. Run make_per_cell_snp_map.pl to generate by chromosome_postion, map of each cell. Outputs txt file with allele at each unique snp position in the genome (if no reads there outputs NA).
8. Generate fasta files of snps for recombinate cells using make_recomb_cell_fasta.pl (Uses only snp positions with at least 200 cells genotyped at that postion, file      generated by recomb_cells_per_pos.pl).
9. Make heatmaps of per cell recombination patterns (compressed over 200,000bp windows) with make_multi_cell_window_maps.pl.

For Bulk RNA-seq analysis (Cryptosporidium):
Overview of what I did and what each code does:
1.	Run mpileup on each sorted bam file (samtools mpileup -f path/to/Genome.fasta path/to/sorted.bam > path/to/output.mpileup)
2.	Ran find_cross_deletions.pl (perl find_cross_deletions.pl 1000 mapped_Cp_isra_RG.mpileup mapped_Cp_tom_RG.mpileup snp_comparison_Cp.txt). The output is a txt file with the: name of the window by chromosome_window, size of the window (I just ran at 1000bp but this can be changed), Input1 number of reads, Input2 number of reads.
3.	Ran paired_end_analysis (perl paired_end_analysis.pl mapped_Cp_tom_RG.bam mapped_Cp_tom_RG_paired_end_analysis.txt). First input is the bam file (it doesn’t need to be sorted) then output. File output is location, what type it is (Insertion, Inversion, or Deletion), and number of reads that support it. 
Note: this is entirely based on the samflag
4.	Check regions of interest in IGV.
