# SB-CAR

Code for insertion site mapping analysis for Ye et al. 2022 (In-Revision at Nature Biomedical Engineering)


Prerequisites:
1. Packages: cutadapt, hisat2, samtools, bedtools
2. Need BBDuk.sh script for quality filtering
3. Build GRCh38 human genome assembly for raw-to-sam.sh
4. R Packages as listed in splink_visualization.R
5. Comparative figures (e.g. MAJESTIC vs LV, MAJESTIC vs MCMC) require reference .bed files that are included in this repository. They were obtained either from Supplementary Data or direct request from two publications: Miskey et al. 2022, Nucleic Acids Research https://academic.oup.com/nar/article/50/5/2807/6533617#340619945 and Querques et al. 2019, Nature Biotechnology https://www.nature.com/articles/s41587-019-0291-z. Data from Querques et al were aligned against the hg19 genome and were thus lifted to the hg38 assembly using the web LiftOver tool at https://genome.ucsc.edu/cgi-bin/hgLiftOver (Minimum ratio of bases that overlap = 0.95).



Quick Workflow: 
1. Use raw-to-sam.sh to convert raw .fastq files to .sam format. Requires bbduk.sh
2. Use sam-to-bed.sh to convert .sam files to .bed format. 
3. Use splink_visualization.R to convert .bed files to GRanges objects, which can then be used to generate a variety of visualizations including karyograms and figures showing frequency of genomic safe harbor and functional gene region insertions. 




