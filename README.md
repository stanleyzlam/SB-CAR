# SB-CAR

Code for insertion site mapping analysis for Ye et al. 2023 (Nature Biomedical Engineering; https://www.nature.com/articles/s41551-023-01058-6). 
The main goal with this analysis was to map the genomic insertion site profile of our new gene-transfer method MAJESTIC and to compare its insertion profile with that of the Minicircle/mRNA gene transfer method. We sequenced three cell samples side-by-side, with three technical replicates each: 1) AAV+mRNA (CAR-Ts generated via Majestic), 2) MC+mRNA (CAR-Ts generated via Minicircle/mRNA electroporation), and 3) AAVonly (control T cells). 


## Quick Workflow: 
1. Use raw-to-sam.sh to convert raw .fastq files to .sam format. Requires bbduk.sh from BBTools
2. Use sam-to-bed.sh to convert .sam files to .bed format. 
3. Use splink_visualization.R to convert .bed files to GRanges objects, which can then be used to generate a variety of visualizations including karyograms and figures showing frequency of genomic safe harbor and functional gene region insertions. 


## Prerequisites:
1. Packages: cutadapt, hisat2, samtools, bedtools
2. Need BBDuk.sh script for quality filtering
3. Build GRCh38 human genome assembly for mapping reads to the genome (needed for raw-to-sam.sh)
4. R Packages as listed in splink_visualization.R
5. Comparative figures (e.g. MAJESTIC vs MC) require reference .bed files that are included in this repository. They were obtained either from Supplementary Data or direct request from two publications: Miskey et al. 2022, Nucleic Acids Research https://academic.oup.com/nar/article/50/5/2807/6533617#340619945 and Querques et al. 2019, Nature Biotechnology https://www.nature.com/articles/s41587-019-0291-z. We thank Csaba Miskey for sharing data from Querques et al, 2019. Data from Querques et al were originally aligned against the hg19 genome and were thus lifted to the hg38 assembly using the web LiftOver tool at https://genome.ucsc.edu/cgi-bin/hgLiftOver (Minimum ratio of bases that overlap = 0.95).


## Set up directories for data processing via raw-to-sam.sh and sam-to-bed.sh

/home/ks2547/project/Splinkerette/grch38/genome --> contains indexed GRCh38 human genome assembly **Please change this directory to where you keep your indexed genome files in your own run-through**

Root directory: "/home/ks2547/project/Splinkerette/11282022_Miseq/Fastq"

(Root) --> where .fastq files are stored

(Root)/processed --> contains all intermediary files during .fastq to .bed processing

(Root)/bed --> contains all output .bed files


## Set up directories for data visualization via splink_visualization.R
Root directory: "/Users/kazushisuzuki/Desktop/MAJESTIC"

(Root)/20221128_Miseq/Fastq/bed --> where sample .bed files are stored

(Root)/hssb_files --> contains comparison .bed files (safe harbors, lentiviral insertion profile, etc.) from Querques et al., 2019

(Root)/random --> contains a .bed file with 1e6 randomly generated sequences

(Root)/functional_gene_regions/merged --> contains functionally annotated gene region .bed files 

(Root)/figures/karyogram/unique --> output for karyogram plots

(Root)/figures/safeharbors/unique --> output for safe harbor plots

(Root)/figures/functionalheatmap/unique --> output for functionally annotated regions plots


## Structure of Visualization R script: 

### I. Libraries
All libraries used in this script can be loaded here. 

### II. Load datasets into GRanges
Here, we take the processed .bed files for the three conditions (AAV+mRNA, MC+mRNA, AAVonly) and we convert them into GRanges objects (https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html). 
Additionally, for comparisons, we also load the following objects into GRanges objects:
1. A set of 1 million randomly mapped sequences generated via the bedtools random function 
2. A set of .bed files provided by Csaba Miskey (Querques et al., 2019). These include a custom compiled genomic safe harbor .bed file and a lentivirus insertion .bed file (data originally from Roth et al., 2018).
An important variable is "keep_redundant". We set "keep_redundant=F" to indicate that for each .bed file, only unique insertion sites should be kept. 

### III. Karyograms
Here, we visualize the GRanges objects as karyograms. We concatenate all of the separate GRanges objects into a list, and then run a loop to generate karyograms for all of the .bed files. 

### IV. Insertion into genomic safe harbors
Here, we determine the extent of overlap between each condition's insertion profile and a set of genomic safe harbors compiled by Querques et al., 2019. To do this, we use the countOverlaps function from the GenomicRanges package. We create a data frame from these ratios, and plot the results as a horizontal bar chart. 

### V. Insertion into functionally annotated regions
Lastly, we determine the extent of overlap between each condition's insertion profile and functional gene regions .bed files, which we obtained from UCSC Table Browser (Gencode v39). We use the countOverlaps function, create a data frame from these ratios, and plot the results as a pheatmap. 



