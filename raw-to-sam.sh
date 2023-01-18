#!/bin/sh

module load miniconda
conda activate Kazushi


cd /home/ks2547/project/Splinkerette/11282022_Miseq/Fastq # ADJUST TO WORKING DIRECTORY
## Please also adjust directory for GRCh38 below: /home/ks2547/project/Splinkerette/grch38/genome

#gunzip *.gz

echo "Removing existing splink_log.txt (count log file), if present"
rm -f -- splink_log.txt ##reset file counts log

echo "sample,fastq,fastq_q26,nonmobilized_removed,SBadaptor_present,mappedreads_q30" >> splink_log.txt

for f in *.fastq;
    do g="${f%.*}"
        echo "Processing $f"
        
        step0=$(wc -l ${g}.fastq)
        reads=$(echo $step0 | cut -d" " -f1) #extract line count
        reads=$((reads/4)) #calculate read count from line count
        step0="${reads}"
        
        echo "Performing initial QC with BBDuk..."
        bbduk.sh in=$g.fastq out=processed/${g}.q26.fastq trimq=26 minlen=80 maq=30 qtrim=rl #QC 
        step1=$(wc -l processed/${g}.q26.fastq)
        reads=$(echo $step1 | cut -d" " -f1) #extract line count
        reads=$((reads/4)) #calculate read count from line count
        step1="${reads}"

        echo "Removing non-integrated AAV using seqs from AAV ITR..."
        cutadapt --report=minimal -g TATAGTCTAGAACGCGTGCG -e 0.1 --overlap 15 --discard-trimmed processed/${g}.q26.fastq > processed/${g}.mobilized.q26.fastq #filter OUT transposon sequences that were not mobilized/inserted in genome (based on position of ITR)
        step2=$(wc -l processed/${g}.mobilized.q26.fastq)
        filename=$(echo $step2 | cut -d" " -f2) #extract filename
        reads=$(echo $step2 | cut -d" " -f1) #extract line count
        reads=$((reads/4)) #calculate read count from line count
        step2="${reads}"
            
        echo "Trimming out the transposon arms and keeping ONLY those seqs that were trimmed"
        cutadapt --report=minimal -g ^ACTTCAACTG -e 0.1 -m 15 --overlap 10 --discard-untrimmed processed/${g}.mobilized.q26.fastq > processed/${g}.kept.mobilized.q26.fastq # trim out transposon arm + keep ONLY seqs with transposon arm
        step3=$(wc -l processed/${g}.kept.mobilized.q26.fastq)
        reads=$(echo $step3 | cut -d" " -f1) #extract line count
        reads=$((reads/4)) #calculate read count from line count
        step3="${reads}"
            
        echo "Removing some bases from the 5p end" 
        cutadapt -u 10 -o processed/${g}.cut5p.kept.mobilized.q26.fastq processed/${g}.kept.mobilized.q26.fastq
        echo "Trimming reads down to fixed length"
        cutadapt -l 30 -o processed/${g}.cut3p.cut5p.kept.mobilized.q26.fastq processed/${g}.cut5p.kept.mobilized.q26.fastq
        echo "Mapping reads using HISAT2..."
        hisat2 -x /home/ks2547/project/Splinkerette/grch38/genome -U processed/${g}.cut3p.cut5p.kept.mobilized.q26.fastq -S processed/${g}.sam

        echo "Filtering out mapped reads with quality less than 30..."
        samtools view -q 30 -h processed/${g}.sam > processed/${g}.q30.sam
        step4=$(wc -l processed/${g}.q30.sam)
        reads=$(echo $step4 | cut -d" " -f1) #extract line count
        step4="${reads}"
        
        echo "$g,$step0,$step1,$step2,$step3,$step4" >> splink_log.txt

        echo "File line counts saved to splink_log.txt"
        
    done
    
