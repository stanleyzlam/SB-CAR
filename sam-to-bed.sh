#!/bin/sh

module load miniconda
conda activate Kazushi
tmp=/tmp 
cd /home/ks2547/project/Splinkerette/11282022_Miseq/Fastq/processed # ADJUST TO WORKING DIRECTORY

for f in *.q30.sam;
    do g="${f%.*}"
        echo "Processing $f"
        echo $g
        
        echo "Converting from SAM to BAM..."
        samtools view -b -h -S $g.sam > $g.bam
        
        echo "Converting from BAM to BED..."
        bamToBed -i $g.bam > $g.bed
        
    done


mv *.q30.bed /home/ks2547/project/Splinkerette/11282022_Miseq/Fastq/bed
cd /home/ks2547/project/Splinkerette/11282022_Miseq/Fastq/bed
for f in *.bed; do 
    mv -- "$f" "${f%.q30.bed}.bed"
done
