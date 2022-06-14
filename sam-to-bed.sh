#!/bin/sh

p=6 # number of threads/cores
tmp=/tmp 
cd /gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/Fastq/processed

for f in *.q30.sam;
    do g="${f%.*}"
        echo "Processing $f"
        echo $g
        
        echo "Converting from SAM to BAM..."
        samtools view -b -h -S $g.sam > $g.bam
        
        echo "Converting from BAM to BED..."
        bamToBed -i $g.bam > $g.bed
        
    done


mv *.q30.bed /gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/Fastq/bed
cd /gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/Fastq/bed
for f in *.bed; do 
    mv -- "$f" "${f%.q30.bed}.bed"
done