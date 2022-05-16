#!/bin/bash

## Modern data processing and mapping using mem

## $1 - Threads
## $2 - Raw reads file
## $3 - Code name for sample
## $4 - Results folder
## $5 - Reference
## $6 - Mitogenome reference

## merge sequencer output fastq files together
#zcat $2/*$3*R1*.fastq.gz > $2/$3_R1.fastq
#zcat $2/*$3*R2*.fastq.gz > $2/$3_R2.fastq

### R1 and R2 adapter trimmer and remove reads shorter than 30
#pigz -d -p $1 $2/$3_R1.fastq* $2/$3_R2.fastq*
#~/Software/skewer/skewer -m pe -l 30 -t $1 $2/$3_1.fastq.gz $2/$3_2.fastq.gz -o $4/$3

## map merged and non merged reads separately
#bwa index $5
bwa mem -M -t $1 $5 $4/$3-trimmed-pair1.fastq $4/$3-trimmed-pair2.fastq > $4/$3.sam
samtools view -F 4 -q 30 -uS -@ $1 $4/$3.sam -o $4/$3.bam
samtools sort -@ $1 -m 2G $4/$3.bam -o $4/$3_map.sort.bam
rm $4/$3.sam $4/$3.bam

## Remove duplicates
samtools rmdup -S $4/$3_map.sort.bam $4/$3_map.rmdup.bam

## Sort, index and produce mapping results
samtools sort -m 2G -@ $1 $4/$3_map.rmdup.bam -o $4/$3_map.rmdup.sort.bam
samtools index $4/$3_map.rmdup.sort.bam
samtools flagstat $4/$3_map.rmdup.sort.bam > $4/$3.mapping.results.txt
samtools depth $4/$3_map.rmdup.sort.bam | awk '{sum+=$3;cnt++}END{print " read depth " sum/cnt " total mapped bp " sum}' >> $4/$3.mapping.results.txt

## realign around indels
# prepare reference (do in advance)
#picard CreateSequenceDictionary REFERENCE=Bowhead_scaffolds.fasta OUTPUT=Bowhead_scaffolds.dict

## Add read groups
picard AddOrReplaceReadGroups I=$4/$3_map.rmdup.sort.bam O=$4/$3_map.rmdup.sort_RG.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$3
samtools index $4/$3_map.rmdup.sort_RG.bam
## realign around indels using GATK
#GenomeAnalysisTK -T RealignerTargetCreator -R $5 -I $4/$3_map.rmdup.sort_RG.bam -o $4/$3.mapped.intervals -nt $1
#GenomeAnalysisTK -T IndelRealigner -R $5 -I $4/$3_map.rmdup.sort_RG.bam -targetIntervals $4/$3.mapped.intervals -o $4/$3_map.rmdup.sort_RG_realigned.bam

### Map and parse for mitogenome (same as above)
#bwa mem -M -t $1 $6 $4/$3-trimmed-pair1.fastq $4/$3-trimmed-pair2.fastq | samtools view -F 4 -q 30 -uS -@ $1 - | samtools sort -@ $1 -o $4/$3_map.mito.sort.bam
#samtools rmdup -S $4/$3_map.mito.sort.bam $4/$3_map.mito.rmdup.bam
#samtools sort -@ $1 $4/$3_map.mito.rmdup.bam -o $4/$3_map.mito.rmdup.sort.bam
#samtools index $4/$3_map.mito.rmdup.sort.bam
#picard AddOrReplaceReadGroups I=$4/$3_map.mito.rmdup.sort.bam O=$4/$3_map.mito.rmdup.sort_RG.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$3
#samtools index $4/$3_map.mito.rmdup.sort_RG.bam


#samtools mpileup -Q 30 -uf $5 $4/$3_map.rmdup.sort.bam | bcftools call -c - | vcfutils.pl vcf2fq -d 10 | gzip > $4/$3_diploid.fq.gz
pigz -p $1 $2/$3*.fastq
rm $4/$3_map.sort.bam $4/$3_map.rmdup.bam $4/$3_map.mito.sort.bam $4/$3_map.mito.rmdup.sort.bam $4/$3_map.rmdup.sort.bam

