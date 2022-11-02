#!/bin/bash
##step1.filter low quality sequence by fastp(v0.20.0)
fastp -i ${sampleID}.R1.fastq.gz -o $od/QC_filter/$sam.clean_R1.fastq.gz -I ${sampleID}.R2.fastq.gz -O $od/QC_filter/$sam.clean_R2.fastq.gz \
-h $od/QC_filter/fastp/$sam.filt.html -j $od/QC_filter/fastp/$sam.filt.json --cut_by_quality3 -W 4 -M 20 -n 5 -c -l 50 -w 3

##step2.filter host sequencing (bwa v2.2.1)
#index of reference genome
bwa-mem2.avx2 index $dbs{host}{genome_filt}
#mapping to reference and extract mapped sequence
bwa-mem2.avx2 mem -t 16 -T 30 $dbs{host}{genome_filt} $od/QC_filter/$sam.clean_R1.fastq.gz $od/QC_filter/$sam.clean_R2.fastq.gz \
-R \'\@RG\\tID:$sam_RG\\tPL:illumina\\tSM:$sam_RG\' | samtools view -S -b -f 4 - | \
bedtools bamtofastq -i /dev/stdin -fq $od/filt_host/$sam.filt_host_R1.fastq -fq2 $od/filt_host/$sam.filt_host_R2.fastq 
#Generate quality reports
fastp -A -G -Q -L -i $od/filt_host/$sam.filt_host_R1.fastq -o $od/filt_host/$sam.final_R1.fastq.gz -I $od/filt_host/$sam.filt_host_R2.fastq -O $od/filt_host/$sam.final_R
2.fastq.gz -h $od/filt_host/fastp/$sam.final.html -j $od/filt_host/fastp/$sam.final.json
