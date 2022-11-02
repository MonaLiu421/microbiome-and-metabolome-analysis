#!/bin/bash
#reads assemble to cpntig (Megahit v1.1.3)
megahit --min-count 2 --k-min 27 --k-max 127 --k-step 20 --num-cpu-threads 32 --min-contig-len 1000 \
-1 /final_reads/$sam.final_R1.fastq.gz \
-2 /final_reads/$sam.final_R2.fastq.gz \
-o /assembl/sample_assembl/$sam \

bowtie2-build /assembl/sample_assembl/$sam/$sam.contigs.fasta /assembl/sample_assembl/$sam/$sam.contigs.fasta

bowtie2 -p 32 --no-unal  -x /assembl/sample_assembl/$sam/$sam.contigs.fasta -1 /final_reads/final_reads/$sam.final_R1.fastq.gz 
-2 /final_reads/final_reads/$sam.final_R2.fastq.gz --un-conc-gz /assembl/unmapped_data/$sam.unmapped_contig.gz 2> /assembl/sample_assembl/$sam/$sam.bowtie2.log 
