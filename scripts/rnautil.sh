#!/bin/bash

if [ "$1" == "align.tx" ]; then
    minimap2 -ax map-ont /dilithium/Data/Nanopore/Analysis/RNAconsort/180319_rnaconsort/mettl3/EEF.refseq.fa /dilithium/Data/Nanopore/Analysis/RNAconsort/NA12878-DirectRNA.pass.dedup.fastq.gz | samtools sort -o dRNA.allEEF.sorted.bam -T reads.tmp
    samtools index dRNA.allEEF.sorted.bam
    minimap2 -ax map-ont /dilithium/Data/Nanopore/Analysis/RNAconsort/180319_rnaconsort/mettl3/EEF.refseq.fa /dilithium/Data/Nanopore/Analysis/RNAconsort/IVT/01_29_18_R94_12878IVT_pass.fastq | samtools sort -o IVT.allEEF.sorted.bam -T reads.tmp
    samtools index IVT.allEEF.sorted.bam
fi

if [ "$1" == "align.wg" ]; then
    minimap2 -ax splice -uf -k14 /mithril/Data/NGS/Reference/human38/GRCH38.fa /dilithium/Data/Nanopore/rna/IVT/01_29_18_R94_12878IVT_pass.fastq | samtools sort -o IVT.GRCh38.sorted.bam -T reads.tmp
    samtools index IVT.GRCh38.sorted.bam
    minimap2 -ax splice -uf -k14 /mithril/Data/NGS/Reference/human38/GRCH38.fa /dilithium/Data/Nanopore/rna/fastq/NA12878-DirectRNA.pass.dedup.fastq.gz | samtools sort -o dRNA.GRCh38.sorted.bam -T reads.tmp
    samtools index dRNA.GRCh38.sorted.bam
fi

if [ "$1" == "eventalign" ]; then
    source activate polish
    ~/nanopolish/nanopolish eventalign --scale-events -n -t 8 --reads /shared/directRNA/NA12878-DirectRNA.pass.dedup.fastq.gz --bam dRNA.GRCh38.sorted.bam --genome GRCh38.fasta >dRNA.EEF2.hg38.scaled.txt
    ~/nanopolish/nanopolish eventalign --scale-events -n -t 8 --reads /shared/IVT/IVT_NA12878.fastq --bam IVT.GRCh38.sorted.bam --genome GRCh38.fasta >IVT.EEF2.hg38.scaled.txt 
fi

if [ "$1" == "pare" ]; then
    awk 'NR == 1; NR > 1 {print $0 | "grep -E 'GGACT'"}' dRNA.EEF2.hg38.scaled.txt | sed 's/T/U/g' > dRNA.EEF2.GGACU.evalign.txt
    awk 'NR == 1; NR > 1 {print $0 | "grep -E 'GGACT'"}' IVT.EEF2.hg38.scaled.txt | sed 's/T/U/g' > IVT.EEF2.GGACU.evalign.txt
fi







    





