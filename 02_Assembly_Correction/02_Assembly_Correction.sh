#!/bin/bash


HICDATAPATH=/home/ddelgadillo/DC-Assembly/Data/HiC
ASMPATH=/home/ddelgadillo/DC-Assembly/01_DC_Canu_Assembly/HiCanuAsm

ln -s $ASMPATH/canu_asm.cntgs_v1.fasta
draft=canu_asm.cntgs_v1.fasta

bwa index -a bwtsw canu_asm.cntgs_v1.fasta  

#Canu draft assembly v2 corrected with raw HiC reads


HiC1=$HICDATAPATH/barrero-meneses-diacol-capiro_S3HiC_R1.fastq.gz
HiC2=$HICDATAPATH/barrero-meneses-diacol-capiro_S3HiC_R2.fastq.gz

bwa mem canu_asm.cntgs_v1.fasta $HiC1 $HiC2 -t 20 -k 20 > alCanuv1PEHiC.sam

samtools sort alCanuv1PEHiC.sam -@ 20 -o alCanuv1PEHiC.bam
samtools index -@ 20 alCanuv1PEHiC.bam
ALLHiC_corrector -m alCanuv1PEHiC.bam -r canu_asm.cntgs_v1.fasta -o canu_asm.cntgs_v1_corrv2.fasta -t 20

cp canu_asm.cntgs_v1_corrv2.fasta $ASMPATH


#HiC 


#Canu draft assembly v3 corrected with HiC valid contacts

HiC1Valid=$HICDATAPATH/bm_pot_R1.fastq
HiC2Valid=$HICDATAPATH/bm_pot_R1.fastq


bwa mem canu_asm.cntgs_v1.fasta $HiC1Valid $HiC2Valid -t 20 -k 20 > alCanuv1PEHiCValid.sam

samtools sort alCanuv1PEHiCValid.sam -@ 20 -o alCanuv1PEHiCValid.bam
samtools index -@ 20 alCanuv1PEHiCValid.bam
ALLHiC_corrector -m alCanuv1PEHiCValid.bam -r canu_asm.cntgs_v1.fasta -o canu_asm.cntgs_v1_corrv3.fasta -t 20

cp canu_asm.cntgs_v1_corrv3.fasta $ASMPATH
