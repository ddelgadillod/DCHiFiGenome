#!/bin/bash

#Flye Assembly

flye --pacbio-hifi /home/ddelgadillo/DC-Assembly/Data/HiFiReads/m64041_hifi_reads.fastq.gz --out-dir Flye_hfPotAssemFlye --genome-size 2.8g --threads 20 --keep-haplotypes --resume-from consensus


#Canu Assembly

canu -p asm -d CanuhfPotAssem_hptyps genomeSize=700m -pacbio-hifi /home/ddelgadillo/DC-Assembly/Data/HiFiReads/m64041_hifi_reads.fastq.gz useGrid=false corOutCoverage=200 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" 


#Hifi ASM Assembly

hifiasm -o HASMhfPotAssem -t20 -D10 /home/ddelgadillo/DC-Assembly/Data/HiFiReads/m64041_hifi_reads.fastq.gz --high-het


#IPA Assembly

ipa local --nthreads 20 --njobs 1 -i /home/ddelgadillo/DC-Assembly/Data/HiFiReads/m64041_hifi_reads.fastq.gz