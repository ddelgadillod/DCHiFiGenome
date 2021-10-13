#!/bin/bash

CONDA=/home/ddelgadillo/miniconda3/bin/activate
ALLHIC=/home/ddelgadillo/Software/ALLHiC
STSPATH=/home/ddelgadillo/DC-Assembly/03_ALLHiC_Phasing
ASMPATH=/home/ddelgadillo/DC-Assembly/01_DC_Canu_Assembly/HiCanuAsm
REFPATH=/home/ddelgadillo/DC-Assembly/Data/Reference
HiCPATH=/home/ddelgadillo/DC-Assembly/Data/HiC

source $CONDA
conda activate HiC

export PATH=$ALLHiC/scripts/:$ALLHiC/bin/:$PATH

#Strategy 0: Phasing and scaffolding with RH89-039-16 v3 potato assembly as reference, Canu draft assembly v1 and raw HiC reads

mkdir -p $STSPATH/Strategy-0
cd $STSPATH/Strategy-0

draft=$ASMPATH/canu_asm.cntgs_v1.fasta
reference=$REFPATH/RH89-039-16_potato_genome_assembly.v3.fa.gz
HiC1=$HiCPATH/barrero-meneses-diacol-capiro_S3HiC_R1.fastq.gz
HiC2=$HiCPATH/barrero-meneses-diacol-capiro_S3HiC_R2.fastq.gz

echo 'Bwa index and samtools faidx to index draft genome assembly' >> ctrl_log

bwa index -a bwtsw $draft  
samtools faidx $draft  

echo 'Aligning Hi-C reads to the draft assembly' >> ctrl_log

bwa aln -t 20 $draft $HiC1 > sample_R1.sai  
bwa aln -t 20 $draft $HiC2 > sample_R2.sai  
bwa sampe $draft sample_R1.sai sample_R2.sai $HiC1 $HiC2 > sample.bwa_aln.sam  

echo 'Filtering SAM file and removing low-quality mapping hits' >> ctrl_log

PreprocessSAMs.pl sample.bwa_aln.sam $draft MBOI
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam  
samtools view -@ 20 -bt $draft.fai sample.clean.sam > sample.clean.bam  

echo 'Make Alle.cntg.table' >> ctrl.log
ref_cds=$REFPATH/RH89-039-16_potato.v3.gene_models.cds.fa
N=4
gmap_build -D . -d DB $draft
gmap -D . -d DB -t 20 -f 2 -n $N  $ref_cds > gmap.gff3
gmap2AlleleTable.pl gmap.gff3
gmap2AlleleTable.pl $REFPATH/RH89-039-16_potato_gene_models.v3.gff3

echo 'Separate homologous groups to reduce scaffolding complexity' >> ctrl_log

cp Allele.ctg.table AlleleTable/Allele.ctg.table.v2

for c in $(seq 1 12)
do 
	sed -i 's/chr$c_1/chr$c/g' Allele.ctg.table.v2
	sed -i 's/chr$c_2/chr$c/g' Allele.ctg.table.v2

done

#chrnJC.list -> Single column plain text file with chromosome names

partition_gmap.pl -g Allele.ctg.table.v2 -r $draft -b sample.clean.bam -d Partition_jc -l chrnJC.list


echo 'Build superscaffolds using ALLHiC pipeline jc' >> ctrl_log 


for c in $(seq 1 12)
do 
  cd $STSPATH/Strategy-0/Partition_jc/chr$c/
  
  cp $STSPATH/Strategy-0/sample.clean.bam .
  cp $STSPATH/Strategy-0/Allele.ctg.table.v2 .
  
  echo 'Prune chr'$c >> ctrl_log 
  ALLHiC_prune -i Allele.ctg.table.v2 -b sample.clean.bam -r seq.fasta
  
  echo 'Partition chr'$c >> ctrl_log 
  ALLHiC_partition -r seq.fasta -b prunning.bam -e GATC -k 8 -m 25
  
  echo 'Rescue chr'$c >> ctrl_log 
  ALLHiC_rescue -r seq.fasta -b sample.clean.bam -c prunning.clusters.txt -i prunning.counts_GATC.txt
  
  echo 'Optimize chr'$c >> ctrl_log 
  allhic extract sample.clean.bam seq.fasta --RE GATC
  for K in {1..8};do allhic optimize group${K}.txt sample.clean.clm;done
  
  echo 'Build chr'$c >> ctrl_log
  ALLHiC_build seq.fasta
  
  #remove large and unnecesary files
  
  rm -rf *.sam
  rm -rf sample.clean.bam
  rm -rf Allele.ctg.table.v2
  rm -rf log.txt
  rm -rf removedb_nonBest.txt
done



#Strategy 1: Phasing and scaffolding with RH89-039-16 v3 potato assembly as reference, Canu draft assembly v2(corrected with raw HiC reads) and raw HiC reads


mkdir $STSPATH/Strategy-1
cd $STSPATH/Strategy-1

draft=/home/ddelgadillo/DC-Assembly/01_DC_Canu_Assembly/canu_asm.cntgs_v1_corrv2.fasta
reference=$REFPATH/RH89-039-16_potato_genome_assembly.v3.fa.gz
HiC1=$HiCPATH/barrero-meneses-diacol-capiro_S3HiC_R1.fastq.gz
HiC2=$HiCPATH/barrero-meneses-diacol-capiro_S3HiC_R2.fastq.gz

echo 'Bwa index and samtools faidx to index draft genome assembly' >> ctrl_log

bwa index -a bwtsw $draft  
samtools faidx $draft  

echo 'Aligning Hi-C reads to the draft assembly' >> ctrl_log

bwa aln -t 20 $draft $HiC1 > sample_R1.sai  
bwa aln -t 20 $draft $HiC2 > sample_R2.sai  
bwa sampe $draft sample_R1.sai sample_R2.sai $HiC1 $HiC2 > sample.bwa_aln.sam  

echo 'Filtering SAM file and removing low-quality mapping hits' >> ctrl_log

PreprocessSAMs.pl sample.bwa_aln.sam $draft MBOI
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam  
samtools view -@ 20 -bt $draft.fai sample.clean.sam > sample.clean.bam  

echo 'Make Alle.cntg.table' >> ctrl.log
ref_cds=$REFPATH/RH89-039-16_potato.v3.gene_models.cds.fa
N=4
gmap_build -D . -d DB $draft
gmap -D . -d DB -t 20 -f 2 -n $N  $ref_cds > gmap.gff3
gmap2AlleleTable.pl gmap.gff3
gmap2AlleleTable.pl $REFPATH/RH89-039-16_potato_gene_models.v3.gff3

echo 'Separate homologous groups to reduce scaffolding complexity' >> ctrl_log

cp Allele.ctg.table AlleleTable/Allele.ctg.table.v2

for c in $(seq 1 12)
do 
	sed -i 's/chr$c_1/chr$c/g' Allele.ctg.table.v2
	sed -i 's/chr$c_2/chr$c/g' Allele.ctg.table.v2

done

#chrnJC.list -> Single column plain text file with chromosome names

partition_gmap.pl -g Allele.ctg.table.v2 -r $draft -b sample.clean.bam -d Partition_jc -l chrnJC.list


echo 'Build superscaffolds using ALLHiC pipeline jc' >> ctrl_log 


for c in $(seq 1 12)
do 
  cd $STSPATH/Strategy-1/Partition_jc/chr$c/
  
  cp $STSPATH/Strategy-1/sample.clean.bam .
  cp $STSPATH/Strategy-1/Allele.ctg.table.v2 .
  
  echo 'Prune chr'$c >> ctrl_log 
  ALLHiC_prune -i Allele.ctg.table.v2 -b sample.clean.bam -r seq.fasta
  
  echo 'Partition chr'$c >> ctrl_log 
  ALLHiC_partition -r seq.fasta -b prunning.bam -e GATC -k 8 -m 25
  
  echo 'Rescue chr'$c >> ctrl_log 
  ALLHiC_rescue -r seq.fasta -b sample.clean.bam -c prunning.clusters.txt -i prunning.counts_GATC.txt
  
  echo 'Optimize chr'$c >> ctrl_log 
  allhic extract sample.clean.bam seq.fasta --RE GATC
  for K in {1..8};do allhic optimize group${K}.txt sample.clean.clm;done
  
  echo 'Build chr'$c >> ctrl_log
  ALLHiC_build seq.fasta
  
  #remove large and unnecesary files
  
  rm -rf *.sam
  rm -rf sample.clean.bam
  rm -rf Allele.ctg.table.v2
  rm -rf log.txt
  rm -rf removedb_nonBest.txt
done



#Strategy 2: Phasing and scaffolding with RH89-039-16 v3 potato assembly as reference, Canu draft assembly v2(corrected with raw HiC reads) and HiC valid contacts

mkdir $STSPATH/Strategy-2
cd $STSPATH/Strategy-2

draft=$ASMPATH/canu_asm.cntgs_v1_corrv2.fasta
reference=$REFPATH/RH89-039-16_potato_genome_assembly.v3.fa.gz
HiC1=$HiCPATH/bm_pot_R1.fastq
HiC2=$HiCPATH/bm_pot_R1.fastq

echo 'Bwa index and samtools faidx to index draft genome assembly' >> ctrl_log

bwa index -a bwtsw $draft  
samtools faidx $draft  

echo 'Aligning Hi-C reads to the draft assembly' >> ctrl_log

bwa aln -t 20 $draft $HiC1 > sample_R1.sai  
bwa aln -t 20 $draft $HiC2 > sample_R2.sai  
bwa sampe $draft sample_R1.sai sample_R2.sai $HiC1 $HiC2 > sample.bwa_aln.sam  

echo 'Filtering SAM file and removing low-quality mapping hits' >> ctrl_log

PreprocessSAMs.pl sample.bwa_aln.sam $draft MBOI
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam  
samtools view -@ 20 -bt $draft.fai sample.clean.sam > sample.clean.bam  

echo 'Make Alle.cntg.table' >> ctrl.log
ref_cds=$REFPATH/RH89-039-16_potato.v3.gene_models.cds.fa
N=4
gmap_build -D . -d DB $draft
gmap -D . -d DB -t 20 -f 2 -n $N  $ref_cds > gmap.gff3
gmap2AlleleTable.pl gmap.gff3
gmap2AlleleTable.pl $REFPATH/RH89-039-16_potato_gene_models.v3.gff3

echo 'Separate homologous groups to reduce scaffolding complexity' >> ctrl_log

cp Allele.ctg.table AlleleTable/Allele.ctg.table.v2

for c in $(seq 1 12)
do 
	sed -i 's/chr$c_1/chr$c/g' Allele.ctg.table.v2
	sed -i 's/chr$c_2/chr$c/g' Allele.ctg.table.v2

done

#chrnJC.list -> Single column plain text file with chromosome names

partition_gmap.pl -g Allele.ctg.table.v2 -r $draft -b sample.clean.bam -d Partition_jc -l chrnJC.list


echo 'Build superscaffolds using ALLHiC pipeline jc' >> ctrl_log 


for c in $(seq 1 12)
do 
  cd $STSPATH/Strategy-2/Partition_jc/chr$c/
  
  cp $STSPATH/Strategy-2/sample.clean.bam .
  cp $STSPATH/Strategy-2/Allele.ctg.table.v2 .
  
  echo 'Prune chr'$c >> ctrl_log 
  ALLHiC_prune -i Allele.ctg.table.v2 -b sample.clean.bam -r seq.fasta
  
  echo 'Partition chr'$c >> ctrl_log 
  ALLHiC_partition -r seq.fasta -b prunning.bam -e GATC -k 8 -m 25
  
  echo 'Rescue chr'$c >> ctrl_log 
  ALLHiC_rescue -r seq.fasta -b sample.clean.bam -c prunning.clusters.txt -i prunning.counts_GATC.txt
  
  echo 'Optimize chr'$c >> ctrl_log 
  allhic extract sample.clean.bam seq.fasta --RE GATC
  for K in {1..8};do allhic optimize group${K}.txt sample.clean.clm;done
  
  echo 'Build chr'$c >> ctrl_log
  ALLHiC_build seq.fasta
  
  #remove large and unnecesary files
  
  rm -rf *.sam
  rm -rf sample.clean.bam
  rm -rf Allele.ctg.table.v2
  rm -rf log.txt
  rm -rf removedb_nonBest.txt
done



#Strategy 3: Phasing and scaffolding with RH89-039-16 v3 potato assembly as reference, Canu draft assembly v1 and HiC valid contacts

mkdir $STSPATH/Strategy-3
cd $STSPATH/Strategy-3

draft=$ASMPATH/canu_asm.cntgs_v1.fasta
reference=$REFPATH/RH89-039-16_potato_genome_assembly.v3.fa.gz
HiC1=$HiCPATH/bm_pot_R1.fastq
HiC2=$HiCPATH/bm_pot_R1.fastq

echo 'Bwa index and samtools faidx to index draft genome assembly' >> ctrl_log

bwa index -a bwtsw $draft  
samtools faidx $draft  

echo 'Aligning Hi-C reads to the draft assembly' >> ctrl_log

bwa aln -t 20 $draft $HiC1 > sample_R1.sai  
bwa aln -t 20 $draft $HiC2 > sample_R2.sai  
bwa sampe $draft sample_R1.sai sample_R2.sai $HiC1 $HiC2 > sample.bwa_aln.sam  

echo 'Filtering SAM file and removing low-quality mapping hits' >> ctrl_log

PreprocessSAMs.pl sample.bwa_aln.sam $draft MBOI
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam  
samtools view -@ 20 -bt $draft.fai sample.clean.sam > sample.clean.bam  

echo 'Make Alle.cntg.table' >> ctrl.log
ref_cds=$REFPATH/RH89-039-16_potato.v3.gene_models.cds.fa
N=4
gmap_build -D . -d DB $draft
gmap -D . -d DB -t 20 -f 2 -n $N  $ref_cds > gmap.gff3
gmap2AlleleTable.pl gmap.gff3
gmap2AlleleTable.pl $REFPATH/RH89-039-16_potato_gene_models.v3.gff3

echo 'Separate homologous groups to reduce scaffolding complexity' >> ctrl_log

cp Allele.ctg.table AlleleTable/Allele.ctg.table.v2

for c in $(seq 1 12)
do 
	sed -i 's/chr$c_1/chr$c/g' Allele.ctg.table.v2
	sed -i 's/chr$c_2/chr$c/g' Allele.ctg.table.v2

done

#chrnJC.list -> Single column plain text file with chromosome names

partition_gmap.pl -g Allele.ctg.table.v2 -r $draft -b sample.clean.bam -d Partition_jc -l chrnJC.list


echo 'Build superscaffolds using ALLHiC pipeline jc' >> ctrl_log 


for c in $(seq 1 12)
do 
  cd $STSPATH/Strategy-3/Partition_jc/chr$c/
  
  cp $STSPATH/Strategy-3/sample.clean.bam .
  cp $STSPATH/Strategy-3/Allele.ctg.table.v2 .
  
  echo 'Prune chr'$c >> ctrl_log 
  ALLHiC_prune -i Allele.ctg.table.v2 -b sample.clean.bam -r seq.fasta
  
  echo 'Partition chr'$c >> ctrl_log 
  ALLHiC_partition -r seq.fasta -b prunning.bam -e GATC -k 8 -m 25
  
  echo 'Rescue chr'$c >> ctrl_log 
  ALLHiC_rescue -r seq.fasta -b sample.clean.bam -c prunning.clusters.txt -i prunning.counts_GATC.txt
  
  echo 'Optimize chr'$c >> ctrl_log 
  allhic extract sample.clean.bam seq.fasta --RE GATC
  for K in {1..8};do allhic optimize group${K}.txt sample.clean.clm;done
  
  echo 'Build chr'$c >> ctrl_log
  ALLHiC_build seq.fasta
  
  #remove large and unnecesary files
  
  rm -rf *.sam
  rm -rf sample.clean.bam
  rm -rf Allele.ctg.table.v2
  rm -rf log.txt
  rm -rf removedb_nonBest.txt
done




#Strategy 4: Phasing and scaffolding with RH89-039-16 v3 potato assembly as reference, Canu draft assembly v3(corrected with HiC valid contacts) and HiC valid contacts

mkdir $STSPATH/Strategy-4
cd $STSPATH/Strategy-4

draft=/home/ddelgadillo/DC-Assembly/01_DC_Canu_Assembly/canu_asm.cntgs_v1_corrv3.fasta
reference=$REFPATH/RH89-039-16_potato_genome_assembly.v3.fa.gz
HiC1=$HiCPATH/bm_pot_R1.fastq
HiC2=$HiCPATH/bm_pot_R1.fastq

echo 'Bwa index and samtools faidx to index draft genome assembly' >> ctrl_log

bwa index -a bwtsw $draft  
samtools faidx $draft  

echo 'Aligning Hi-C reads to the draft assembly' >> ctrl_log

bwa aln -t 20 $draft $HiC1 > sample_R1.sai  
bwa aln -t 20 $draft $HiC2 > sample_R2.sai  
bwa sampe $draft sample_R1.sai sample_R2.sai $HiC1 $HiC2 > sample.bwa_aln.sam  

echo 'Filtering SAM file and removing low-quality mapping hits' >> ctrl_log

PreprocessSAMs.pl sample.bwa_aln.sam $draft MBOI
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam  
samtools view -@ 20 -bt $draft.fai sample.clean.sam > sample.clean.bam  

echo 'Make Alle.cntg.table' >> ctrl.log
ref_cds=$REFPATH/RH89-039-16_potato.v3.gene_models.cds.fa
N=4
gmap_build -D . -d DB $draft
gmap -D . -d DB -t 20 -f 2 -n $N  $ref_cds > gmap.gff3
gmap2AlleleTable.pl gmap.gff3
gmap2AlleleTable.pl $REFPATH/RH89-039-16_potato_gene_models.v3.gff3

echo 'Separate homologous groups to reduce scaffolding complexity' >> ctrl_log

cp Allele.ctg.table AlleleTable/Allele.ctg.table.v2

for c in $(seq 1 12)
do 
	sed -i 's/chr$c_1/chr$c/g' Allele.ctg.table.v2
	sed -i 's/chr$c_2/chr$c/g' Allele.ctg.table.v2

done

#chrnJC.list -> Single column plain text file with chromosome names

partition_gmap.pl -g Allele.ctg.table.v2 -r $draft -b sample.clean.bam -d Partition_jc -l chrnJC.list


echo 'Build superscaffolds using ALLHiC pipeline jc' >> ctrl_log 


for c in $(seq 1 12)
do 
  cd $STSPATH/Strategy-4/Partition_jc/chr$c/
  
  cp $STSPATH/Strategy-4/sample.clean.bam .
  cp $STSPATH/Strategy-4/Allele.ctg.table.v2 .
  
  echo 'Prune chr'$c >> ctrl_log 
  ALLHiC_prune -i Allele.ctg.table.v2 -b sample.clean.bam -r seq.fasta
  
  echo 'Partition chr'$c >> ctrl_log 
  ALLHiC_partition -r seq.fasta -b prunning.bam -e GATC -k 8 -m 25
  
  echo 'Rescue chr'$c >> ctrl_log 
  ALLHiC_rescue -r seq.fasta -b sample.clean.bam -c prunning.clusters.txt -i prunning.counts_GATC.txt
  
  echo 'Optimize chr'$c >> ctrl_log 
  allhic extract sample.clean.bam seq.fasta --RE GATC
  for K in {1..8};do allhic optimize group${K}.txt sample.clean.clm;done
  
  echo 'Build chr'$c >> ctrl_log
  ALLHiC_build seq.fasta
  
  #remove large and unnecesary files
  
  rm -rf *.sam
  rm -rf sample.clean.bam
  rm -rf Allele.ctg.table.v2
  rm -rf log.txt
  rm -rf removedb_nonBest.txt
done



#Strategy 5: Phasing and scaffolding with DMv6 potato assembly as reference, Canu draft assembly v1 and HiC valid contacts

mkdir $STSPATH/Strategy-5
cd $STSPATH/Strategy-5

draft=$ASMPATH/canu_asm.cntgs_v1.fasta
reference=$REFPATH/DMv6.fa.gz
HiC1=$HiCPATH/bm_pot_R1.fastq
HiC2=$HiCPATH/bm_pot_R1.fastq

echo 'Bwa index and samtools faidx to index draft genome assembly' >> ctrl_log

bwa index -a bwtsw $draft  
samtools faidx $draft  

echo 'Aligning Hi-C reads to the draft assembly' >> ctrl_log

bwa aln -t 20 $draft $HiC1 > sample_R1.sai  
bwa aln -t 20 $draft $HiC2 > sample_R2.sai  
bwa sampe $draft sample_R1.sai sample_R2.sai $HiC1 $HiC2 > sample.bwa_aln.sam  

echo 'Filtering SAM file and removing low-quality mapping hits' >> ctrl_log

PreprocessSAMs.pl sample.bwa_aln.sam $draft MBOI
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam  
samtools view -@ 20 -bt $draft.fai sample.clean.sam > sample.clean.bam  

echo 'Make Alle.cntg.table' >> ctrl.log
ref_cds=$REFPATH/DMv6.cds.fa
N=4
gmap_build -D . -d DB $draft
gmap -D . -d DB -t 20 -f 2 -n $N  $ref_cds > gmap.gff3
gmap2AlleleTable.pl gmap.gff3
gmap2AlleleTable.pl $REFPATH/DMv6.gff3

echo 'Separate homologous groups to reduce scaffolding complexity' >> ctrl_log

cp Allele.ctg.table AlleleTable/Allele.ctg.table.v2

for c in $(seq 1 12)
do 
	sed -i 's/chr$c_1/chr$c/g' Allele.ctg.table.v2
	sed -i 's/chr$c_2/chr$c/g' Allele.ctg.table.v2

done

#chrnJC.list -> Single column plain text file with chromosome names

partition_gmap.pl -g Allele.ctg.table.v2 -r $draft -b sample.clean.bam -d Partition_jc -l chrnJC.list


echo 'Build superscaffolds using ALLHiC pipeline jc' >> ctrl_log 


for c in $(seq 1 12)
do 
  cd $STSPATH/Strategy-5/Partition_jc/chr$c/
  
  cp $STSPATH/Strategy-5/sample.clean.bam .
  cp $STSPATH/Strategy-5/Allele.ctg.table.v2 .
  
  echo 'Prune chr'$c >> ctrl_log 
  ALLHiC_prune -i Allele.ctg.table.v2 -b sample.clean.bam -r seq.fasta
  
  echo 'Partition chr'$c >> ctrl_log 
  ALLHiC_partition -r seq.fasta -b prunning.bam -e GATC -k 8 -m 25
  
  echo 'Rescue chr'$c >> ctrl_log 
  ALLHiC_rescue -r seq.fasta -b sample.clean.bam -c prunning.clusters.txt -i prunning.counts_GATC.txt
  
  echo 'Optimize chr'$c >> ctrl_log 
  allhic extract sample.clean.bam seq.fasta --RE GATC
  for K in {1..8};do allhic optimize group${K}.txt sample.clean.clm;done
  
  echo 'Build chr'$c >> ctrl_log
  ALLHiC_build seq.fasta
  
  #remove large and unnecesary files
  
  rm -rf *.sam
  rm -rf sample.clean.bam
  rm -rf Allele.ctg.table.v2
  rm -rf log.txt
  rm -rf removedb_nonBest.txt
done



#Strategy 6: Phasing and scaffolding with DMv6 potato assembly as reference, Canu draft assembly v2(corrected with raw HiC reads) and HiC valid contacts

mkdir $STSPATH/Strategy-6
cd $STSPATH/Strategy-6

draft=/home/ddelgadillo/DC-Assembly/01_DC_Canu_Assembly/canu_asm.cntgs_v1_corrv2.fasta
reference=$REFPATH/DMv6.fa.gz
HiC1=$HiCPATH/bm_pot_R1.fastq
HiC2=$HiCPATH/bm_pot_R1.fastq

echo 'Bwa index and samtools faidx to index draft genome assembly' >> ctrl_log

bwa index -a bwtsw $draft  
samtools faidx $draft  

echo 'Aligning Hi-C reads to the draft assembly' >> ctrl_log

bwa aln -t 20 $draft $HiC1 > sample_R1.sai  
bwa aln -t 20 $draft $HiC2 > sample_R2.sai  
bwa sampe $draft sample_R1.sai sample_R2.sai $HiC1 $HiC2 > sample.bwa_aln.sam  

echo 'Filtering SAM file and removing low-quality mapping hits' >> ctrl_log

PreprocessSAMs.pl sample.bwa_aln.sam $draft MBOI
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam  
samtools view -@ 20 -bt $draft.fai sample.clean.sam > sample.clean.bam  

echo 'Make Alle.cntg.table' >> ctrl.log
ref_cds=$REFPATH/DMv6.cds.fa
N=4
gmap_build -D . -d DB $draft
gmap -D . -d DB -t 20 -f 2 -n $N  $ref_cds > gmap.gff3
gmap2AlleleTable.pl gmap.gff3
gmap2AlleleTable.pl $REFPATH/DMv6.gff3

echo 'Separate homologous groups to reduce scaffolding complexity' >> ctrl_log

cp Allele.ctg.table AlleleTable/Allele.ctg.table.v2

for c in $(seq 1 12)
do 
	sed -i 's/chr$c_1/chr$c/g' Allele.ctg.table.v2
	sed -i 's/chr$c_2/chr$c/g' Allele.ctg.table.v2

done

#chrnJC.list -> Single column plain text file with chromosome names

partition_gmap.pl -g Allele.ctg.table.v2 -r $draft -b sample.clean.bam -d Partition_jc -l chrnJC.list


echo 'Build superscaffolds using ALLHiC pipeline jc' >> ctrl_log 


for c in $(seq 1 12)
do 
  cd $STSPATH/Strategy-6/Partition_jc/chr$c/
  
  cp $STSPATH/Strategy-6/sample.clean.bam .
  cp $STSPATH/Strategy-6/Allele.ctg.table.v2 .
  
  echo 'Prune chr'$c >> ctrl_log 
  ALLHiC_prune -i Allele.ctg.table.v2 -b sample.clean.bam -r seq.fasta
  
  echo 'Partition chr'$c >> ctrl_log 
  ALLHiC_partition -r seq.fasta -b prunning.bam -e GATC -k 8 -m 25
  
  echo 'Rescue chr'$c >> ctrl_log 
  ALLHiC_rescue -r seq.fasta -b sample.clean.bam -c prunning.clusters.txt -i prunning.counts_GATC.txt
  
  echo 'Optimize chr'$c >> ctrl_log 
  allhic extract sample.clean.bam seq.fasta --RE GATC
  for K in {1..8};do allhic optimize group${K}.txt sample.clean.clm;done
  
  echo 'Build chr'$c >> ctrl_log
  ALLHiC_build seq.fasta
  
  #remove large and unnecesary files
  
  rm -rf *.sam
  rm -rf sample.clean.bam
  rm -rf Allele.ctg.table.v2
  rm -rf log.txt
  rm -rf removedb_nonBest.txt
done



#Strategy 7: Phasing and scaffolding with DMv6 potato assembly as reference, Canu draft assembly v3(corrected with HiC valid contacts) and HiC valid contacts

mkdir $STSPATH/Strategy-7
cd $STSPATH/Strategy-7

draft=/home/ddelgadillo/DC-Assembly/01_DC_Canu_Assembly/canu_asm.cntgs_v1_corrv3.fasta
reference=$REFPATH/DMv6_genome.fa
HiC1=$HiCPATH/bm_pot_R1.fastq
HiC2=$HiCPATH/bm_pot_R1.fastq

echo 'Bwa index and samtools faidx to index draft genome assembly' >> ctrl_log

bwa index -a bwtsw $draft  
samtools faidx $draft  

echo 'Aligning Hi-C reads to the draft assembly' >> ctrl_log

bwa aln -t 20 $draft $HiC1 > sample_R1.sai  
bwa aln -t 20 $draft $HiC2 > sample_R2.sai  
bwa sampe $draft sample_R1.sai sample_R2.sai $HiC1 $HiC2 > sample.bwa_aln.sam  

echo 'Filtering SAM file and removing low-quality mapping hits' >> ctrl_log

PreprocessSAMs.pl sample.bwa_aln.sam $draft MBOI
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam  
samtools view -@ 20 -bt $draft.fai sample.clean.sam > sample.clean.bam  

echo 'Make Alle.cntg.table' >> ctrl.log
ref_cds=$REFPATH/DMv6.cds.fa
N=4
gmap_build -D . -d DB $draft
gmap -D . -d DB -t 20 -f 2 -n $N  $ref_cds > gmap.gff3
gmap2AlleleTable.pl gmap.gff3
gmap2AlleleTable.pl $REFPATH/DMv6.gff3

echo 'Separate homologous groups to reduce scaffolding complexity' >> ctrl_log

cp Allele.ctg.table AlleleTable/Allele.ctg.table.v2

for c in $(seq 1 12)
do 
	sed -i 's/chr$c_1/chr$c/g' Allele.ctg.table.v2
	sed -i 's/chr$c_2/chr$c/g' Allele.ctg.table.v2

done

#chrnJC.list -> Single column plain text file with chromosome names

partition_gmap.pl -g Allele.ctg.table.v2 -r $draft -b sample.clean.bam -d Partition_jc -l chrnJC.list


echo 'Build superscaffolds using ALLHiC pipeline jc' >> ctrl_log 


for c in $(seq 1 12)
do 
  cd $STSPATH/Strategy-7/Partition_jc/chr$c/
  
  cp $STSPATH/Strategy-7/sample.clean.bam .
  cp $STSPATH/Strategy-7/Allele.ctg.table.v2 .
  
  echo 'Prune chr'$c >> ctrl_log 
  ALLHiC_prune -i Allele.ctg.table.v2 -b sample.clean.bam -r seq.fasta
  
  echo 'Partition chr'$c >> ctrl_log 
  ALLHiC_partition -r seq.fasta -b prunning.bam -e GATC -k 8 -m 25
  
  echo 'Rescue chr'$c >> ctrl_log 
  ALLHiC_rescue -r seq.fasta -b sample.clean.bam -c prunning.clusters.txt -i prunning.counts_GATC.txt
  
  echo 'Optimize chr'$c >> ctrl_log 
  allhic extract sample.clean.bam seq.fasta --RE GATC
  for K in {1..8};do allhic optimize group${K}.txt sample.clean.clm;done
  
  echo 'Build chr'$c >> ctrl_log
  ALLHiC_build seq.fasta
  
  #remove large and unnecesary files
  
  rm -rf *.sam
  rm -rf sample.clean.bam
  rm -rf Allele.ctg.table.v2
  rm -rf log.txt
  rm -rf removedb_nonBest.txt
done
