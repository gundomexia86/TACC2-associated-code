#!/bin/bash
## Step 1: star alignment
export SAMPLE=B109-NC_L2_375375
export DATA=/data/xiatl/RNA-seq-LZR/cleandata/
export RESULT=/data/xiatl/RNA-seq-LZR/result/${SAMPLE}
mkdir ${RESULT}

gunzip ${DATA}/${SAMPLE}.R1.clean.fastq.gz
gunzip ${DATA}/${SAMPLE}.R2.clean.fastq.gz
STAR --genomeDir /data/data_genome/human_hg38/ --readFilesIn ${DATA}/${SAMPLE}.R1.clean.fastq \
${DATA}/${SAMPLE}.R2.clean.fastq \
--runThreadN 4 --outSAMtype SAM --outSAMstrandField intronMotif \
--outFileNamePrefix ${RESULT}/${SAMPLE}_
samtools view -bS ${RESULT}/${SAMPLE}_Aligned.out.sam >${RESULT}/${SAMPLE}_Aligned.out.bam
samtools sort ${RESULT}/${SAMPLE}_Aligned.out.bam -o ${RESULT}/${SAMPLE}_Aligned_sort.bam
gzip ${DATA}/${SAMPLE}.R1.clean.fastq
gzip ${DATA}/${SAMPLE}.R2.clean.fastq

## Step 2: cuffquant
#!/bin/bash
export SAMPLE=B109-NC_L2_375375
export RESULT=/data/xiatl/RNA-seq-LZR/result/${SAMPLE}

cuffquant -p 4 --no-update-check -o ${RESULT}/${SAMPLE}_quant \
/data/data_genome/human_hg38/human.gtf \
--library-type ff-unstranded \
${RESULT}/${SAMPLE}_Aligned_sort.bam

## Step 3: cuffnorm
#!/bin/bash
export RESULT=/data/xiatl/RNA-seq-LZR/result/
mkdir ${RESULT}/Norm

cuffnorm -p 4 -L B109-NC,B109-T1,B109-T2,B150-NC,B150-T1,B150-T2,EC18-NC,EC18-T1,EC18-T2,TE3-NC,TE3-T1,TE3-T2 \
--no-update-check /data/data_genome/human_hg38/human.gtf \
-o ${RESULT}/Norm \
${RESULT}/B109-NC_L2_375375/B109-NC_L2_375375_quant/abundances.cxb \
${RESULT}/B109-T1_L2_376376/B109-T1_L2_376376_quant/abundances.cxb \
${RESULT}/B109-T4_L2_377377/B109-T2_L2_377377_quant/abundances.cxb \
${RESULT}/B150-NC_L2_369369/B150-NC_L2_369369_quant/abundances.cxb \
${RESULT}/B150-T1_L2_370370/B150-T1_L2_370370_quant/abundances.cxb \
${RESULT}/B150-T4_L2_371371/B150-T2_L2_371371_quant/abundances.cxb \
${RESULT}/EC18-NC_L2_365365/EC18-NC_L2_365365_quant/abundances.cxb \
${RESULT}/EC18-T1_L2_366366/EC18-T1_L2_366366_quant/abundances.cxb \
${RESULT}/EC18-T4_L2_368368/EC18-T2_L2_368368_quant/abundances.cxb \
${RESULT}/TE3-NC_L2_372372/TE3-NC_L2_372372_quant/abundances.cxb \
${RESULT}/TE3-T1_L2_373373/TE3-T1_L2_373373_quant/abundances.cxb \
${RESULT}/TE3-T4_L2_374374/TE3-T2_L2_374374_quant/abundances.cxb

