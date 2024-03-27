#!/bin/bash

##QC SIA sequencing data and query against viral database. Parse blastx outputs in MEGAN to get read counts for viral taxa. Repeat for each sample.##
gunzip Sample1_R1_001.fastq.gz;
gunzip Sample1_R2_001.fastq.gz;
cutadapt -a file:illumina_dna_prep_adapters.fa -A file:illumina_dna_prep_adapters.fa -o Sample1_cutadapt_adaptTrim_R1_2022.fastq -p Sample1_cutadapt_adaptTrim_R2_2022.fastq Sample1_R1_001.fastq Sample1_R2_001.fastq -j 28 1> Sample1_cutadapt_adaptTrim_2022.log.txt 2>&1;
cutadapt -a GTTTCCCAGTCACGATC -a GATCGTGACTGGGAAAC -A GTTTCCCAGTCACGATC -A GATCGTGACTGGGAAAC -o Sample1_cutadapt_primerTrimR_R1_2022.fastq -p Sample1_cutadapt_primerTrimR_R2_2022.fastq Sample1_cutadapt_adaptTrim_R1_2022.fastq Sample1_cutadapt_adaptTrim_R2_2022.fastq -j 28 1> Sample1_cutadapt_primerTrimR_2022.log.txt 2>&1;
cutadapt -g GTTTCCCAGTCACGATC -g GATCGTGACTGGGAAAC -G GTTTCCCAGTCACGATC -G GATCGTGACTGGGAAAC -o Sample1_cutadapt_primerTrimRL_R1_2022.fastq -p Sample1_cutadapt_primerTrimRL_R2_2022.fastq Sample1_cutadapt_primerTrimR_R1_2022.fastq Sample1_cutadapt_primerTrimR_R2_2022.fastq -j 28 1> Sample1_cutadapt_primerTrimRL_2022.log.txt 2>&1;
bbduk.sh in=Sample1_cutadapt_primerTrimRL_R1_2022.fastq in2=Sample1_cutadapt_primerTrimRL_R2_2022.fastq out=Sample1-cutadapt_QC_R1_2022.fastq out2=Sample1-cutadapt_QC_R2_2022.fastq qtrim=rl trimq=30 minlength=75 minavgquality=20 removeifeitherbad=f otm=t tpe=t overwrite=t 1> Sample1-cutadapt_QC_2022.log.txt 2>&1;
bbduk.sh in=Sample1-cutadapt_QC_R1_2022.fastq in2=Sample1-cutadapt_QC_R2_2022.fastq ref=phix174_ill.ref.fa.gz out=Sample1-R1-cutadapt_phixRemoved_2022.fastq out2=Sample1-R2-cutadapt_phixRemoved_2022.fastq k=31 hdist=1 overwrite=t 1> Sample1_cutadapt_phixRemoved_2022.log.txt 2>&1;
bbmap.sh minid=.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 -Xmx64g path=Human_GRCh38 in=Sample1-R1-cutadapt_phixRemoved_2022.fastq in2=Sample1-R2-cutadapt_phixRemoved_2022.fastq outu=Sample1_cutadapt_hostRemoved_2022.fastq outm=Sample1_cutadapt_hostMatched_2022.fastq 1>Sample1_cutadapt_hostRemoval.log_2022.txt 2>&1;
dedupe.sh in=Sample1_cutadapt_hostRemoved_2022.fastq out=Sample1_cutadapt_firstDeduplication_2022.fastq outd=Sample1_cutadapt_firstDuplication_2022.fastq csf=dedupe.cluster.stats overwrite=t minidentity=99 1> Sample1_cutadapt_firstDeduplication.log_2022.txt 2>&1;
bbmerge.sh in=Sample1_cutadapt_firstDeduplication_2022.fastq out=Sample1_cutadapt_firstDeduplicationMerged_2022.fastq outu=Sample1_cutadapt_firstDeduplicationUnMerged_2022.fastq 1>Sample1_cutadapt_firstDeduplicationMerged.log_2022.txt 2>&1;
cat Sample1_cutadapt_firstDeduplicationMerged_2022.fastq Sample1_cutadapt_firstDeduplicationUnMerged_2022.fastq > Sample1_cutadapt_firstDeduplicationMerged_UnMerged_2022.fastq;
dedupe.sh in=Sample1_cutadapt_firstDeduplicationMerged_UnMerged_2022.fastq out=Sample1_cutadapt_secondDeduplication_2022.fastq outd=Sample1_cutadapt_secondDuplication_2022.fastq csf=dedupe.cluster.stats overwrite=t minidentity=100 ac=f 1> Sample1_cutadapt_secondDeduplication.log_2022.txt 2>&1;
bbduk.sh in=Sample1_cutadapt_secondDeduplication_2022.fastq out=Sample1_cutadapt_secondDeduplication_filtered_2022.fastq minlength=75 overwrite=t 1>Sample1_cutadapt_secondDeduplication_filtered.log_2022.txt 2>&1; 
sed -n '1~4s/^@/>/p;2~4p' Sample1_cutadapt_secondDeduplication_filtered_2022.fastq > Sample1_cutadapt_secondDeduplication_filtered_2022.fasta;
blastx -db RefSeqPlusNeighbor2023 -query Sample1_cutadapt_secondDeduplication_filtered_2022.fasta -evalue 1e-3 -mt_mode 1 -num_threads 8 -out Sample1_blastx.out;
