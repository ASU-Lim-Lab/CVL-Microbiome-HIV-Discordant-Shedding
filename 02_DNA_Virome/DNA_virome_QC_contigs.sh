#!/bin/bash

##QC MDA sequencing data. Repeat for each sample.##
gunzip Sample1_R1_001.fastq.gz;
gunzip Sample1_R2_001.fastq.gz;
cutadapt -a file:illumina_dna_prep_adapters.fa -A file:illumina_dna_prep_adapters.fa -o Sample1_cutadapt_adaptTrim_R1_2022.fastq -p Sample1_cutadapt_adaptTrim_R2_2022.fastq Sample1_R1_001.fastq Sample1_R2_001.fastq -j 28 1> Sample1_cutadapt_adaptTrim_2022.log.txt 2>&1;
bbduk.sh in=Sample1_cutadapt_adaptTrim_R1_2022.fastq in2=Sample1_cutadapt_adaptTrim_R2_2022.fastq out=Sample1-cutadapt_QC_R1_2022.fastq out2=Sample1-cutadapt_QC_R2_2022.fastq qtrim=rl trimq=30 minlength=75 minavgquality=20 removeifeitherbad=f otm=t tpe=t overwrite=t 1> Sample1-cutadapt_QC_2022.log.txt 2>&1;
bbduk.sh in=Sample1-cutadapt_QC_R1_2022.fastq in2=Sample1-cutadapt_QC_R2_2022.fastq ref=phix174_ill.ref.fa.gz out=Sample1-R1-cutadapt_phixRemoved_2022.fastq out2=Sample1-R2-cutadapt_phixRemoved_2022.fastq k=31 hdist=1 overwrite=t 1> Sample1_cutadapt_phixRemoved_2022.log.txt 2>&1;
bbmap.sh minid=.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 -Xmx64g path=Human_GRCh38 in=Sample1-R1-cutadapt_phixRemoved_2022.fastq in2=Sample1-R2-cutadapt_phixRemoved_2022.fastq outu=Sample1_cutadapt_hostRemoved_2022.fastq outm=Sample1_cutadapt_hostMatched_2022.fastq 1>Sample1_cutadapt_hostRemoval.log_2022.txt 2>&1;
dedupe.sh in=Sample1_cutadapt_hostRemoved_2022.fastq out=Sample1_cutadapt_firstDeduplication_2022.fastq outd=Sample1_cutadapt_firstDuplication_2022.fastq csf=dedupe.cluster.stats overwrite=t minidentity=99 1> Sample1_cutadapt_firstDeduplication.log_2022.txt 2>&1;
bbmerge.sh in=Sample1_cutadapt_firstDeduplication_2022.fastq out=Sample1_cutadapt_firstDeduplicationMerged_2022.fastq outu=Sample1_cutadapt_firstDeduplicationUnMerged_2022.fastq 1>Sample1_cutadapt_firstDeduplicationMerged.log_2022.txt 2>&1;
cat Sample1_cutadapt_firstDeduplicationMerged_2022.fastq Sample1_cutadapt_firstDeduplicationUnMerged_2022.fastq > Sample1_cutadapt_firstDeduplicationMerged_UnMerged_2022.fastq;
dedupe.sh in=Sample1_cutadapt_firstDeduplicationMerged_UnMerged_2022.fastq out=Sample1_cutadapt_secondDeduplication_2022.fastq outd=Sample1_cutadapt_secondDuplication_2022.fastq csf=dedupe.cluster.stats overwrite=t minidentity=100 ac=f 1> Sample1_cutadapt_secondDeduplication.log_2022.txt 2>&1;
bbduk.sh in=Sample1_cutadapt_secondDeduplication_2022.fastq out=Sample1_cutadapt_secondDeduplication_filtered_2022.fastq minlength=75 overwrite=t 1>Sample1_cutadapt_secondDeduplication_filtered.log_2022.txt 2>&1; 
sed -n '1~4s/^@/>/p;2~4p' Sample1_cutadapt_secondDeduplication_filtered_2022.fastq > Sample1_cutadapt_secondDeduplication_filtered_2022.fasta;

##Build contigs. Input: host-removed reads from line 9. Repeat for each sample.##
metaspades.py -o Sample1 --12 Sample1_cutadapt_hostRemoved_2022.fastq --phred-offset 33 -t 88 -m 1500;
awk '/^>/{print ">I82757_Contig_" ++i; next}{print}' < Sample1/contigs.fasta > Sample1/Sample1_RenamedContigs.fasta;

##Concatenate sample contigs into one file and sequencing control contigs into a separate file. Then deduplicate contigs with CD-HIT-EST.##
cd-hit-est -i MDA_Sample_Contigs.fasta -o MDA_SampleContigs_cd-hit.fasta -c 0.95 -n 10 -G 0 -aS 0.95 -g 1 -r 1 -M 20000 -d 0 -T 28 1> MDA_samples_cd-hit.log.txt 2>&1
cd-hit-est -i MDA_Control_Contigs.fasta -o MDA_ControlContigs_cd-hit.fasta -c 0.95 -n 10 -G 0 -aS 0.95 -g 1 -r 1 -M 20000 -d 0 -T 28 1> MDA_RCA_controls_cd-hit.log.txt 2>&1

#Concatenate deduplicated sample and control contigs into one file and filter to minimum length 500.##
bbduk.sh in=MDA_Sample_Control_Contigs_cd-hit.fasta out=MDA_Sample_Control_Contigs_min500.fasta minlen=500 1> MDA_Sample_Control_Contigs_length_filter.log.txt 2>&1;

##Identify candidate viral contigs with VirSorter2.##
virsorter run -w MDA_vs2 -i MDA_Sample_Control_Contigs_min500.fasta --include-groups "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae" -j 80 all 1> MDA_vs2.log.txt 2>&1

##Identify candidate viral contigs with Cenote-Taker2.##
python run_cenote-taker2.py -c MDA_Sample_Control_Contigs_min500.fasta -r MDA_ct2 -m 32 -t 88 -p false -db standard --minimum_length_linear 500 2>&1 | tee MDA_ct2_log.txt

##Query candidate viral contigs against NCBI NT database to identify human sequences. Remove contigs with hits to human sequences if >=95% identity and query coverage.##
blastn -task megablast -db ntDB2020 -query MDA_candidate_viral_contigs.fasta -evalue 1e-10 -num_threads 56 -max_target_seqs 1 -outfmt "6 qseqid sseqid staxids evalue bitscore pident nident qcovs length mismatch qlen slen" -out MDA_candidate_viral_contigs_megablast_NT.out 

##Use CheckV to check for false positive bacterial sequences.##
checkv end_to_end MDA_candidate_viruses_human_removed.fasta MDA_checkv -t 88 -d checkv-db-v1.5 1> MDA_checkv.log.txt 2>&1

##Keep contigs if complete, high quality, medium quality, low-quality with >=3 viral genes and <=1 host gene, or proviruses. Discard contigs with high kmer frequency (CheckV) or classified as bacteria (Cenote-Taker2).##

##Query viral contigs against reference virus database to assign taxonomy. Use taxonomy of best database hit. Remove contigs classified as Mimiviridae.##
blastx -db RefSeqPlusNeighbor2023 -query MDA_viruses.fasta -evalue 1e-3 -num_threads 88 -max_target_seqs 1 -outfmt "6 qseqid sseqid staxids evalue bitscore pident nident qcovs length mismatch qlen slen" -out MDA_viruses_blastx.out

##Build Bowtie2 database with viral contigs.##
bowtie2-build MDA_final_viruses.fasta MDA_final_viruses_bt2DB 1> MDA_final_viruses_bt2DB.log.txt 2>&1

##Map QC'ed reads to viral contigs. Repeat for each sample and join the outputs into a feature table, which will be the input for decontam (see DNA virome diversity R script).##
bowtie2 -x MDA_final_viruses_bt2DB -U Sample1_cutadapt_secondDeduplication_filtered_2022.fasta -p 88 -f -S Sample1readsMapped_bt2.sam;
samtools view -h -F 0x904 Sample1readsMapped_bt2.sam > Sample1primaryMapped_bt2.sam;
samtools view -h -f 0x4 Sample1readsMapped_bt2.sam > Sample1unmapped_bt2.sam;
samtools view -S -b Sample1primaryMapped_bt2.sam > Sample1primaryMapped_bt2.bam;
samtools sort Sample1primaryMapped_bt2.bam > Sample1primaryMapped_sorted_bt2.bam;
samtools index Sample1primaryMapped_sorted_bt2.bam;
samtools idxstats Sample1primaryMapped_sorted_bt2.bam > Sample1counts_bt2.txt;

