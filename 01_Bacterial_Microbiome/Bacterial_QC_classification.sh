#!/bin/bash

##QC bacterial sequencing data. Repeat for each sample.##
bbduk.sh in=Sample1_R1.fastq in2=Sample1_R2.fastq ref=illumina_dna_prep_adapters.fa out=Sample1-adaptTrimQC_R1_2021.fastq out2=Sample1-adaptTrimQC_R2_2021.fastq k=19 hdist=1 ktrim=r qtrim=rl mink=11 trimq=30 minlength=75 minavgquality=20 removeifeitherbad=f otm=t tpe=t overwrite=t 1> Sample1-adaptTrimQC_2021.log.txt 2>&1;
bbduk.sh in=Sample1-adaptTrimQC_R1_2021.fastq in2=Sample1-adaptTrimQC_R2_2021.fastq ref=phix174_ill.ref.fa.gz out=Sample1-R1-phixRemoved_2021.fastq out2=Sample1-R2-phixRemoved_2021.fastq k=31 hdist=1 overwrite=t 1> Sample1phixRemoved_2021.log.txt 2>&1;
bbmap.sh minid=.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 -Xmx64g path=Human_GRCh38 in=Sample1-R1-phixRemoved_2021.fastq in2=Sample1-R2-phixRemoved_2021.fastq outu=Sample1hostRemoved_2021.fastq outm=Sample1hostMatched_2021.fastq 1>Sample1hostRemoval.log_2021.txt 2>&1;
dedupe.sh in=Sample1hostRemoved_2021.fastq out=Sample1firstDeduplication_2021.fastq outd=Sample1firstDuplication_2021.fastq csf=dedupe.cluster.stats overwrite=t minidentity=99 1> Sample1firstDeduplication.log_2021.txt 2>&1;
bbmerge.sh in=Sample1firstDeduplication_2021.fastq out=Sample1firstDeduplicationMerged_2021.fastq outu=Sample1firstDeduplicationUnMerged_2021.fastq 1>Sample1firstDeduplicationMerged.log_2021.txt 2>&1;
cat Sample1firstDeduplicationMerged_2021.fastq Sample1firstDeduplicationUnMerged_2021.fastq > Sample1firstDeduplicationMerged_UnMerged_2021.fastq;
dedupe.sh in=Sample1firstDeduplicationMerged_UnMerged_2021.fastq out=Sample1secondDeduplication_2021.fastq outd=Sample1secondtDuplication_2021.fastq csf=dedupe.cluster.stats overwrite=t minidentity=100 ac=f 1> Sample1secondDeduplication.log_2021.txt 2>&1;
bbduk.sh in=Sample1secondDeduplication_2021.fastq out=Sample1secondDeduplication_filtered_2021.fastq minlength=75 overwrite=t 1>Sample1secondDeduplication_filtered.log_2021.txt 2>&1; 
sed -n '1~4s/^@/>/p;2~4p' Sample1secondDeduplication_filtered_2021.fastq > Sample1secondDeduplication_filtered_2021.fasta;

##Convert previously built Kraken2 database to KrakenUniq database. Archaeal, bacterial, viral, plasmid, human and fungal libraries downloaded and built using Kraken2 in January 2022. Convert to KrakenUniq database following instructions at https://github.com/larssnip/HumGut.##

##Classify clean bacterial reads using KrakenUniq. Repeat for remaining samples.##
krakenuniq --preload --threads 15 --db krakenuniqDB --report-file Sample1_krakenuniq_report.txt --output Sample1_krakenuniq_output.txt Sample1secondDeduplication_filtered_2021.fasta;
krakenuniq --threads 15 --db krakenuniqDB --report-file Sample2_krakenuniq_report.txt --output Sample2_krakenuniq_output.txt Sample2secondDeduplication_filtered_2021.fasta;

##Convert to metaphlan (mpa) format. Use kreport2mpa.py by Jennifer Lu (https://github.com/jenniferlu717/KrakenTools), modified to accept krakenuniq report as input and output read counts and unique k-mer counts.##
python kreport2mpa.py -r Sample1_krakenuniq_report.txt -o mpa/Sample1_krakenuniq_mpa.txt --read_count --intermediate-ranks --display-header;

##For taxa in each sample, convert unique k-mer counts into unique k-mers per million QC reads. In each sample, mask read counts for taxa with <1800 unique k-mers per million reads.##
##Combine masked mpa files for each sample into a single feature table using combine_mpa.py by Jennifer Lu (https://github.com/jenniferlu717/KrakenTools).##
##Parse feature table to species-level and drop false positives/non-bacterial taxa using parse_krakenuniqreports_mpa.py. The parsed feature table is used as input for downstream taxonomy-based analyses (see bacterial diversity R script).##

##Functional profiling with HUMAnN 3 on unmerged, length-filtered paired reads. Repeat for each sample.##
bbduk.sh in=Sample1firstDeduplication_2021.fastq out=Sample1firstDeduplication_filtered_2021.fastq minlength=75 overwrite=t 1>Sample1firstDeduplication_filtered.log_2021.txt 2>&1;
humann --input Sample1firstDeduplication_filtered_2021.fastq --output humann_v3.6/Sample1 --metaphlan-options "--bowtie2db metaphlanDB" --threads 28 1> Sample1_humann3_log.txt 2>&1;
humann_renorm_table --input humann_v3.6/Sample1/Sample1firstDeduplication_filtered_2021_genefamilies.tsv --units relab --output humann_v3.6/Sample1/Sample1_genefamilies_norm_relab.tsv --update-snames 1> Sample1_humann_genefamilies_renorm.log.txt 2>&1;
humann_renorm_table --input humann_v3.6/Sample1/Sample1firstDeduplication_filtered_2021_pathabundance.tsv --units relab --output humann_v3.6/Sample1/Sample1_pathabundance_norm_relab.tsv --update-snames 1> Sample1_humann_pathabundance_renorm.log.txt 2>&1;

##Join individual functional profiles into feature tables. Use the feature tables for downstream functional profiling analyses.##
humann_join_tables -i humann_v3.6/ -o humann_v3.6/humann3_genefamilies_relab_joined.tsv --file_name _genefamilies_norm_relab.tsv -s
humann_join_tables -i humann_v3.6/ -o humann_v3.6/humann3_pathabundance_relab_joined.tsv --file_name _pathabundance_norm_relab.tsv -s
humann_join_tables -i humann_v3.6/ -o humann_v3.6/humann3_genefamilies_RPK_joined.tsv --file_name _genefamilies.tsv -s
humann_join_tables -i humann_v3.6/ -o humann_v3.6/humann3_pathabundance_RPK_joined.tsv --file_name _pathabundance.tsv -s