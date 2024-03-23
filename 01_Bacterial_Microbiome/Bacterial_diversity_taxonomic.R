##Contamination analysis with decontam. Input: feature table from KrakenUniq, taxa with <1800 unique k-mers per million QC reads masked, parsed to species level.##
library(decontam)
data <- read.delim("parsed_krakenuniq_mpa_masked.txt", row.names=1)
metadata<-read.delim("metadata_decontam.txt", row.names = 1)
data<-as.matrix(t(data))
contam.default<-isContaminant(data, method = 'prevalence', neg =metadata$neg_ctrl, threshold=0.1)
clean.default <- t(data[,!contam.default$contaminant=="TRUE"])
clean.default <- clean.default[,!metadata$any_ctrl=="TRUE"]
clean.default <- clean.default[rowSums(clean.default) > 0,]
write.table(clean.default,"krakenUniq_clean_species_default.txt", sep="\t", row.names = TRUE, col.names = NA)

##K-means clustering. Input: clean, species-level feature table with read counts in each sample expressed as relative abundance.##
library(factoextra)
library(cluster)
data<-read.delim("krakenUniq_clean_species_relativeAbundance.txt", header = TRUE, row.names = 1)
data.transposed<-t(data)
opt.clust<-fviz_nbclust(data.transposed, kmeans, k.max = 25, method = "wss") #Find optimal number of clusters by within sum of squares.
opt.clust
set.seed(100)
gap_stat <- clusGap(data.transposed, FUN = kmeans, nstart = 25, K.max = 25, B = 50) #Find optimal number of clusters by gap statistic.
fviz_gap_stat(gap_stat)
sil.kmeans <- fviz_nbclust(data.transposed, kmeans, k.max = 25, method = "silhouette") #Find optimal number of clusters by silhouette score.
sil.kmeans
set.seed(100)
kmean6<-kmeans(data.transposed, 6, iter.max = 10, nstart = 25) #K-means with 6 clusters.
kmean6

##Multinomial logistic regression to test association of clusters with sample metadata. Input: sample metadata including k-means cluster assignment.##
library(mclogit)
data<-read.delim("kmeans_clustering_mblogit.txt", header = TRUE, sep = "\t")
data$Cluster<-as.factor(data$Cluster)
sub1<-subset(data,Cluster =='C2'|Cluster =='C3'|Cluster =='C4'|Cluster =='C5'|Cluster =='C6')
sub2<-subset(data,Cluster =='C3'|Cluster =='C4'|Cluster =='C5'|Cluster =='C6')
sub3<-subset(data,Cluster =='C4'|Cluster =='C5'|Cluster =='C6')
sub4<-subset(data,Cluster =='C5'|Cluster =='C6')
  #By sample type:
set.seed(100)
mbl_type <- mblogit(Cluster~Sample_type,  data = data)
set.seed(100)
mbl_type_sub1 <- mblogit(Cluster~Sample_type,  data = sub1)
set.seed(100)
mbl_type_sub2 <- mblogit(Cluster~Sample_type,  data = sub2)
set.seed(100)
mbl_type_sub3 <- mblogit(Cluster~Sample_type,  data = sub3)
set.seed(100)
mbl_type_sub4 <- mblogit(Cluster~Sample_type,  data = sub4)
summary_type <- summary(mbl_type)
summary_type_sub1 <- summary(mbl_type_sub1)
summary_type_sub2 <- summary(mbl_type_sub2)
summary_type_sub3 <- summary(mbl_type_sub3)
summary_type_sub4 <- summary(mbl_type_sub4)
type_c1ref <- as.data.frame(summary_type$coefficients)
type_c1ref$reference <- c("C1")
type_c1ref$comparison <- rownames(type_c1ref)
type_c2ref <- as.data.frame(summary_type_sub1$coefficients)
type_c2ref$reference <- c("C2")
type_c2ref$comparison <- rownames(type_c2ref)
type_c3ref <- as.data.frame(summary_type_sub2$coefficients)
type_c3ref$reference <- c("C3")
type_c3ref$comparison <- rownames(type_c3ref)
type_c4ref <- as.data.frame(summary_type_sub3$coefficients)
type_c4ref$reference <- c("C4")
type_c4ref$comparison <- rownames(type_c4ref)
type_c5ref <- as.data.frame(summary_type_sub4$coefficients)
type_c5ref$reference <- c("C5")
type_c5ref$comparison <- rownames(type_c5ref)
statsdata_type <- rbind(type_c1ref,type_c2ref,type_c3ref,type_c4ref,type_c5ref)
statsdata_type <- statsdata_type[!grepl("Intercept",statsdata_type$comparison),]
statsdata_type$p.adjusted <- p.adjust(statsdata_type$`Pr(>|z|)`, method = "BH")
  #By time and shedder/non-shedder status with interaction term:
set.seed(100)
mbl <- mblogit(Cluster~Shedder_or_control*Days_since_visit1+Shedder_or_control+Days_since_visit1,  data = data, random=~1|PatientID)
set.seed(100)
mbl_sub1 <- mblogit(Cluster~Shedder_or_control*Days_since_visit1+Shedder_or_control+Days_since_visit1,  data = sub1, random=~1|PatientID)
set.seed(100)
mbl_sub2 <- mblogit(Cluster~Shedder_or_control*Days_since_visit1+Shedder_or_control+Days_since_visit1,  data = sub2, random=~1|PatientID)
set.seed(100)
mbl_sub3 <- mblogit(Cluster~Shedder_or_control*Days_since_visit1+Shedder_or_control+Days_since_visit1,  data = sub3, random=~1|PatientID)
set.seed(100)
mbl_sub4 <- mblogit(Cluster~Shedder_or_control*Days_since_visit1+Shedder_or_control+Days_since_visit1,  data = sub4, random=~1|PatientID)
summary(mbl)
summary(mbl_sub1)
summary(mbl_sub2)
summary(mbl_sub3)
summary(mbl_sub4)
  #By time and shedder/non-shedder status, no interaction term:
set.seed(100)
mbl_noInt <- mblogit(Cluster~Shedder_or_control+Days_since_visit1,  data = data, random=~1|PatientID)
set.seed(100)
mbl_noInt_sub1 <- mblogit(Cluster~Shedder_or_control+Days_since_visit1,  data = sub1, random=~1|PatientID)
set.seed(100)
mbl_noInt_sub2 <- mblogit(Cluster~Shedder_or_control+Days_since_visit1,  data = sub2, random=~1|PatientID)
set.seed(100)
mbl_noInt_sub3 <- mblogit(Cluster~Shedder_or_control+Days_since_visit1,  data = sub3, random=~1|PatientID)
set.seed(100)
mbl_noInt_sub4 <- mblogit(Cluster~Shedder_or_control+Days_since_visit1,  data = sub4, random=~1|PatientID)
summary <- summary(mbl_noInt)
summary_sub1 <- summary(mbl_noInt_sub1)
summary_sub2 <- summary(mbl_noInt_sub2)
summary_sub3 <- summary(mbl_noInt_sub3)
summary_sub4 <- summary(mbl_noInt_sub4)
c1ref <- as.data.frame(summary$coefficients)
c1ref$reference <- c("C1")
c1ref$comparison <- rownames(c1ref)
c2ref <- as.data.frame(summary_sub1$coefficients)
c2ref$reference <- c("C2")
c2ref$comparison <- rownames(c2ref)
c3ref <- as.data.frame(summary_sub2$coefficients)
c3ref$reference <- c("C3")
c3ref$comparison <- rownames(c3ref)
c4ref <- as.data.frame(summary_sub3$coefficients)
c4ref$reference <- c("C4")
c4ref$comparison <- rownames(c4ref)
c5ref <- as.data.frame(summary_sub4$coefficients)
c5ref$reference <- c("C5")
c5ref$comparison <- rownames(c5ref)
statsdata <- rbind(c1ref,c2ref,c3ref,c4ref,c5ref)
statsdata <- statsdata[!grepl("Intercept",statsdata$comparison),]
statsdata$p.adjusted <- p.adjust(statsdata$`Pr(>|z|)`, method = "BH")
  #By discordant timepoint, shedders only:
shed<-subset(data,Shedder_or_control =='Shedder')
shed_sub1<-subset(shed,Cluster =='C2'|Cluster =='C3'|Cluster =='C4'|Cluster =='C5'|Cluster =='C6')
#shed.sub2<-subset(shed,Cluster =='C3'|Cluster =='C4'|Cluster =='C5'|Cluster =='C6') #No shedders in cluster 3.
shed_sub3<-subset(shed,Cluster =='C4'|Cluster =='C5'|Cluster =='C6')
shed_sub4<-subset(shed,Cluster =='C5'|Cluster =='C6')
set.seed(100)
mbl_shed <- mblogit(Cluster~sample_DSC,  data = shed, random=~1|PatientID)
set.seed(100)
mbl_shed_sub1 <- mblogit(Cluster~sample_DSC,  data = shed_sub1, random=~1|PatientID)
set.seed(100)
mbl_shed_sub3 <- mblogit(Cluster~sample_DSC,  data = shed_sub3, random=~1|PatientID)
set.seed(100)
mbl_shed_sub4 <- mblogit(Cluster~sample_DSC,  data = shed_sub4, random=~1|PatientID)
sum_shed <- summary(mbl_shed)
sum_shed1 <- summary(mbl_shed_sub1)
sum_shed3 <- summary(mbl_shed_sub3)
sum_shed4 <- summary(mbl_shed_sub4)
shed_c1ref <- as.data.frame(sum_shed$coefficients)
shed_c1ref$reference <- c("C1")
shed_c1ref$comparison <- rownames(shed_c1ref)
shed_c2ref <- as.data.frame(sum_shed1$coefficients)
shed_c2ref$reference <- c("C2")
shed_c2ref$comparison <- rownames(shed_c2ref)
shed_c4ref <- as.data.frame(sum_shed3$coefficients)
shed_c4ref$reference <- c("C4")
shed_c4ref$comparison <- rownames(shed_c4ref)
shed_c5ref <- as.data.frame(sum_shed4$coefficients)
shed_c5ref$reference <- c("C5")
shed_c5ref$comparison <- rownames(shed_c5ref)
statsdata_shed <- rbind(shed_c1ref,shed_c2ref,shed_c4ref,shed_c5ref)
statsdata_shed <- statsdata_shed[!grepl("Intercept",statsdata_shed$comparison),]
statsdata_shed$p.adjusted <- p.adjust(statsdata_shed$`Pr(>|z|)`, method = "BH")

##Calculate Shannon index. Input: clean, species-level feature table with read counts in each sample normalized to species reads per 100,000 QC reads.##
library(vegan)
data<-read.delim("krakenUniq_clean_species_normalized.txt", row.names = 1)
dataTransposed<-t(data)
shannon<-diversity(dataTransposed, index = 'shannon')

##Calculate Bray-Curtis and Sorensen dissimilarity. Input: clean, species-level feature table with read counts in each sample expressed as relative abundance.##
library(vegan)
data<-read.delim("krakenUniq_clean_species_relativeAbundance.txt", row.names = 1)
data.transposed<-t(data)
dis.bc <- vegdist(data.transposed, method = "bray")
data.transposed[data.transposed > 0] <- 1
dis.sor <- vegdist(data.transposed, method = "bray")

##Mann-Whitney tests of Bray-Curtis dissimilarity within and between groups. Input: table with Bray-Curtis dissimilarity and metadata for sample pairs.##
stats.data <- read.delim("BC_long_format_metadata.txt", header = TRUE, sep = "\t")
w_b.ctrls.all<-subset(stats.data,comparison_all =='Control.Control.within'|comparison_all =='Control.Control.between')
w_b.shed.all<-subset(stats.data,comparison_all =='Shedder.Shedder.within'|comparison_all =='Shedder.Shedder.between')
w.ctrl_shed.all<-subset(stats.data,comparison_all =='Control.Control.within'|comparison_all =='Shedder.Shedder.within')
b.ctrl_shed.all<-subset(stats.data,comparison_all =='Control.Control.between'|comparison_all =='Shedder.Shedder.between')
w.shed.DSC_nonDSC<-subset(stats.data,comparison_DSC_vs_notDSC =='Shedder.within.nonDSC.nonDSC'|comparison_DSC_vs_notDSC =='Shedder.within.DSC.nonDSC')
w.ctrl.DSC_nonDSC<-subset(stats.data,comparison_DSC_vs_notDSC =='Control.within.nonDSC.nonDSC'|comparison_DSC_vs_notDSC =='Control.within.DSC.nonDSC')
mw.w_b.ctrls.all<-wilcox.test(bray_curtis~comparison_all, data = w_b.ctrls.all, alternative = c("two.sided"))
mw.w_b.shed.all<-wilcox.test(bray_curtis~comparison_all, data = w_b.shed.all, alternative = c("two.sided"))
mw.w.ctrl_shed.all<-wilcox.test(bray_curtis~comparison_all, data = w.ctrl_shed.all, alternative = c("two.sided"))
mw.b.ctrl_shed.all<-wilcox.test(bray_curtis~comparison_all, data = b.ctrl_shed.all, alternative = c("two.sided"))
mw.w.shed.DSC_nonDSC<-wilcox.test(bray_curtis~comparison_DSC_vs_notDSC, data = w.shed.DSC_nonDSC, alternative = c("two.sided"))
mw.w.ctrl.DSC_nonDSC<-wilcox.test(bray_curtis~comparison_DSC_vs_notDSC, data = w.ctrl.DSC_nonDSC, alternative = c("two.sided"))
bc.test<-c("mw.w_b.ctrls.all","mw.w_b.shed.all","mw.w.ctrl_shed.all","mw.b.ctrl_shed.all","mw.w.shed.DSC_nonDSC","mw.w.ctrl.DSC_nonDSC")
p.val<-c(mw.w_b.ctrls.all$p.value,mw.w_b.shed.all$p.value,mw.w.ctrl_shed.all$p.value,mw.b.ctrl_shed.all$p.value,mw.w.shed.DSC_nonDSC$p.value,mw.w.ctrl.DSC_nonDSC$p.value)
bc.pvals<-data.frame(bc.test,p.val)
pval.adjust <- p.adjust(bc.pvals$p.val, method = "BH")

##Mann-Whitney tests of Sorensen dissimilarity within and between groups. Input: table with Sorensen dissimilarity and metadata for sample pairs.##
stats.data.sor <- read.delim("sorensen_long_format_metadata.txt", header = TRUE, sep = "\t")
w_b.ctrls.all.sor<-subset(stats.data.sor,comparison_all =='Control.Control.within'|comparison_all =='Control.Control.between')
w_b.shed.all.sor<-subset(stats.data.sor,comparison_all =='Shedder.Shedder.within'|comparison_all =='Shedder.Shedder.between')
w.ctrl_shed.all.sor<-subset(stats.data.sor,comparison_all =='Control.Control.within'|comparison_all =='Shedder.Shedder.within')
b.ctrl_shed.all.sor<-subset(stats.data.sor,comparison_all =='Control.Control.between'|comparison_all =='Shedder.Shedder.between')
w.shed.DSC_nonDSC.sor<-subset(stats.data.sor,comparison_DSC_vs_notDSC =='Shedder.within.nonDSC.nonDSC'|comparison_DSC_vs_notDSC =='Shedder.within.DSC.nonDSC')
w.ctrl.DSC_nonDSC.sor<-subset(stats.data.sor,comparison_DSC_vs_notDSC =='Control.within.nonDSC.nonDSC'|comparison_DSC_vs_notDSC =='Control.within.DSC.nonDSC')
mw.w_b.ctrls.all.sor<-wilcox.test(sorensen~comparison_all, data = w_b.ctrls.all.sor, alternative = c("two.sided"))
mw.w_b.shed.all.sor<-wilcox.test(sorensen~comparison_all, data = w_b.shed.all.sor, alternative = c("two.sided"))
mw.w.ctrl_shed.all.sor<-wilcox.test(sorensen~comparison_all, data = w.ctrl_shed.all.sor, alternative = c("two.sided"))
mw.b.ctrl_shed.all.sor<-wilcox.test(sorensen~comparison_all, data = b.ctrl_shed.all.sor, alternative = c("two.sided"))
mw.w.shed.DSC_nonDSC.sor<-wilcox.test(sorensen~comparison_DSC_vs_notDSC, data = w.shed.DSC_nonDSC.sor, alternative = c("two.sided"))
mw.w.ctrl.DSC_nonDSC.sor<-wilcox.test(sorensen~comparison_DSC_vs_notDSC, data = w.ctrl.DSC_nonDSC.sor, alternative = c("two.sided"))
sor.test<-c("mw.w_b.ctrls.all.sor","mw.w_b.shed.all.sor","mw.w.ctrl_shed.all.sor","mw.b.ctrl_shed.all.sor","mw.w.shed.DSC_nonDSC.sor","mw.w.ctrl.DSC_nonDSC.sor")
pval.initial<-c(mw.w_b.ctrls.all.sor$p.value,mw.w_b.shed.all.sor$p.value,mw.w.ctrl_shed.all.sor$p.value,mw.b.ctrl_shed.all.sor$p.value,mw.w.shed.DSC_nonDSC.sor$p.value,mw.w.ctrl.DSC_nonDSC.sor$p.value)
sor.pvals<-data.frame(sor.test,pval.initial)
pval.adjust.sor <- p.adjust(sor.pvals$pval.initial, method = "BH")

##PCoA. Input: clean, species-level feature table with read counts in each sample expressed as relative abundance.##
library(vegan)
library(phyloseq)
library(ggplot2)
feature_table<-read.delim("krakenUniq_clean_species_relativeAbundance.txt", row.names = 1)
feature_table <- as.data.frame(feature_table)
dim(feature_table)
OTU=otu_table(feature_table,taxa_are_rows =TRUE)
sample_names(OTU)
metadata<-import_qiime_sample_data("metadata.txt")
sample_names(metadata)
identical(colnames(feature_table),rownames(metadata))
physeq =             phyloseq(OTU,metadata)
  #All samples:
set.seed(100)
ord <- ordinate(physeq,
                method = "PCoA",
                distance = "bray")
PCoA.BC  <- plot_ordination(physeq = physeq,
                            ordination = ord,
                            shape = "Shedder_or_control",
                            color = "Days_since_visit1",
                            axes = c(1,2),
                            title='Bray Curtis PCoA') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_color_gradientn(colors = c("#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a"), breaks = c(0,180,360,540,720,750))
PCoA.BC
  #Shedders only:
physeq_shed <- subset_samples(physeq, Shedder_or_control=="Shedder")
set.seed(100)
ord_shed <- ordinate(physeq_shed,
                     method = "PCoA",
                     distance = "bray")
PCoA.BC.shed  <- plot_ordination(physeq = physeq_shed,
                                 ordination = ord_shed,
                                 shape = "DSC_all",
                                 color = "DSC_all",
                                 axes = c(1,2),
                                 title='Bray Curtis PCoA') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_shape_manual(values = c(17,17))
PCoA.BC.shed

##PERMANOVA. Input: clean, species-level feature table with read counts in each sample expressed as relative abundance.##
library(vegan)
data<-read.delim("krakenUniq_clean_species_relativeAbundance.txt", row.names = 1)
data.transposed<-t(data)
dis.bc <- vegdist(data.transposed, method = "bray")
dis.bc2<-as.matrix(dis.bc)
metadata <- read.delim("metadata.txt", header = TRUE, sep = '\t')
  #All samples by sample type, weighted:
set.seed(100)
adonis_type <- adonis2(dis.bc ~ Sample_type, data = metadata, permutations = 999, by = "margin")
adonis_type
  #All samples by sample type, unweighted:
data.transposed[data.transposed > 0] <- 1
dis.sor <- vegdist(data.transposed, method = "bray")
dis.sor2 <-as.matrix(dis.sor)
set.seed(100)
adonis_type_sor <- adonis2(dis.sor ~ Sample_type, data = metadata, permutations = 999, by = "margin")
adonis_type_sor

  #All samples by time and shedder/non-shedder status with interaction term, weighted:
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata$Patient_factor))
set.seed(100)
adonis_int <- adonis2(dis.bc ~ Days_since_visit1*Shedder_or_control + Days_since_visit1 + Shedder_or_control, data = metadata, permutations = perm, by = "margin")
adonis_int
  #All samples, no interaction term, weighted:
set.seed(100)
adonis_noInt <- adonis2(dis.bc ~ Days_since_visit1 + Shedder_or_control, data = metadata, permutations = perm, by = "margin")
adonis_noInt
  #Shedders only, by discordant shedding timepoint, weighted:
sub_metadata <- metadata[metadata$Shedder_or_control=="Shedder",]
sub_BC <- dis.bc2[,colnames(dis.bc2) %in% sub_metadata$Concat_bacterial_ID]
sub_BC2 <- sub_BC[rownames(sub_BC) %in% sub_metadata$Concat_bacterial_ID,]
dis.bc.shed <- as.dist(sub_BC2)
set.seed(100)
perm2 <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_metadata$Patient_factor))
set.seed(100)
adonis_DSC <- adonis2(dis.bc.shed ~ DSC_all, data = sub_metadata, permutations = perm2, by = "margin")
adonis_DSC

##MaAsLin2. Input: clean, species-level feature table with read counts in each sample expressed as relative abundance.##
library(Maaslin2)
feature_table<-read.delim("krakenUniq_clean_species_relativeAbundance.txt", row.names = 1)
feature_table<-as.data.frame(t(feature_table))
metadata<-read.delim("metadata.txt", row.names = 1)
  #Species relative abundance vs time and shedder/nonshedder status with interaction term:
metadata$shedder_time = (metadata$Shedder_or_control == "Shedder") * metadata$Days_since_visit1
metadata$control_time = (metadata$Shedder_or_control == "Control") * metadata$Days_since_visit1
fit_data_interaction= Maaslin2(
  input_data = feature_table, 
  input_metadata = metadata,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_interaction",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1", "Shedder_or_control","shedder_time"),
  random_effects = c("Patient_factor"))
fit_data_interaction2= Maaslin2(
  input_data = feature_table, 
  input_metadata = metadata,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_interaction2",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1", "Shedder_or_control","control_time"),
  random_effects = c("Patient_factor"))
interaction_all_results <- read.delim("maaslin2_interaction/all_results.tsv")
interaction_subset <- subset(interaction_all_results,metadata=="shedder_time")
interaction_subset$qval.readjust <- p.adjust(interaction_subset$pval, method = "BH")
interaction2_all_results <- read.delim("maaslin2_interaction2/all_results.tsv")
interaction2_subset <- subset(interaction2_all_results,metadata=="control_time")
interaction2_subset$qval.readjust <- p.adjust(interaction2_subset$pval, method = "BH")
  #Species relative abundance by time and shedder/nonshedder status, no interaction term:
fit_data_full_model= Maaslin2(
  input_data = feature_table, 
  input_metadata = metadata,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_full_model",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1","Shedder_or_control"),
  random_effects = c("Patient_factor"))
full_model_results <- read.delim("maaslin2_full_model/all_results.tsv")
full_model_time <- subset(full_model_results,metadata=="Days_since_visit1")
full_model_time$qval.readjust <- p.adjust(full_model_time$pval, method = "BH")
full_model_shed_ctrl <- subset(full_model_results,metadata=="Shedder_or_control")
full_model_shed_ctrl$qval.readjust <- p.adjust(full_model_shed_ctrl$pval, method = "BH")
  #Species relative abundance by shedding timepoint, shedders only:
identical(rownames(feature_table),rownames(metadata))
feature_table$Shedder_or_control<-metadata$Shedder_or_control
shedders <- subset(feature_table,Shedder_or_control=="Shedder")
shedders <- shedders[,!names(shedders) %in% c('Shedder_or_control')]
fit_data_DSC= Maaslin2(
  input_data = shedders, 
  input_metadata = metadata,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_DSC",
  transform = "none",
  normalization = "none",
  fixed_effects = c("DSC_all"),
  random_effects = c("Patient_factor"))
