library(phyloseq)
library(ggplot2)
library(doParallel)
library(foreach)
library(ape)
library(dplyr)
library(reshape2)
library(coin)
library(vegan)
library(Maaslin2)
library(nlme)

##Read in anellovirus maximum likelihood tree from IQ-TREE and root using zetatorquevirus reference sequence as outgroup.##
table = read.delim("RCA_anello_relative_abundance_withRef.txt", row.names = 1, sep = "\t")
tree <- read_tree("anello_ORF1.treefile", errorIfNULL = TRUE)
rooted <- root(tree, outgroup = c("Zeta-AB041961-ORF1"), resolve.root = TRUE)
OTU = otu_table(table, taxa_are_rows = TRUE)
physeq = phyloseq(OTU, rooted)
allTaxa = taxa_names(physeq)
cvlTaxa <- allTaxa[!(allTaxa %in% "Zeta-AB041961-ORF1")]
cvlPhyseq = prune_taxa(cvlTaxa, physeq)
cvlPhyseq
is.rooted(phy_tree(cvlPhyseq))
write.tree(rooted,"anello_ORF1_rerooted.treefile")
write.tree(phy_tree(cvlPhyseq),"cvlPhyseq.treefile")

##Calculate weighted UniFrac distance between samples using rooted tree and anellovirus contig relative abundance.##
set.seed(100)
dist <- UniFrac(cvlPhyseq, weighted = TRUE, normalized = TRUE, parallel = TRUE, fast = TRUE)
matrix <- as.matrix(dist)
write.table(matrix, "anello_UniFrac_weighted.txt", quote = FALSE, sep = "\t", col.names = NA)
matrix[lower.tri(matrix, diag = TRUE)] <- NA
matrix.long<-melt(matrix)
matrix.long <- matrix.long[!is.na(matrix.long$value),]
write.table(matrix.long, "anello_UniFrac_weighted_long_format.txt", sep="\t", row.names = FALSE)

##Calculate unweighted UniFrac distance.##
set.seed(100)
dist2 <- UniFrac(cvlPhyseq, weighted = FALSE, normalized = TRUE, parallel = TRUE, fast = TRUE)
matrix2 <- as.matrix(dist2)
write.table(matrix2, "anello_UniFrac_unweighted.txt", quote = FALSE, sep = "\t", col.names = NA)
matrix2[lower.tri(matrix2, diag = TRUE)] <- NA
matrix2.long<-melt(matrix2)
matrix2.long <- matrix2.long[!is.na(matrix2.long$value),]
write.table(matrix2.long, "anello_UniFrac_unweighted_long_format.txt", sep="\t", row.names = FALSE)

##Permutation tests of weighted anellovirus UniFrac distance.##
stats.data <- read.delim("anello_UniFrac_weighted_metadata.txt", header = TRUE, sep = "\t")
stats.data$comparison <- as.factor(stats.data$comparison)
stats.data$block <- c("block1")
stats.data$block <- as.factor(stats.data$block)
w_b.ctrls<-subset(stats.data,comparison =='Control.Control.within'|comparison =='Control.Control.between')
w_b.shed<-subset(stats.data,comparison =='Shedder.Shedder.within'|comparison =='Shedder.Shedder.between')
w.ctrl_shed<-subset(stats.data,comparison =='Control.Control.within'|comparison =='Shedder.Shedder.within')
b.ctrl_shed<-subset(stats.data,comparison =='Control.Control.between'|comparison =='Shedder.Shedder.between')
set.seed(100)
perm.w_b.ctrls<-independence_test(weighted_UniFrac~comparison|block, data = w_b.ctrls, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w_b.shed<-independence_test(weighted_UniFrac~comparison|block, data = w_b.shed, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w.ctrl_shed<-independence_test(weighted_UniFrac~comparison|block, data = w.ctrl_shed, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.b.ctrl_shed<-independence_test(weighted_UniFrac~comparison|block, data = b.ctrl_shed, alternative = c("two.sided"), distribution=approximate(nresample=100000))
weighted.test<-c("perm.w_b.ctrls","perm.w_b.shed","perm.w.ctrl_shed","perm.b.ctrl_shed")
p.val<-c(pvalue(perm.w_b.ctrls),pvalue(perm.w_b.shed),pvalue(perm.w.ctrl_shed),pvalue(perm.b.ctrl_shed))
weighted.pvals<-data.frame(weighted.test,p.val)
#Benjamini-Hochberg correction for multiple comparisons.#
pval.adjust <- p.adjust(weighted.pvals$p.val, method = "BH")
weighted.adjust<-cbind(weighted.pvals,pval.adjust)
write.table(weighted.adjust, "anello_weight_unifrac_permutation_adjust1.txt",row.names = FALSE,sep = "\t")

##Permutation tests of unweighted anellovirus UniFrac distance.##
stats.data2 <- read.delim("anello_UniFrac_unweighted_metadata.txt", header = TRUE, sep = "\t")
stats.data2$comparison <- as.factor(stats.data2$comparison)
stats.data2$block <- c("block1")
stats.data2$block <- as.factor(stats.data2$block)
w_b.ctrls.2<-subset(stats.data2,comparison =='Control.Control.within'|comparison =='Control.Control.between')
w_b.shed.2<-subset(stats.data2,comparison =='Shedder.Shedder.within'|comparison =='Shedder.Shedder.between')
w.ctrl_shed.2<-subset(stats.data2,comparison =='Control.Control.within'|comparison =='Shedder.Shedder.within')
b.ctrl_shed.2<-subset(stats.data2,comparison =='Control.Control.between'|comparison =='Shedder.Shedder.between')
set.seed(100)
perm.w_b.ctrls.unw<-independence_test(unweighted_UniFrac~comparison|block, data = w_b.ctrls.2, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w_b.shed.unw<-independence_test(unweighted_UniFrac~comparison|block, data = w_b.shed.2, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w.ctrl_shed.unw<-independence_test(unweighted_UniFrac~comparison|block, data = w.ctrl_shed.2, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.b.ctrl_shed.unw<-independence_test(unweighted_UniFrac~comparison|block, data = b.ctrl_shed.2, alternative = c("two.sided"), distribution=approximate(nresample=100000))
unweighted.test<-c("perm.w_b.ctrls.unw","perm.w_b.shed.unw","perm.w.ctrl_shed.unw","perm.b.ctrl_shed.unw")
p.val.unw<-c(pvalue(perm.w_b.ctrls.unw),pvalue(perm.w_b.shed.unw),pvalue(perm.w.ctrl_shed.unw),pvalue(perm.b.ctrl_shed.unw))
unweighted.pvals<-data.frame(unweighted.test,p.val.unw)
#Benjamini-Hochberg correction for multiple comparisons.#
pval.adjust.unw <- p.adjust(unweighted.pvals$p.val.unw, method = "BH")
unweighted.adjust<-cbind(unweighted.pvals,pval.adjust.unw)
write.table(unweighted.adjust, "anello_unweight_permutation_adjust1.txt",row.names = FALSE,sep = "\t")

##PCoA, weighted, all samples.##
metadata<-import_qiime_sample_data("Virome_Metadata_anello.txt")
sample_names(metadata)
identical(sample_names(OTU),sample_names(metadata))
cvlPhyseqMeta = merge_phyloseq(cvlPhyseq,metadata)
cvlPhyseqMeta
set.seed(100)
ord <- ordinate(cvlPhyseqMeta,
                method = "PCoA",
                distance = "unifrac",
                weighted = TRUE)
PCoA  <- plot_ordination(physeq = cvlPhyseqMeta,
                         ordination = ord,
                         shape = "Shedder_or_control", # metadata variable
                         color = "Days_since_visit1", # metadata variable
                         axes = c(1,2),
                         title='Anellovirus Weighted UniFrac PCoA') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_color_gradientn(colors = c("#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a"), breaks = c(0,180,360,540,720,750))
PCoA

##PCoA, unweighted, all samples.##
set.seed(100)
ord_unw <- ordinate(cvlPhyseqMeta,
                    method = "PCoA",
                    distance = "unifrac",
                    weighted = FALSE)
PCoA_unw <- plot_ordination(physeq = cvlPhyseqMeta,
                            ordination = ord_unw,
                            shape = "Shedder_or_control", # metadata variable
                            color = "Days_since_visit1", # metadata variable
                            axes = c(1,2),
                            title='Anellovirus Unweighted UniFrac PCoA') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_color_gradientn(colors = c("#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a"), breaks = c(0,180,360,540,720,750))
PCoA_unw

##PCoA, weighted, subset by shedder/control status.##
physeq_ctrl <- subset_samples(cvlPhyseqMeta, Shedder_or_control=="Control")
set.seed(100)
ord_ctrl <- ordinate(physeq_ctrl,
                     method = "PCoA",
                     distance = "unifrac",
                     weighted = TRUE)
PCoA_w_ctrl  <- plot_ordination(physeq = physeq_ctrl,
                                ordination = ord_ctrl,
                                shape = "Shedder_or_control", # metadata variable
                                color = "Days_since_visit1", # metadata variable
                                axes = c(1,2),
                                title='Anellovirus Weighted UniFrac PCoA Controls') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_color_gradientn(colors = c("#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a"), breaks = c(0,180,360,540,720,750))
PCoA_w_ctrl

physeq_shed <- subset_samples(cvlPhyseqMeta, Shedder_or_control=="Shedder")
set.seed(100)
ord_shed <- ordinate(physeq_shed,
                     method = "PCoA",
                     distance = "unifrac",
                     weighted = TRUE)
PCoA_w_shed  <- plot_ordination(physeq = physeq_shed,
                                ordination = ord_shed,
                                shape = "DSC_all", # metadata variable
                                color = "Days_since_visit1", # metadata variable
                                axes = c(1,2),
                                title='Anellovirus Weighted UniFrac PCoA Shedders') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_color_gradientn(colors = c("#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a"), breaks = c(0,180,360,540,720,750))+
  scale_shape_manual(values = c(17,17))
PCoA_w_shed

##PCoA, unweighted, subset to controls or shedders##
set.seed(100)
ord_ctrl_unw <- ordinate(physeq_ctrl,
                         method = "PCoA",
                         distance = "unifrac",
                         weighted = FALSE)
PCoA_unw_ctrl <- plot_ordination(physeq = physeq_ctrl,
                                 ordination = ord_ctrl_unw,
                                 shape = "Shedder_or_control", # metadata variable
                                 color = "Days_since_visit1", # metadata variable
                                 axes = c(1,2),
                                 title='Anellovirus Unweighted UniFrac PCoA Controls') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_color_gradientn(colors = c("#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a"), breaks = c(0,180,360,540,720,750))
PCoA_unw_ctrl

set.seed(100)
ord_shed_unw <- ordinate(physeq_shed,
                         method = "PCoA",
                         distance = "unifrac",
                         weighted = FALSE)
PCoA_unw_shed <- plot_ordination(physeq = physeq_shed,
                                 ordination = ord_shed_unw,
                                 shape = "DSC_all", # metadata variable
                                 color = "Days_since_visit1", # metadata variable
                                 axes = c(1,2),
                                 title='Anellovirus Unweighted UniFrac PCoA Shedders') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_color_gradientn(colors = c("#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a"), breaks = c(0,180,360,540,720,750))
PCoA_unw_shed

##PERMANOVA##
adonis.metadata <- read.delim("Virome_Metadata_anello.txt", header = TRUE, sep = '\t')

##All samples, weighted UniFrac, with interaction term.##
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = adonis.metadata$Patient_factor))
set.seed(100)
adonis_int <- adonis2(matrix ~ Days_since_visit1*Shedder_or_control + Days_since_visit1 + Shedder_or_control, data = adonis.metadata, permutations = perm, by = "margin")
adonis_int
write.table(adonis_int, "anello_weighted_UniFrac_adonis2_allSamples_interaction.txt", sep = "\t")
##No interaction term.##
set.seed(100)
adonis_noInt <- adonis2(matrix ~ Days_since_visit1 + Shedder_or_control, data = adonis.metadata, permutations = perm, by = "margin")
adonis_noInt
write.table(adonis_noInt, "anello_weighted_UniFrac_adonis2_allSamples_noInteraction.txt", sep = "\t")

##All samples, unweighted UniFrac, with interaction term.##
set.seed(100)
adonis_unw_int <- adonis2(matrix2 ~ Days_since_visit1*Shedder_or_control + Days_since_visit1 + Shedder_or_control, data = adonis.metadata, permutations = perm, by = "margin")
adonis_unw_int
write.table(adonis_unw_int, "anello_unweighted_UniFrac_adonis2_allSamples_interaction.txt", sep = "\t")
##No interaction term.##
set.seed(100)
adonis_unw <- adonis2(matrix2 ~ Days_since_visit1 + Shedder_or_control, data = adonis.metadata, permutations = perm, by = "margin")
adonis_unw
write.table(adonis_unw, "anello_unweighted_UniFrac_adonis2_allSamples_noInteraction.txt", sep = "\t")

metadata.ctrls <- adonis.metadata[adonis.metadata$Shedder_or_control=="Control",]

##Subset weighted UniFrac, controls.##
ctrls <- matrix[,colnames(matrix) %in% metadata.ctrls$RCA_ID]
ctrls <- ctrls[rownames(ctrls) %in% metadata.ctrls$RCA_ID,]

##Controls over time.##
set.seed(100)
perm.ctrls <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata.ctrls$Patient_factor))
set.seed(100)
adonis_ctrl <- adonis2(ctrls ~ Days_since_visit1, data = metadata.ctrls, permutations = perm.ctrls, by = "margin")
adonis_ctrl
write.table(adonis_ctrl, "anello_weighted_UniFrac_adonis2_controls.txt",sep = "\t")

metadata.shed <- adonis.metadata[adonis.metadata$Shedder_or_control=="Shedder",]

##Subset weighted UniFrac, shedders.##
shed <- matrix[,colnames(matrix) %in% metadata.shed$RCA_ID]
shed <- shed[rownames(shed) %in% metadata.shed$RCA_ID,]

##Shedders over time.##
set.seed(100)
perm.shed <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata.shed$Patient_factor))
set.seed(100)
adonis_shed <- adonis2(shed ~ Days_since_visit1, data = metadata.shed, permutations = perm.shed, by = "margin")
adonis_shed
write.table(adonis_shed, "anello_weighted_UniFrac_adonis2_shedders.txt",sep = "\t")

##Shedders, DSC vs non-DSC.##
set.seed(100)
adonis_DSC <- adonis2(shed ~ DSC_all, data = metadata.shed, permutations = perm.shed, by = "margin")
adonis_DSC
write.table(adonis_DSC, "anello_weighted_UniFrac_adonis2_shedders_DSC.txt",sep = "\t")

##Subset unweighted UniFrac, controls.##
ctrls2 <- matrix2[,colnames(matrix2) %in% metadata.ctrls$RCA_ID]
ctrls2 <- ctrls2[rownames(ctrls2) %in% metadata.ctrls$RCA_ID,]

##Controls over time, unweighted.##
set.seed(100)
adonis_ctrl2 <- adonis2(ctrls2 ~ Days_since_visit1, data = metadata.ctrls, permutations = perm.ctrls, by = "margin")
adonis_ctrl2
write.table(adonis_ctrl2, "anello_unweighted_UniFrac_adonis2_controls.txt",sep = "\t")

##Subset unweighted UniFrac, shedders.##
shed2 <- matrix2[,colnames(matrix2) %in% metadata.shed$RCA_ID]
shed2 <- shed2[rownames(shed2) %in% metadata.shed$RCA_ID,]

##Shedders over time, unweighted.##
set.seed(100)
adonis_shed2 <- adonis2(shed2 ~ Days_since_visit1, data = metadata.shed, permutations = perm.shed, by = "margin")
adonis_shed2
write.table(adonis_shed2, "anello_unweighted_UniFrac_adonis2_shedders.txt",sep = "\t")

##Shedders, DSC vs non-DSC, unweighted.##
set.seed(100)
adonis_DSC2 <- adonis2(shed2 ~ DSC_all, data = metadata.shed, permutations = perm.shed, by = "margin")
adonis_DSC2
write.table(adonis_DSC2, "anello_unweighted_UniFrac_adonis2_shedders_DSC.txt",sep = "\t")

##Differential abundance analysis with MaAsLin2.##
anello_genus<-read.delim("RCA_clean_default_relAbund_anello_genus.txt", row.names = 1)
anello_genus<-as.data.frame(t(anello_genus))
metadata<-read.delim("Virome_Metadata_032323.txt", row.names = 1)

##Make variable to test for interaction.##
metadata$shedder_time = (metadata$Shedder_or_control == "Shedder") * metadata$Days_since_visit1

fit_data_genus_shedInteract= Maaslin2(
  input_data = anello_genus, 
  input_metadata = metadata,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2/anello_only/genus_shedInteract",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1", "Shedder_or_control","shedder_time"),
  random_effects = c("Patient_factor"))

##Read in results, subset interaction variable results, and recalculate q-values.##
shedInteract_all_results <- read.delim("maaslin2/anello_only/genus_shedInteract/all_results.tsv")
shedInteract_subset <- subset(shedInteract_all_results,metadata=="shedder_time")
shedInteract_subset$qval.readjust <- p.adjust(shedInteract_subset$pval, method = "BH")
write.table(shedInteract_subset,"maaslin2/anello_only/genus_shedInteract/shedInteract_subset_results.tsv",sep = "\t",col.names = NA)

##Full model.##
fit_data_full_model= Maaslin2(
  input_data = anello_genus, 
  input_metadata = metadata,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2/anello_only/genus_full_model",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1","Shedder_or_control"),
  random_effects = c("Patient_factor"))

##Subset variables to readjust q values.##
full_model_results <- read.delim("maaslin2/anello_only/genus_full_model/all_results.tsv")
full_model_time <- subset(full_model_results,metadata=="Days_since_visit1")
full_model_time$qval.readjust <- p.adjust(full_model_time$pval, method = "BH")
write.table(full_model_time,"maaslin2/anello_only/genus_full_model/time_all_results_readjusted.tsv",sep = "\t",col.names = NA)
full_model_shed_ctrl <- subset(full_model_results,metadata=="Shedder_or_control")
full_model_shed_ctrl$qval.readjust <- p.adjust(full_model_shed_ctrl$pval, method = "BH")
write.table(full_model_shed_ctrl,"maaslin2/anello_only/genus_full_model/shed_ctrl_all_results_readjusted.tsv",sep = "\t",col.names = NA)

##Subset shedders and non-shedders.##
shed_metadata <- subset(metadata,Shedder_or_control=="Shedder")
shedders_genus <- anello_genus[rownames(anello_genus) %in% rownames(shed_metadata),]
ctrl_metadata <- subset(metadata,Shedder_or_control=="Control")
controls_genus <- anello_genus[rownames(anello_genus) %in% rownames(ctrl_metadata),]

##Time, non-shedders, genus.##
fit_data_time_genus_ctrl= Maaslin2(
  input_data = controls_genus, 
  input_metadata = ctrl_metadata,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2/anello_only/genus_time_ctrls",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1"),
  random_effects = c("Patient_factor"))

##Time, shedders, genus.##
fit_data_time_genus_shed= Maaslin2(
  input_data = shedders_genus, 
  input_metadata = shed_metadata,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2/anello_only/genus_time_shed",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1"),
  random_effects = c("Patient_factor"))

##Linear mixed modeling of CD4 counts vs anellovirus NGS abundance.##
cd4.data <- read.delim("anello_cd4.txt") #Input: anellovirus NGS relative abundance and CD4 counts.
ctrl.cd4 <- subset(cd4.data, Shedder_or_control=="Control")
shed.cd4 <- subset(cd4.data, Shedder_or_control=="Shedder")

#Anelloviridae vs CD4 with interaction term.##
anello_int <- lme(Anelloviridae_relAbund ~ CD4_count + Shedder_or_control + (CD4_count*Shedder_or_control), random =~1|Patient_factor, data = cd4.data)
summary(anello_int)

#Anelloviridae vs CD4 no interaction term. Use this since the interaction term is not significant.#
anello_no_int <- lme(Anelloviridae_relAbund ~ CD4_count, random =~1|Patient_factor, data = cd4.data)
anello_sum <- summary(anello_no_int)
anello_sum$tTable

anello_ctrls <- lme(Anelloviridae_relAbund ~ CD4_count, random =~1|Patient_factor, data = ctrl.cd4)
anello_ctrl_sum <- summary(anello_ctrls)
anello_ctrl_sum$tTable

anello_shed <- lme(Anelloviridae_relAbund ~ CD4_count, random =~1|Patient_factor, data = shed.cd4)
summary(anello_shed)

#Alphatorquevirus vs CD4 with interaction term.##
alpha_int <- lme(Alphatorque_relAbund ~ CD4_count + Shedder_or_control + (CD4_count*Shedder_or_control), random =~1|Patient_factor, data = cd4.data)
summary(alpha_int)

#Alphatorquevirus vs CD4 no interaction term. Use this since the interaction term is not significant.#
alpha_no_int <- lme(Alphatorque_relAbund ~ CD4_count, random =~1|Patient_factor, data = cd4.data)
alpha_sum <- summary(alpha_no_int)
alpha_sum$tTable

alpha_ctrls <- lme(Alphatorque_relAbund ~ CD4_count, random =~1|Patient_factor, data = ctrl.cd4)
alpha_sum_ctrl <- summary(alpha_ctrls)
alpha_sum_ctrl$tTable

alpha_shed <- lme(Alphatorque_relAbund ~ CD4_count, random =~1|Patient_factor, data = shed.cd4)
alpha_sum_shed <- summary(alpha_shed)
alpha_sum_shed$tTable

#Betatorquevirus vs CD4 with interaction term.##
beta_int <- lme(Betatorque_relAbund ~ CD4_count + Shedder_or_control + (CD4_count*Shedder_or_control), random =~1|Patient_factor, data = cd4.data)
summary(beta_int)

#Betatorquevirus vs CD4 no interaction term. Use this since the interaction term is not significant.#
beta_no_int <- lme(Betatorque_relAbund ~ CD4_count, random =~1|Patient_factor, data = cd4.data)
summary(beta_no_int)

#Gammatorquevirus vs CD4 with interaction term.##
gamma_int <- lme(Gammatorque_relAbund ~ CD4_count + Shedder_or_control + (CD4_count*Shedder_or_control), random =~1|Patient_factor, data = cd4.data)
summary(gamma_int)

#Gammatorquevirus vs CD4 no interaction term. Use this since the interaction term is not significant.#
gamma_no_int <- lme(Gammatorque_relAbund ~ CD4_count, random =~1|Patient_factor, data = cd4.data)
summary(gamma_no_int)

##Linear mixed modeling of alphatorquevirus qPCR copy numbers.##
ATV.data <- read.delim("Virome_Metadata_ATV.txt") #Input: alphatorquevirus copy numbers and sample metadata.

#ATV copy numbers vs anellovirus NGS abundance.#
anello <- lm(ATV_copies_ul_TNA ~ Anelloviridae_relAbund, data = ATV.data)
summary(anello)

alpha <- lm(ATV_copies_ul_TNA ~ Alphatorque_relAbund, data = ATV.data)
summary(alpha)

beta <- lm(ATV_copies_ul_TNA ~ Betatorque_relAbund, data = ATV.data)
summary(beta)

gamma <- lm(ATV_copies_ul_TNA ~ Gammatorque_relAbund, data = ATV.data)
summary(gamma)

#ATV copies over time with interaction term.##
ATV_time_int <- lme(ATV_copies_ul_TNA ~ Days_since_visit1 + Shedder_or_control + (Days_since_visit1*Shedder_or_control), random =~1|Patient_factor, data = ATV.data)
summary(ATV_time_int)

ATV_time_no_int <- lme(ATV_copies_ul_TNA ~ Days_since_visit1 + Shedder_or_control, random =~1|Patient_factor, data = ATV.data)
summary(ATV_time_no_int)

ATV_time_only <- lme(ATV_copies_ul_TNA ~ Days_since_visit1, random =~1|Patient_factor, data = ATV.data)
summary(ATV_time_only)

ATV_shed_ctrl_only <- lme(ATV_copies_ul_TNA ~ Shedder_or_control, random =~1|Patient_factor, data = ATV.data)
summary(ATV_shed_ctrl_only)

#Subset samples that have CD4 data available.#
ATV.data.cd4 <- ATV.data[!is.na(ATV.data$CD4_count),]

#ATV copies by CD4 count with interaction term.##
ATV_cd4_int <- lme(ATV_copies_ul_TNA ~ CD4_count + Shedder_or_control + (CD4_count*Shedder_or_control), random =~1|Patient_factor, data = ATV.data.cd4)
summary(ATV_cd4_int)

ATV_cd4_no_int <- lme(ATV_copies_ul_TNA ~ CD4_count + Shedder_or_control, random =~1|Patient_factor, data = ATV.data.cd4)
summary(ATV_cd4_no_int)

ATV_cd4_only <- lme(ATV_copies_ul_TNA ~ CD4_count, random =~1|Patient_factor, data = ATV.data.cd4)
summary(ATV_cd4_only)

ATV_shed_ctrl_subset <- lme(ATV_copies_ul_TNA ~ Shedder_or_control, random =~1|Patient_factor, data = ATV.data.cd4)
summary(ATV_shed_ctrl_subset)

#Non-shedders only.#
cd4_ctrl <- subset(ATV.data.cd4, Shedder_or_control=="Control")
ATV_cd4_ctrl <- lme(ATV_copies_ul_TNA ~ CD4_count, random =~1|Patient_factor, data = cd4_ctrl)
summary(ATV_cd4_ctrl)

#Shedders only.#
cd4_shed <- subset(plot.data.cd4, Shedder_or_control=="Shedder")
ATV_cd4_shed <- lme(ATV_copies_ul_TNA ~ CD4_count, random =~1|Patient_factor, data = cd4_shed)
summary(ATV_cd4_shed)

