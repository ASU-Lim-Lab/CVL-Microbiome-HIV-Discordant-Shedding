library(phyloseq)
library(ggplot2)
library(doParallel)
library(foreach)
library(ape)
library(dplyr)
library(reshape2)
library(vegan)

##Read in papillomavirus maximum likelihood tree from IQ-TREE and root using mupapillomavirus reference sequence as outgroup.##
table = read.delim("RCA_papilloma_relAbund_withRef.txt", row.names = 1, sep = "\t")
tree <- read_tree("papilloma_E1E2L2L1.treefile", errorIfNULL = TRUE)
rooted <- root(tree, outgroup = c("Mu-NC_038525-E1-E2-L2-L1"), resolve.root = TRUE)
OTU = otu_table(table, taxa_are_rows = TRUE)
physeq = phyloseq(OTU, rooted)
allTaxa = taxa_names(physeq)
cvlTaxa <- allTaxa[!(allTaxa %in% "Mu-NC_038525-E1-E2-L2-L1")]
cvlPhyseq = prune_taxa(cvlTaxa, physeq)
cvlPhyseq
is.rooted(phy_tree(cvlPhyseq))
write.tree(rooted,"papilloma_E1E2L2L1_rerooted.treefile")
write.tree(phy_tree(cvlPhyseq),"papilloma_cvlPhyseq.treefile")

##Calculate weighted UniFrac distance between samples using rooted tree and papillomavirus contig relative abundance.##
set.seed(100)
dist <- UniFrac(cvlPhyseq, weighted = TRUE, normalized = TRUE, parallel = TRUE, fast = TRUE)
matrix <- as.matrix(dist)
write.table(matrix, "papilloma_UniFrac_weighted_070323.txt", quote = FALSE, sep = "\t", col.names = NA)
matrix[lower.tri(matrix, diag = TRUE)] <- NA
matrix.long<-melt(matrix)
matrix.long <- matrix.long[!is.na(matrix.long$value),]
write.table(matrix.long, "papilloma_UniFrac_weighted_long_format.txt", sep="\t", row.names = FALSE)

set.seed(100)
dist2 <- UniFrac(cvlPhyseq, weighted = FALSE, normalized = TRUE, parallel = TRUE, fast = TRUE)
matrix2 <- as.matrix(dist2)
write.table(matrix2, "papilloma_UniFrac_unweighted.txt", quote = FALSE, sep = "\t", col.names = NA)
matrix2[lower.tri(matrix2, diag = TRUE)] <- NA
matrix2.long<-melt(matrix2)
matrix2.long <- matrix2.long[!is.na(matrix2.long$value),]
write.table(matrix2.long, "papilloma_UniFrac_unweighted_long_format.txt", sep="\t", row.names = FALSE)

#Permutation tests of weighted papillomavirus UniFrac distance.#
stats.data <- read.delim("papilloma_UniFrac_weighted_metadata.txt", header = TRUE, sep = "\t")
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
write.table(weighted.adjust, "papilloma_weight_unifrac_permutation_adjust.txt",row.names = FALSE,sep = "\t")

#Permutation tests of unweighted papillomavirus UniFrac distance.#
stats.data2 <- read.delim("papilloma_UniFrac_unweighted_metadata.txt", header = TRUE, sep = "\t")
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
pval.adjust.unw <- p.adjust(unweighted.pvals$p.val.unw, method = "BH")
unweighted.adjust<-cbind(unweighted.pvals,pval.adjust.unw)
write.table(unweighted.adjust, "papilloma_unweight_permutation_adjust.txt",row.names = FALSE,sep = "\t")

##PCoA, weighted, all samples.##
metadata<-import_qiime_sample_data("Virome_Metadata_papilloma.txt")
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
                         title='Papillomavirus Weighted UniFrac PCoA') +
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
                            title='Papillomavirus Unweighted UniFrac PCoA') +
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
                                title='Papillomavirus Weighted UniFrac PCoA Controls') +
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
                                title='Papillomavirus Weighted UniFrac PCoA Shedders') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_color_gradientn(colors = c("#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a"), breaks = c(0,180,360,540,720,750))
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
                                 title='Papillomavirus Unweighted UniFrac PCoA Controls') +
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
                                 color = "DSC_all", # metadata variable
                                 axes = c(1,2),
                                 title='Papillomavirus Unweighted UniFrac PCoA Shedders') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_shape_manual(values = c(17,17))
PCoA_unw_shed

##PERMANOVA##
adonis.metadata <- read.delim("Virome_Metadata_papilloma.txt", header = TRUE, sep = '\t')

##All samples, weighted UniFrac, with interaction term.##
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = adonis.metadata$Patient_factor))
set.seed(100)
adonis_int <- adonis2(matrix ~ Days_since_visit1*Shedder_or_control + Days_since_visit1 + Shedder_or_control, data = adonis.metadata, permutations = perm, by = "margin")
adonis_int
write.table(adonis_int, "papilloma_weighted_UniFrac_adonis2_allSamples_interaction.txt", sep = "\t")
##No interaction term.##
set.seed(100)
adonis_noInt <- adonis2(matrix ~ Days_since_visit1 + Shedder_or_control, data = adonis.metadata, permutations = perm, by = "margin")
adonis_noInt
write.table(adonis_noInt, "papilloma_weighted_UniFrac_adonis2_allSamples_noInteraction.txt", sep = "\t")

##All samples, unweighted UniFrac, with interaction term.##
set.seed(100)
adonis_unw_int <- adonis2(matrix2 ~ Days_since_visit1*Shedder_or_control + Days_since_visit1 + Shedder_or_control, data = adonis.metadata, permutations = perm, by = "margin")
adonis_unw_int
write.table(adonis_unw_int, "papilloma_unweighted_UniFrac_adonis2_allSamples_interaction.txt", sep = "\t")
##No interaction term.##
set.seed(100)
adonis_unw <- adonis2(matrix2 ~ Days_since_visit1 + Shedder_or_control, data = adonis.metadata, permutations = perm, by = "margin")
adonis_unw
write.table(adonis_unw, "papilloma_unweighted_UniFrac_adonis2_allSamples_noInteraction.txt", sep = "\t")

metadata.ctrls <- adonis.metadata[adonis.metadata$Shedder_or_control=="Control",]

##Subset weighted UniFrac, non-shedders.##
ctrls <- matrix[,colnames(matrix) %in% metadata.ctrls$RCA_ID]
ctrls <- ctrls[rownames(ctrls) %in% metadata.ctrls$RCA_ID,]

##Non-shedders over time.##
set.seed(100)
perm.ctrls <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata.ctrls$Patient_factor))
set.seed(100)
adonis_ctrl <- adonis2(ctrls ~ Days_since_visit1, data = metadata.ctrls, permutations = perm.ctrls, by = "margin")
adonis_ctrl
write.table(adonis_ctrl, "papilloma_weighted_UniFrac_adonis2_controls.txt",sep = "\t")

metadata.shed <- adonis.metadata[adonis.metadata$Shedder_or_control=="Shedder",]

##Subset weighted UniFrac, non-shedders.##
shed <- matrix[,colnames(matrix) %in% metadata.shed$RCA_ID]
shed <- shed[rownames(shed) %in% metadata.shed$RCA_ID,]

##Shedders over time.##
set.seed(100)
perm.shed <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata.shed$Patient_factor))
set.seed(100)
adonis_shed <- adonis2(shed ~ Days_since_visit1, data = metadata.shed, permutations = perm.shed, by = "margin")
adonis_shed
write.table(adonis_shed, "papilloma_weighted_UniFrac_adonis2_shedders.txt",sep = "\t")

##Shedders, DSC vs non-DSC.##
set.seed(100)
adonis_DSC <- adonis2(shed ~ DSC_all, data = metadata.shed, permutations = perm.shed, by = "margin")
adonis_DSC
write.table(adonis_DSC, "papilloma_weighted_UniFrac_adonis2_shedders_DSC.txt",sep = "\t")

##Subset unweighted UniFrac, non-shedders.##
ctrls2 <- matrix2[,colnames(matrix2) %in% metadata.ctrls$RCA_ID]
ctrls2 <- ctrls2[rownames(ctrls2) %in% metadata.ctrls$RCA_ID,]

##Non-shedders over time, unweighted.##
set.seed(100)
adonis_ctrl2 <- adonis2(ctrls2 ~ Days_since_visit1, data = metadata.ctrls, permutations = perm.ctrls, by = "margin")
adonis_ctrl2
write.table(adonis_ctrl2, "papilloma_unweighted_UniFrac_adonis2_controls.txt",sep = "\t")

##Subset unweighted UniFrac, shedders.##
shed2 <- matrix2[,colnames(matrix2) %in% metadata.shed$RCA_ID]
shed2 <- shed2[rownames(shed2) %in% metadata.shed$RCA_ID,]

##Shedders over time, unweighted.##
set.seed(100)
adonis_shed2 <- adonis2(shed2 ~ Days_since_visit1, data = metadata.shed, permutations = perm.shed, by = "margin")
adonis_shed2
write.table(adonis_shed2, "papilloma_unweighted_UniFrac_adonis2_shedders.txt",sep = "\t")

##Shedders, DSC vs non-DSC, unweighted.##
set.seed(100)
adonis_DSC2 <- adonis2(shed2 ~ DSC_all, data = metadata.shed, permutations = perm.shed, by = "margin")
adonis_DSC2
write.table(adonis_DSC2, "papilloma_unweighted_UniFrac_adonis2_shedders_DSC.txt",sep = "\t")
