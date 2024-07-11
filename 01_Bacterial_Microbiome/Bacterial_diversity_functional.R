##Contamination analysis, pathways. Input: pathway abundance table from HUMAnN3, with abundance expressed as reads per kilobase (RPK).##
library(decontam)
pathabund <- read.delim("humann3_pathabundance_RPK_joined.tsv", row.names = 1, header = TRUE)
colnames(pathabund) <- sub("f.*","",colnames(pathabund))
metadata <- read.delim("metadata.txt", header = TRUE)
identical(colnames(pathabund),metadata$Concat_bacterial_ID)
unstratified <- pathabund[!grepl("g__.*s__|unclassified|UNMAPPED",rownames(pathabund)),]
unstratified <- as.matrix(t(unstratified))
contam<-isContaminant(unstratified, method = 'prevalence', neg =metadata$neg_ctrl, threshold=0.1) #No contaminants identified.
unstratified <- t(unstratified)
no_ctrls <- unstratified[,!metadata$any_ctrl=="TRUE"]
no_ctrls <- no_ctrls[rowSums(no_ctrls) > 0,]
no_ctrls <- as.data.frame(no_ctrls)
relab2 <- function(x){
  return(sweep(x,2,colSums(x),"/"))
}
no_zeroSums <- no_ctrls[,colSums(no_ctrls) > 0]
path_relab <- relab2(no_zeroSums)
write.table(path_relab, "humann_unstratified_pathway_relativeAbundance.txt", sep = "\t", row.names = TRUE, col.names = NA)

##Heatmap of pathway relative abundance. Input: clean unstratified pathway relative abundance from line 18.##
library(gplots)
library(RColorBrewer)
plot_integ <- as.matrix(path_relab[!rownames(path_relab) %in% "UNINTEGRATED",])
palette <- colorRampPalette(c("#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15"))(n = 10)
plot_prev <- plot_integ[rowSums(plot_integ>0)>12,] #Drop rows with 12 or fewer non-zero values (10% prevalence threshold).
plot_prev_na <- plot_prev
plot_prev_na[plot_prev_na == 0] <- NA
heatmap_prev <- heatmap.2(plot_prev, Rowv = TRUE, Colv = TRUE, trace = "none", col = palette)
heatmap_prev2 <- heatmap.2(plot_prev_na, Rowv = heatmap_prev$rowDendrogram, Colv = heatmap_prev$colDendrogram, trace = "none", col = palette, na.color = "#d3d3d3")

##Maaslin2, pathways. Input: clean unstratified pathway relative abundance from line 18.##
library(Maaslin2)
maaslin_met <-read.delim("metadata_humann_path.txt", row.names = 1)
maaslin_met$shedder_time = (maaslin_met$Shedder_or_control == "Shedder") * maaslin_met$Days_since_visit1
maaslin_met$control_time = (maaslin_met$Shedder_or_control == "Control") * maaslin_met$Days_since_visit1
  #Pathway relative abundance vs time and shedder/nonshedder status with interaction term:
interaction= Maaslin2(
  input_data = path_relab, 
  input_metadata = maaslin_met,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_interaction",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1", "Shedder_or_control","shedder_time"),
  random_effects = c("Patient_factor"))
interaction2= Maaslin2(
  input_data = path_relab, 
  input_metadata = maaslin_met,
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
  #Pathway relative abundance vs time and shedder/nonshedder status, no interaction term:
full_model= Maaslin2(
  input_data = path_relab, 
  input_metadata = maaslin_met,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_path_relab_shed_ctrl_time",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1", "Shedder_or_control"),
  random_effects = c("Patient_factor"))
full_model_results <- read.delim("maaslin2_path_relab_shed_ctrl_time/all_results.tsv")
full_model_time <- subset(full_model_results,metadata=="Days_since_visit1")
full_model_time$qval.readjust <- p.adjust(full_model_time$pval, method = "BH")
full_model_shed_ctrl <- subset(full_model_results,metadata=="Shedder_or_control")
full_model_shed_ctrl$qval.readjust <- p.adjust(full_model_shed_ctrl$pval, method = "BH")
  #Pathway relative abundance by shedding timepoint, shedders only:
shed <- subset(maaslin_met,Shedder_or_control == "Shedder")
path_relab_shed <- path_relab[,colnames(path_relab) %in% shed$I.number]
path_relab_DSC= Maaslin2(
  input_data = path_relab_shed, 
  input_metadata = shed,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_path_relab_DSC",
  transform = "none",
  normalization = "none",
  fixed_effects = c("DSC_all"),
  random_effects = c("Patient_factor"))
  #Pathway relative abundance by microbiome community cluster, shedders and nonshedders:
c2ref <- subset(maaslin_met,Cluster=="C2"|Cluster=="C3"|Cluster=="C4"|Cluster=="C5"|Cluster=="C6")
c3ref <- subset(maaslin_met,Cluster=="C3"|Cluster=="C4"|Cluster=="C5"|Cluster=="C6")
c4ref <- subset(maaslin_met,Cluster=="C4"|Cluster=="C5"|Cluster=="C6")
c5ref <- subset(maaslin_met,Cluster=="C5"|Cluster=="C6")
c2ref_df <- path_relab[,colnames(path_relab) %in% c2ref$I.number]
c3ref_df <- path_relab[,colnames(path_relab) %in% c3ref$I.number]
c4ref_df <- path_relab[,colnames(path_relab) %in% c4ref$I.number]
c5ref_df <- path_relab[,colnames(path_relab) %in% c5ref$I.number]
path_relab_c1ref= Maaslin2(
  input_data = path_relab, 
  input_metadata = maaslin_met,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_path_relab_c1ref",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Cluster"),
  reference = c("Cluster,C1"),
  random_effects = c("Patient_factor"))
path_relab_c2ref= Maaslin2(
  input_data = c2ref_df, 
  input_metadata = c2ref,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_path_relab_c2ref",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Cluster"),
  reference = c("Cluster,C2"),
  random_effects = c("Patient_factor"))
path_relab_c3ref= Maaslin2(
  input_data = c3ref_df, 
  input_metadata = c3ref,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_path_relab_c3ref",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Cluster"),
  reference = c("Cluster,C3"),
  random_effects = c("Patient_factor"))
path_relab_c4ref= Maaslin2(
  input_data = c4ref_df, 
  input_metadata = c4ref,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_path_relab_c4ref",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Cluster"),
  reference = c("Cluster,C4"),
  random_effects = c("Patient_factor"))
path_relab_c5ref= Maaslin2(
  input_data = c5ref_df, 
  input_metadata = c5ref,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_path_relab_c5ref",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Cluster"),
  reference = c("Cluster,C5"),
  random_effects = c("Patient_factor"))
c1results <- read.delim("maaslin2_path_relab_c1ref/all_results.tsv")
c1results$reference <- c("C1")
c2results <- read.delim("maaslin2_path_relab_c2ref/all_results.tsv")
c2results$reference <- c("C2")
c3results <- read.delim("maaslin2_path_relab_c3ref/all_results.tsv")
c3results$reference <- c("C3")
c4results <- read.delim("maaslin2_path_relab_c4ref/all_results.tsv")
c4results$reference <- c("C4")
c5results <- read.delim("maaslin2_path_relab_c5ref/all_results.tsv")
c5results$reference <- c("C5")
allresults <- rbind(c1results,c2results,c3results,c4results,c5results)
allresults$qval.readjust <- p.adjust(allresults$pval, method = "BH")

##Bray-Curtis dissimilarity, pathways. Input: clean unstratified pathway relative abundance from line 18.##
library(vegan)
path_relab_t <- t(path_relab)
dis.bc <- vegdist(path_relab_t, method = "bray")

##Permutation tests of pathway Bray-Curtis dissimilarity within and between groups. Input: table with Bray-Curtis dissimilarity and metadata for sample pairs.##
library(coin)
stats.data <- read.delim("path_relab_BC_long_format_metadata.txt", header = TRUE, sep = "\t")
stats.data$comparison_all <- as.factor(stats.data$comparison_all)
stats.data$comparison_DSC_vs_notDSC <- as.factor(stats.data$comparison_DSC_vs_notDSC)
stats.data$block <- c("block1")
stats.data$block <- as.factor(stats.data$block)
w_b.ctrls<-subset(stats.data,comparison_all =='Control.Control.within'|comparison_all =='Control.Control.between')
w_b.shed<-subset(stats.data,comparison_all =='Shedder.Shedder.within'|comparison_all =='Shedder.Shedder.between')
w.ctrl_shed<-subset(stats.data,comparison_all =='Control.Control.within'|comparison_all =='Shedder.Shedder.within')
b.ctrl_shed<-subset(stats.data,comparison_all =='Control.Control.between'|comparison_all =='Shedder.Shedder.between')
w.shed.DSC<-subset(stats.data,comparison_DSC_vs_notDSC =='Shedder.within.nonDSC.nonDSC'|comparison_DSC_vs_notDSC =='Shedder.within.DSC.nonDSC')
w.ctrl.DSC<-subset(stats.data,comparison_DSC_vs_notDSC =='Control.within.nonDSC.nonDSC'|comparison_DSC_vs_notDSC =='Control.within.DSC.nonDSC')
set.seed(100)
perm.w_b.ctrls<-independence_test(bray_curtis~comparison_all|block, data = w_b.ctrls, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w_b.shed<-independence_test(bray_curtis~comparison_all|block, data = w_b.shed, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w.ctrl_shed<-independence_test(bray_curtis~comparison_all|block, data = w.ctrl_shed, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.b.ctrl_shed<-independence_test(bray_curtis~comparison_all|block, data = b.ctrl_shed, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w.shed.DSC<-independence_test(bray_curtis~comparison_DSC_vs_notDSC|block, data = w.shed.DSC, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w.ctrl.DSC<-independence_test(bray_curtis~comparison_DSC_vs_notDSC|block, data = w.ctrl.DSC, alternative = c("two.sided"), distribution=approximate(nresample=100000))
bc.test<-c("perm.w_b.ctrls","perm.w_b.shed","perm.w.ctrl_shed","perm.b.ctrl_shed","perm.w.shed.DSC","perm.w.ctrl.DSC")
p.val<-c(pvalue(perm.w_b.ctrls),pvalue(perm.w_b.shed),pvalue(perm.w.ctrl_shed),pvalue(perm.b.ctrl_shed),pvalue(perm.w.shed.DSC),pvalue(perm.w.ctrl.DSC))
bc.pvals<-data.frame(bc.test,p.val)
##Benjamini-Hochberg correction for multiple comparisons.##
pval.adjust <- p.adjust(bc.pvals$p.val, method = "BH")
bc.adjusted<-cbind(bc.pvals,pval.adjust)

##PCoA, pathways. Input: clean unstratified pathway relative abundance.##
library(phyloseq)
library(ggplot2)
path_relab <- read.delim("humann_unstratified_pathway_relativeAbundance.txt", sep = "\t", row.names = 1)
OTU=otu_table(path_relab,taxa_are_rows =TRUE)
sample_names(OTU)
metadata2<-import_qiime_sample_data("Master_Metadata_humann_path.txt")
sample_names(metadata2)
identical(sample_names(OTU),sample_names(metadata2))
physeq =             phyloseq(OTU,metadata2)
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
  scale_shape_manual(values = c(17, 17))
PCoA.BC.shed

##PERMANOVA, pathways. Input: Bray-Curtis dissimilarity matrix from line 173.##
library(vegan)
metadata3 <- read.delim("metadata_humann_path.txt", header = TRUE, sep = '\t')
dis.bc3<-as.matrix(dis.bc)
#By time and shedder/non-shedder status, with interaction term:
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata3$Patient_factor))
set.seed(100)
adonis_int <- adonis2(dis.bc ~ Days_since_visit1*Shedder_or_control + Days_since_visit1 + Shedder_or_control, data = metadata3, permutations = perm, by = "margin")
adonis_int
#By time and shedder/non-shedder status, no interaction term.##
set.seed(100)
adonis_noInt <- adonis2(dis.bc ~ Days_since_visit1 + Shedder_or_control, data = metadata3, permutations = perm, by = "margin")
adonis_noInt
#By shedding timepoint, shedders only:
sub_metadata <- metadata3[metadata3$Shedder_or_control=="Shedder",]
sub_BC <- dis.bc3[,colnames(dis.bc3) %in% sub_metadata$Concat_bacterial_ID]
sub_BC2 <- sub_BC[rownames(sub_BC) %in% sub_metadata$Concat_bacterial_ID,]
dis.bc.shed <- as.dist(sub_BC2)
set.seed(100)
perm2 <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_metadata$Patient_factor))
set.seed(100)
adonis_DSC <- adonis2(dis.bc.shed ~ DSC_all, data = sub_metadata, permutations = perm2, by = "margin")
adonis_DSC

##Contamination analysis, gene families. Input: gene family table from HUMAnN3, with abundance expressed as reads per kilobase (RPK).##
library(decontam)
genes <- read.delim("humann3_genefamilies_RPK_joined.tsv", row.names = 1, header = TRUE)
colnames(genes) <- sub("f.*","",colnames(genes))
metadata <- read.delim("metadata.txt", header = TRUE)
unstrat_genes <- genes[!grepl("g__.*s__|unclassified|UNMAPPED",rownames(genes)),]
unstrat_genes <- as.matrix(t(unstrat_genes))
contam_genes<-isContaminant(unstrat_genes, method = 'prevalence', neg =metadata$neg_ctrl, threshold=0.1)
clean <- t(unstrat_genes[,!contam_genes$contaminant=="TRUE"])
no_ctrls_genes <- clean[,!metadata$any_ctrl=="TRUE"]
no_ctrls_genes <- no_ctrls_genes[rowSums(no_ctrls_genes) > 0,]
no_ctrls_genes <- as.data.frame(no_ctrls_genes)
relab2 <- function(x){
  return(sweep(x,2,colSums(x),"/"))
}
no_zeroSums_genes <- no_ctrls_genes[,colSums(no_ctrls_genes) > 0]
gene_relab <- relab2(no_zeroSums_genes)
write.table(gene_relab, "humann_clean_unstratified_genefam_relativeAbundance.txt", sep = "\t", row.names = TRUE, col.names = NA)

##MaAsLin2, gene families. Input: clean gene family relative abundance from line 280.##
library(Maaslin2)
maaslin_met <-read.delim("Master_Metadata_021023.txt", row.names = 1)
maaslin_met$shedder_time = (maaslin_met$Shedder_or_control == "Shedder") * maaslin_met$Days_since_visit1
maaslin_met$control_time = (maaslin_met$Shedder_or_control == "Control") * maaslin_met$Days_since_visit1
  #Gene family relative abundance vs time and shedder/nonshedder status, with interaction term:
interact_gene= Maaslin2(
  input_data = gene_relab, 
  input_metadata = maaslin_met,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_gene_interaction",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1", "Shedder_or_control","shedder_time"),
  random_effects = c("Patient_factor"),
  plot_scatter = FALSE)
interact_gene2= Maaslin2(
  input_data = gene_relab, 
  input_metadata = maaslin_met,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_gene_interaction2",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1", "Shedder_or_control","control_time"),
  random_effects = c("Patient_factor"),
  plot_scatter = FALSE)
interaction_all_results <- read.delim("maaslin2_gene_interaction/all_results.tsv")
interaction_subset <- subset(interaction_all_results,metadata=="shedder_time")
interaction_subset$qval.readjust <- p.adjust(interaction_subset$pval, method = "BH")
interaction2_all_results <- read.delim("maaslin2_gene_interaction2/all_results.tsv")
interaction2_subset <- subset(interaction2_all_results,metadata=="control_time")
interaction2_subset$qval.readjust <- p.adjust(interaction2_subset$pval, method = "BH")
  #By time and shedder/non-shedder status, no interaction term:
full_model_gene= Maaslin2(
  input_data = gene_relab, 
  input_metadata = maaslin_met,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_gene_relab_shed_ctrl_time",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Days_since_visit1", "Shedder_or_control"),
  random_effects = c("Patient_factor"),
  plot_scatter = FALSE)
full_model_gene_results <- read.delim("maaslin2_gene_relab_shed_ctrl_time/all_results.tsv")
full_model_gene_time <- subset(full_model_gene_results,metadata=="Days_since_visit1")
full_model_gene_time$qval.readjust <- p.adjust(full_model_gene_time$pval, method = "BH")
full_model_gene_shed_ctrl <- subset(full_model_gene_results,metadata=="Shedder_or_control")
full_model_gene_shed_ctrl$qval.readjust <- p.adjust(full_model_gene_shed_ctrl$pval, method = "BH")
  #By shedding timepoint, shedders only:
shed <- subset(maaslin_met,Shedder_or_control == "Shedder")
shed$I.number <- rownames(shed)
gene_relab_shed <- gene_relab[,colnames(gene_relab) %in% shed$I.number]
gene_relab_DSC= Maaslin2(
  input_data = gene_relab_shed, 
  input_metadata = shed,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_gene_relab_DSC",
  transform = "none",
  normalization = "none",
  fixed_effects = c("DSC_all"),
  random_effects = c("Patient_factor"),
  plot_scatter = FALSE)
  #Gene family relative abundance by microbiome community cluster, shedders and nonshedders:
maaslin_met$I.number <- rownames(maaslin_met)
c2ref <- subset(maaslin_met,Cluster=="C2"|Cluster=="C3"|Cluster=="C4"|Cluster=="C5"|Cluster=="C6")
c3ref <- subset(maaslin_met,Cluster=="C3"|Cluster=="C4"|Cluster=="C5"|Cluster=="C6")
c4ref <- subset(maaslin_met,Cluster=="C4"|Cluster=="C5"|Cluster=="C6")
c5ref <- subset(maaslin_met,Cluster=="C5"|Cluster=="C6")
c2ref_df_gene <- gene_relab[,colnames(gene_relab) %in% c2ref$I.number]
c3ref_df_gene <- gene_relab[,colnames(gene_relab) %in% c3ref$I.number]
c4ref_df_gene <- gene_relab[,colnames(gene_relab) %in% c4ref$I.number]
c5ref_df_gene <- gene_relab[,colnames(gene_relab) %in% c5ref$I.number]
gene_relab_c1ref= Maaslin2(
  input_data = gene_relab, 
  input_metadata = maaslin_met,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_gene_relab_c1ref",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Cluster"),
  reference = c("Cluster,C1"),
  random_effects = c("Patient_factor"),
  plot_scatter = FALSE)
gene_relab_c2ref= Maaslin2(
  input_data = c2ref_df_gene, 
  input_metadata = c2ref,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_gene_relab_c2ref",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Cluster"),
  reference = c("Cluster,C2"),
  random_effects = c("Patient_factor"),
  plot_scatter = FALSE)
gene_relab_c3ref= Maaslin2(
  input_data = c3ref_df_gene, 
  input_metadata = c3ref,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_gene_relab_c3ref",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Cluster"),
  reference = c("Cluster,C3"),
  random_effects = c("Patient_factor"),
  plot_scatter = FALSE)
gene_relab_c4ref= Maaslin2(
  input_data = c4ref_df_gene, 
  input_metadata = c4ref,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_gene_relab_c4ref",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Cluster"),
  reference = c("Cluster,C4"),
  random_effects = c("Patient_factor"),
  plot_scatter = FALSE)
gene_relab_c5ref= Maaslin2(
  input_data = c5ref_df_gene, 
  input_metadata = c5ref,
  min_abundance = 0,
  min_prevalence = 0.1,
  output = "maaslin2_gene_relab_c5ref",
  transform = "none",
  normalization = "none",
  fixed_effects = c("Cluster"),
  reference = c("Cluster,C5"),
  random_effects = c("Patient_factor"),
  plot_scatter = FALSE)
c1results_gene <- read.delim("maaslin2_gene_relab_c1ref/all_results.tsv")
c1results_gene$reference <- c("C1")
c2results_gene <- read.delim("maaslin2_gene_relab_c2ref/all_results.tsv")
c2results_gene$reference <- c("C2")
c3results_gene <- read.delim("maaslin2_gene_relab_c3ref/all_results.tsv")
c3results_gene$reference <- c("C3")
c4results_gene <- read.delim("maaslin2_gene_relab_c4ref/all_results.tsv")
c4results_gene$reference <- c("C4")
c5results_gene <- read.delim("maaslin2_gene_relab_c5ref/all_results.tsv")
c5results_gene$reference <- c("C5")
allresults_gene <- rbind(c1results_gene,c2results_gene,c3results_gene,c4results_gene,c5results_gene)
allresults_gene$qval.readjust <- p.adjust(allresults_gene$pval, method = "BH")

##Bray-Curtis dissimilarity, gene families. Input: clean unstratified gene family relative abundance.##
library(vegan)
gene_relab <- read.delim("humann_clean_unstratified_genefam_relativeAbundance.txt", row.names = 1)
gene_relab_t <- t(gene_relab)
dis_bc_gene <- vegdist(gene_relab_t, method = "bray")

##Permutation tests of gene family Bray-Curtis dissimilarity within and between groups. Input: table with Bray-Curtis dissimilarity and metadata for sample pairs.##
library(coin)
stats.data <- read.delim("gene_relab_BC_long_format_metadata.txt", header = TRUE, sep = "\t")
stats.data$comparison_all <- as.factor(stats.data$comparison_all)
stats.data$comparison_DSC_vs_notDSC <- as.factor(stats.data$comparison_DSC_vs_notDSC)
stats.data$block <- c("block1")
stats.data$block <- as.factor(stats.data$block)
w_b.ctrls<-subset(stats.data,comparison_all =='Control.within'|comparison_all =='Control.between')
w_b.shed<-subset(stats.data,comparison_all =='Shedder.within'|comparison_all =='Shedder.between')
w.ctrl_shed<-subset(stats.data,comparison_all =='Control.within'|comparison_all =='Shedder.within')
b.ctrl_shed<-subset(stats.data,comparison_all =='Control.between'|comparison_all =='Shedder.between')
w.shed.DSC<-subset(stats.data,comparison_DSC_vs_notDSC =='Shedder.nonDSC.nonDSC'|comparison_DSC_vs_notDSC =='Shedder.DSC.nonDSC')
w.ctrl.DSC<-subset(stats.data,comparison_DSC_vs_notDSC =='Control.nonDSC.nonDSC'|comparison_DSC_vs_notDSC =='Control.DSC.nonDSC')
set.seed(100)
perm.w_b.ctrls<-independence_test(bray_curtis~comparison_all|block, data = w_b.ctrls, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w_b.shed<-independence_test(bray_curtis~comparison_all|block, data = w_b.shed, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w.ctrl_shed<-independence_test(bray_curtis~comparison_all|block, data = w.ctrl_shed, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.b.ctrl_shed<-independence_test(bray_curtis~comparison_all|block, data = b.ctrl_shed, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w.shed.DSC<-independence_test(bray_curtis~comparison_DSC_vs_notDSC|block, data = w.shed.DSC, alternative = c("two.sided"), distribution=approximate(nresample=100000))
perm.w.ctrl.DSC<-independence_test(bray_curtis~comparison_DSC_vs_notDSC|block, data = w.ctrl.DSC, alternative = c("two.sided"), distribution=approximate(nresample=100000))
bc.test<-c("perm.w_b.ctrls","perm.w_b.shed","perm.w.ctrl_shed","perm.b.ctrl_shed","perm.w.shed.DSC","perm.w.ctrl.DSC")
p.val<-c(pvalue(perm.w_b.ctrls),pvalue(perm.w_b.shed),pvalue(perm.w.ctrl_shed),pvalue(perm.b.ctrl_shed),pvalue(perm.w.shed.DSC),pvalue(perm.w.ctrl.DSC))
bc.pvals<-data.frame(bc.test,p.val)
##Benjamini-Hochberg correction for multiple comparisons.##
pval.adjust <- p.adjust(bc.pvals$p.val, method = "BH")
bc.adjusted<-cbind(bc.pvals,pval.adjust)

##PCoA, gene families. Input: clean unstratified gene family relative abundance.##
library(phyloseq)
library(ggplot2)
OTU_gene=otu_table(gene_relab,taxa_are_rows =TRUE)
sample_names(OTU_gene)
metadata2<-import_qiime_sample_data("metadata.txt")
sample_names(metadata2)
identical(sample_names(OTU_gene),sample_names(metadata2))
physeq_gene =             phyloseq(OTU_gene,metadata2)
  #All samples:
set.seed(100)
ord_gene <- ordinate(physeq_gene,
                method = "PCoA",
                distance = "bray")

PCoA.BC.gene  <- plot_ordination(physeq = physeq_gene,
                            ordination = ord_gene,
                            shape = "Shedder_or_control",
                            color = "Days_since_visit1",
                            axes = c(1,2),
                            title='Bray Curtis PCoA') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_color_gradientn(colors = c("#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a"), breaks = c(0,180,360,540,720,750))
PCoA.BC.gene
  #Shedders only:
physeq_shed_gene <- subset_samples(physeq_gene, Shedder_or_control=="Shedder")
set.seed(100)
ord_shed_gene <- ordinate(physeq_shed_gene,
                     method = "PCoA",
                     distance = "bray")
PCoA.BC.shed.gene  <- plot_ordination(physeq = physeq_shed_gene,
                                 ordination = ord_shed_gene,
                                 shape = "DSC_all",
                                 color = "DSC_all",
                                 axes = c(1,2),
                                 title='Bray Curtis PCoA') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 9)+
  scale_shape_manual(values = c(17,17))
PCoA.BC.shed.gene

##PERMANOVA, gene families. Input: Bray-Curtis dissimilarity matrix from line 436.##
metadata3 <- read.delim("metadata.txt", header = TRUE, sep = '\t')
dis.bc3.gene<-as.matrix(dis_bc_gene)
  #By time and shedder/non-shedder status, with interaction term:
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata3$Patient_factor))
set.seed(100)
adonis_int <- adonis2(dis_bc_gene ~ Days_since_visit1*Shedder_or_control + Days_since_visit1 + Shedder_or_control, data = metadata3, permutations = perm, by = "margin")
adonis_int
  #By time and shedder/non-shedder status, no interaction term:
set.seed(100)
adonis_noInt <- adonis2(dis_bc_gene ~ Days_since_visit1 + Shedder_or_control, data = metadata3, permutations = perm, by = "margin")
adonis_noInt
  #By shedding timepoint, shedders only:
sub_metadata <- metadata3[metadata3$Shedder_or_control=="Shedder",]
sub_BC_gene <- dis.bc3.gene[,colnames(dis.bc3.gene) %in% sub_metadata$Concat_bacterial_ID]
sub_BC2_gene <- sub_BC_gene[rownames(sub_BC_gene) %in% sub_metadata$Concat_bacterial_ID,]
dis.bc.shed.gene <- as.dist(sub_BC2_gene)
identical(sub_metadata$Concat_bacterial_ID, rownames(sub_BC2))
set.seed(100)
perm2 <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_metadata$Patient_factor))
set.seed(100)
adonis_DSC_gene <- adonis2(dis.bc.shed.gene ~ DSC_all, data = sub_metadata, permutations = perm2, by = "margin")
adonis_DSC_gene

