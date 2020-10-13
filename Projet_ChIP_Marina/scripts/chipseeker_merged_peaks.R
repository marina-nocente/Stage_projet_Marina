##### Script Marina pour l'annotation des peaks avec ChIPseeker #####


# https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/12_annotation_functional_analysis.html

# http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#abstract


#### Intallation du package ChIPseeker et de l'annotation des genes #####

# via un environnement conda et un fichier de config "config_R.yaml" pour gerer les versions et les packages.


##### Importation des libraries: #####
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # Annotation package for TxDb object(s)
library(clusterProfiler)
library(annotables)
library(org.Mm.eg.db) # Genome wide annotation for Mouse

print("les libraries ont bien été chargées")


##### Loading des data #####
# As input we need to provide the names of our BED files in a list format.

### les replicats de Pol2, Chd8 et Oct4 ont été mergés. On n'avait malheureusement pas de duplicat pour CTCF et le réplicats 2 de TBP était très mauvais donc j'ai gardé tous les peaks pour CTCF et uniquement le réplicat 1 pour TBP. ###
samplefiles <- list.files("/shared/home/mnocente/macs2/merge/final/modif_final", pattern= ".narrowPeak",full.names = T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("Chd8", "Oct4", "Pol2", "CTCF", "TBP_rep1")


### TBP rep1 uniquement ###
samplefiles_TBP <- list.files("/home/mnocente/Bureau/", pattern= ".narrowPeak",full.names = T)
samplefiles_TBP <- as.list(samplefiles_TBP)
names(samplefiles_TBP) <- c("TBP_rep1")


### TBP_rep1_Pol2_Chd8 uniquement ###
samplefiles_TBP_rep1_Pol2_Chd8 <- list.files("/home/mnocente/Bureau/", pattern= ".narrowPeak",full.names = T)
samplefiles_TBP_rep1_Pol2_Chd8 <- as.list(samplefiles_TBP_rep1_Pol2_Chd8)
names(samplefiles_TBP_rep1_Pol2_Chd8) <- c("TBP_rep1_Pol2_Chd8")

print("les fichiers ont bien été chargés")
print(samplefiles_TBP_rep1_Pol2_Chd8)

### Les intersections de facteurs ###
samplefiles_intersections <- list.files("/shared/home/mnocente/macs2/merge/final/modif_final/intersection", pattern= ".narrowPeak",full.names = T)
samplefiles_intersections <- as.list(samplefiles_intersections)
names(samplefiles_intersections) <- c("Chd8_Oct4_intersection", "Oct4_CTCF_intersection", "TBP_rep1_Pol2_Chd8_intersection")

##### Assign annotation db #####

# We need to assign annotation databases generated from UCSC to a variable:

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
print("les annotations ont bien été chargés")

# ChIPseeker et sa fonction annotatePeak permettent d'annoter les pics avec le gène et la région génomique les plus proches de là où se trouve le pic.
# The annotatePeak function by default uses the TSS method, and provides parameters to specify a max distance cutoff.


##### Get annotations #####
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)

print(peakAnnoList) # annotation information is stored in the peakAnnoList

#########
peakAnnoList_TBP <- lapply(samplefiles_TBP, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)

print(peakAnnoList_TBP)

#########
peakAnnoList_TBP_rep1_Pol2_Chd8 <- lapply(samplefiles_TBP_rep1_Pol2_Chd8, annotatePeak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)

print(peakAnnoList_TBP_rep1_Pol2_Chd8)

#########
peakAnnoList_intersection <- lapply(samplefiles_intersections, annotatePeak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)
print(peakAnnoList_intersection)



##### Visualisation des annotations #####
############## Visualisation des annotations pour les échantillons CTCF, TBP_rep1, Pol2, Chd8 et Oct4 ##############
## Barchart (multiple samples for comparison)
png(file = "Barchart_comparison_sample_annotation_samples_mergePol2_mergeChd8_mergeOct4_CTCF_TBPrep1.png")
plotAnnoBar(peakAnnoList)
dev.off()

## Distribution of TF-binding loci relative to TSS
png(file = "Distribution_of_TF-binding_loci_relative_to_TSS_samples_mergePol2_mergeChd8_mergeOct4_CTCF_TBPrep1.png")
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")
dev.off()


############## Visualisation des annotations pour les intersections des échantillons ##############
## Barchart (multiple samples for comparison)
png(file = "Barchart_comparison_sample_annotation_samples_intersections.png")
plotAnnoBar(peakAnnoList_intersection)
dev.off()

## Distribution of TF-binding loci relative to TSS
png(file = "Distribution_of_TF-binding_loci_relative_to_TSS_samples_intersections.png")
plotDistToTSS(peakAnnoList_intersection, title="Distribution of transcription factor-binding loci \n relative to TSS")
dev.off()



### Get annotation data frame
library(dplyr)
Oct4_annot <- as.data.frame(peakAnnoList[["Oct4"]]@anno)
TBP_rep1_annot <- as.data.frame(peakAnnoList[["TBP_rep1"]]@anno)
TBP_annot <- as.data.frame(peakAnnoList[["TBP"]]@anno)
CTCF_annot <- as.data.frame(peakAnnoList[["CTCF"]]@anno)
Pol2_annot <- as.data.frame(peakAnnoList[["Pol2"]]@anno)
Chd8_annot <- as.data.frame(peakAnnoList[["Chd8"]]@anno)

Chd8_Oct4_intersection_annot <- as.data.frame(peakAnnoList_intersection[["Chd8_Oct4_intersection"]]@anno)
TBP_rep1_Pol2_Chd8_intersection_annot <- as.data.frame(peakAnnoList_intersection[["TBP_rep1_Pol2_Chd8_intersection"]]@anno)
Oct4_CTCF_intersection_annot <- as.data.frame(peakAnnoList_intersection[["Oct4_CTCF_intersection"]]@anno)


# dans ce dataframe, la derniere colonne correspond à la distance au TSS.
# Dans ces dataframes, on doit voir des colonnes correspondant au fichier BED d'entrée et des colonnes supplémentaires contenant le ou les gènes les plus proches, la distance entre le pic et le TSS du gène le plus proche, la région génomique du pic et d'autres informations. 
# Ordre de priorite : Promoter, 5’ UTR, 3’ UTR, Exon, Intron, Downstream (defined as the downstream of gene end), Intergenic.
# Il manque juste les symboles de gènes répertoriés dans le tableau.



###############################################################
##### Functional enrichment analysis #####
###############################################################

#### Single sample analysis ####

library(dplyr)

### La liste de gènes des annotations Oct4 en entree pour une analyse d'enrichissement GO.

## Run GO enrichment analysis 
entrezids_Oct4 <- Oct4_annot$geneId %>% as.character() %>% unique()

ego_Oct4 <- enrichGO(gene = entrezids_Oct4, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

## Output results from GO analysis to a table
#cluster_summary_Oct4 <- data.frame(ego_Oct4)
#write.csv(cluster_summary_Oct4, "results/clusterProfiler_Oct4.csv")

## Dotplot visualization
png(file = "dotplotfunctional_enrichment_Oct4.png", width = 750, height = 1000)
dotplot(ego_Oct4, showCategory=20)
dev.off()



### La liste de gènes des annotations CTCF en entree pour une analyse d'enrichissement GO.

## Run GO enrichment analysis 
entrezids_CTCF <- CTCF_annot$geneId %>% as.character() %>% unique()

ego_CTCF <- enrichGO(gene = entrezids_CTCF, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

## Output results from GO analysis to a table
#cluster_summary_CTCF <- data.frame(ego_CTCF)
#write.csv(cluster_summary_CTCF, "results/clusterProfiler_CTCF.csv")

## Dotplot visualization
png(file = "dotplotfunctional_enrichment_CTCF.png", width = 750, height = 1000)
dotplot(ego_CTCF, showCategory=20)
dev.off()



### La liste de gènes des annotations TBP en entree pour une analyse d'enrichissement GO.

## Run GO enrichment analysis 
entrezids_TBP_rep1 <- TBP_rep1_annot$geneId %>% as.character() %>% unique()

ego_TBP_rep1 <- enrichGO(gene = entrezids_TBP_rep1, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

## Output results from GO analysis to a table
#cluster_summary_TBP_rep1 <- data.frame(ego_TBP_rep1)
#write.csv(cluster_summary_TBP_rep1, "results/clusterProfiler_TBP.csv")

## Dotplot visualization
png(file = "dotplotfunctional_enrichment_TBP_rep1.png", width = 750, height = 1000)
dotplot(ego_TBP_rep1, showCategory=20)
dev.off()



### La liste de gènes des annotations Pol2 en entree pour une analyse d'enrichissement GO.

## Run GO enrichment analysis 
entrezids_Pol2 <- Pol2_annot$geneId %>% as.character() %>% unique()

ego_Pol2 <- enrichGO(gene = entrezids_Pol2, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

## Output results from GO analysis to a table
#cluster_summary_Pol2 <- data.frame(ego_Pol2)
#write.csv(cluster_summary_Pol2, "results/clusterProfiler_Pol2.csv")

## Dotplot visualization
png(file = "dotplotfunctional_enrichment_Pol2.png", width = 750, height = 1000)
dotplot(ego_Pol2, showCategory=20)
dev.off()



### La liste de gènes des annotations Chd8 en entree pour une analyse d'enrichissement GO.

## Run GO enrichment analysis 
entrezids_Chd8 <- Chd8_annot$geneId %>% as.character() %>% unique()

ego_Chd8 <- enrichGO(gene = entrezids_Chd8, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

## Output results from GO analysis to a table
#cluster_summary_Chd8 <- data.frame(ego_Chd8)
#write.csv(cluster_summary_Chd8, "results/clusterProfiler_Chd8.csv")

## Dotplot visualization
png(file = "dotplotfunctional_enrichment_Chd8.png", width = 750, height = 1000)
dotplot(ego_Chd8, showCategory=20)
dev.off()




### La liste de gènes des annotations Chd8_Oct4_merge en entree pour une analyse d'enrichissement GO.

## Run GO enrichment analysis 
entrezids_Chd8_Oct4_intersection <- Chd8_Oct4_intersection_annot$geneId %>% as.character() %>% unique()

ego_Chd8_Oct4_intersection <- enrichGO(gene = entrezids_Chd8_Oct4_intersection, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)


## Dotplot visualization
png(file = "dotplotfunctional_enrichment_Chd8_Oct4_intersection.png", width = 750, height = 1000)
dotplot(ego_Chd8_Oct4_intersection, showCategory=20)
dev.off()


### La liste de gènes des annotations TBP_rep1_Pol2_Chd8 en entree pour une analyse d'enrichissement GO.

## Run GO enrichment analysis 
entrezids_TBP_rep1_Pol2_Chd8_intersection <- TBP_rep1_Pol2_Chd8_intersection_annot$geneId %>% as.character() %>% unique()

ego_TBP_rep1_Pol2_Chd8_intersection <- enrichGO(gene = entrezids_TBP_rep1_Pol2_Chd8_intersection, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)


## Dotplot visualization
png(file = "dotplotfunctional_enrichment_TBP_rep1_Pol2_Chd8_intersection.png", width = 750, height = 1000)
dotplot(ego_TBP_rep1_Pol2_Chd8_intersection, showCategory=20)
dev.off()


### La liste de gènes des annotations Oct4_CTCF_merge en entree pour une analyse d'enrichissement GO.

## Run GO enrichment analysis 
entrezids_Oct4_CTCF_intersection <- Oct4_CTCF_intersection_annot$geneId %>% as.character() %>% unique()

ego_Oct4_CTCF_intersection <- enrichGO(gene = entrezids_Oct4_CTCF_intersection, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)


## Dotplot visualization
png(file = "dotplotfunctional_enrichment_Oct4_CTCF_intersection.png", width = 750, height = 1000)
dotplot(ego_Oct4_CTCF_intersection, showCategory=20)
dev.off()


### Multiple samples ###

# Notre ensemble de données se compose de peaks de plusieurs facteurs différents, il serait donc utile de comparer les résultats d'enrichissement fonctionnel.
## Create a list with genes from each sample
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

## Run KEGG analysis
compKEGG <- compareCluster(geneCluster = genes, 
                           fun = "enrichKEGG",
                           organism = "mouse",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")

png(file = "KEGG_Pathway_Enrichment_Analysis_samples_mergePol2_mergeChd8_mergeOct4_CTCF_TBPrep1.png", width = 1100, height = 1100)
dotplot(compKEGG, showCategory = 20, title = "KEGG_Pathway_Enrichment_Analysis")
dev.off()


### Intersections:
genes_intersections = lapply(peakAnnoList_intersection, function(i) as.data.frame(i)$geneId)

## Run KEGG analysis
compKEGG_intersections <- compareCluster(geneCluster = genes_intersections, 
                           fun = "enrichKEGG",
                           organism = "mouse",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")

png(file = "KEGG_Pathway_Enrichment_Analysis_samples_intersections.png", width = 1100, height = 1100)
dotplot(compKEGG_intersections, showCategory = 20, title = "KEGG_Pathway_Enrichment_Analysis")
dev.off()


#### Overlap of peaks and annotated genes ####
png(file = "Vennplot_mergePol2_mergeChd8_mergeOct4_CTCF_TBPrep1.png", width = 500, height = 500)
genesOverlap= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genesOverlap)
dev.off()

png(file = "Vennplot_Oct4_CTCF.png", width = 500, height = 500)
peakAnnoList_Oct4_CTCF <- c(peakAnnoList$Oct4, peakAnnoList$CTCF)
names(peakAnnoList_Oct4_CTCF) <- c("Oct4", "CTCF")
peakAnnoList_Oct4_CTCF
genesOverlap_Oct4_CTCF= lapply(peakAnnoList_Oct4_CTCF, function(i) as.data.frame(i)$geneId)
vennplot(genesOverlap_Oct4_CTCF)
dev.off()


png(file = "Vennplot_TBP_rep1_Pol2_Chd8.png", width = 500, height = 500)
peakAnnoList_TBP_rep1_Pol2_Chd8 <- c(peakAnnoList$Chd8, peakAnnoList$TBP_rep1, peakAnnoList$Pol2)
names(peakAnnoList_TBP_rep1_Pol2_Chd8) <- c("Chd8", "TBP", "Pol2")
peakAnnoList_TBP_rep1_Pol2_Chd8
genesOverlap_TBP_rep1_Pol2_Chd8= lapply(peakAnnoList_TBP_rep1_Pol2_Chd8, function(i) as.data.frame(i)$geneId)
vennplot(genesOverlap_TBP_rep1_Pol2_Chd8)
dev.off()

png(file = "Vennplot_Chd8_Oct4.png", width = 500, height = 500)
peakAnnoList_Chd8_Oct4 <- c(peakAnnoList$Chd8, peakAnnoList$Oct4)
names(peakAnnoList_Chd8_Oct4) <- c("Chd8", "Oct4")
peakAnnoList_Chd8_Oct4
genesOverlap_Chd8_Oct4= lapply(peakAnnoList_Chd8_Oct4, function(i) as.data.frame(i)$geneId)
vennplot(genesOverlap_Chd8_Oct4)
dev.off()


######################################################################
#### Comptage du nombre de bases sous les peaks de chaque facteur ####
######################################################################

#### Oct4
Oct4_avant2000 <- subset(Oct4_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
Oct4_apres2000 <- subset(Oct4_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
Oct4_distal <- rbind(Oct4_avant2000, Oct4_apres2000) 
head(Oct4_distal)

Oct4_prox <- subset(Oct4_annot, distanceToTSS >= -400) # on garde les peaks entre -400 pb et +100 pb autour du TSS
Oct4_prox <- subset(Oct4_prox, distanceToTSS <= 100) 
head(Oct4_prox)

nb_pb_peaks_Oct_distal <- sum(Oct4_distal$width)
nb_pb_peaks_Oct_distal

nb_pb_peaks_Oct_prox <- sum(Oct4_prox$width)
nb_pb_peaks_Oct_prox

nb_pb_peaks_Oct_total <- sum(Oct4_annot$width)
nb_pb_peaks_Oct_total

nb_pb_peaks_Oct_non_distal <- nb_pb_peaks_Oct_total - nb_pb_peaks_Oct_distal
nb_pb_peaks_Oct_non_distal
nb_pb_peaks_Oct_non_prox <- nb_pb_peaks_Oct_total - nb_pb_peaks_Oct_prox
nb_pb_peaks_Oct_non_prox

# On exporte les fichiers contenant les peaks d'Oct4 distaux et proximaux
write.table(x = Oct4_distal, file = "peaks_Oct4_distaux.tsv", sep = "\t")
write.table(x = Oct4_prox, file = "peaks_Oct4_proximaux.tsv", sep = "\t")



#### TBP_rep1
TBP_rep1_avant2000 <- subset(TBP_rep1_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
TBP_rep1_apres2000 <- subset(TBP_rep1_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
TBP_rep1_distal <- rbind(TBP_rep1_avant2000, TBP_rep1_apres2000) 
head(TBP_rep1_distal)

TBP_rep1_prox <- subset(TBP_rep1_annot, distanceToTSS >= -400) # on garde les peaks entre -400 pb et +100 pb autour du TSS
TBP_rep1_prox <- subset(TBP_rep1_prox, distanceToTSS <= 100) 
head(TBP_rep1_prox)

nb_pb_peaks_TBP_rep1_distal <- sum(TBP_rep1_distal$width)
nb_pb_peaks_TBP_rep1_distal

nb_pb_peaks_TBP_rep1_prox <- sum(TBP_rep1_prox$width)
nb_pb_peaks_TBP_rep1_prox

nb_pb_peaks_TBP_rep1_total <- sum(TBP_rep1_annot$width)
nb_pb_peaks_TBP_rep1_total

nb_pb_peaks_TBP_rep1_non_distal <- nb_pb_peaks_TBP_rep1_total - nb_pb_peaks_TBP_rep1_distal
nb_pb_peaks_TBP_rep1_non_distal
nb_pb_peaks_TBP_rep1_non_prox <- nb_pb_peaks_TBP_rep1_total - nb_pb_peaks_TBP_rep1_prox
nb_pb_peaks_TBP_rep1_non_prox

# On exporte les fichiers contenant les peaks d'Oct4 distaux et proximaux
write.table(x = TBP_rep1_distal, file = "peaks_TBP_rep1_distaux.tsv", sep = "\t")
write.table(x = TBP_rep1_prox, file = "peaks_TBP_rep1_proximaux.tsv", sep = "\t")


#### TBP_rep1_Pol2_Chd8
TBP_rep1_Pol2_Chd8_avant2000 <- subset(TBP_rep1_Pol2_Chd8_intersection_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
TBP_rep1_Pol2_Chd8_apres2000 <- subset(TBP_rep1_Pol2_Chd8_intersection_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
TBP_rep1_Pol2_Chd8_distal <- rbind(TBP_rep1_Pol2_Chd8_avant2000, TBP_rep1_Pol2_Chd8_apres2000) 
head(TBP_rep1_Pol2_Chd8_distal)

TBP_rep1_Pol2_Chd8_prox <- subset(TBP_rep1_Pol2_Chd8_intersection_annot, distanceToTSS >= -400) # on garde les peaks entre -400 pb et +100 pb autour du TSS
TBP_rep1_Pol2_Chd8_prox <- subset(TBP_rep1_Pol2_Chd8_prox, distanceToTSS <= 100) 
head(TBP_rep1_Pol2_Chd8_prox)

nb_pb_peaks_TBP_rep1_Pol2_Chd8_distal <- sum(TBP_rep1_Pol2_Chd8_distal$width)
nb_pb_peaks_TBP_rep1_Pol2_Chd8_distal

nb_pb_peaks_TBP_rep1_Pol2_Chd8_prox <- sum(TBP_rep1_Pol2_Chd8_prox$width)
nb_pb_peaks_TBP_rep1_Pol2_Chd8_prox

nb_pb_peaks_TBP_rep1_Pol2_Chd8_total <- sum(TBP_rep1_Pol2_Chd8_intersection_annot$width)
nb_pb_peaks_TBP_rep1_Pol2_Chd8_total


nb_pb_peaks_TBP_rep1_Pol2_Chd8_non_distal <- nb_pb_peaks_TBP_rep1_Pol2_Chd8_total - nb_pb_peaks_TBP_rep1_Pol2_Chd8_distal
nb_pb_peaks_TBP_rep1_Pol2_Chd8_non_distal
nb_pb_peaks_TBP_rep1_Pol2_Chd8_non_prox <- nb_pb_peaks_TBP_rep1_Pol2_Chd8_total - nb_pb_peaks_TBP_rep1_Pol2_Chd8_prox
nb_pb_peaks_TBP_rep1_Pol2_Chd8_non_prox

# On exporte les fichiers contenant les peaks d'Oct4 distaux et proximaux
write.table(x = TBP_rep1_Pol2_Chd8_distal, file = "peaks_TBP_rep1_Pol2_Chd8_distaux.tsv", sep = "\t")
write.table(x = TBP_rep1_Pol2_Chd8_prox, file = "peaks_TBP_rep1_Pol2_Chd8_proximaux.tsv", sep = "\t")



#### CTCF
CTCF_avant2000 <- subset(CTCF_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
CTCF_apres2000 <- subset(CTCF_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
CTCF_distal <- rbind(CTCF_avant2000, CTCF_apres2000) 
head(CTCF_distal)

CTCF_prox <- subset(CTCF_annot, distanceToTSS >= -400) # on garde les peaks entre -400 pb et +100 pb autour du TSS
CTCF_prox <- subset(CTCF_prox, distanceToTSS <= 100) 
head(CTCF_prox)

nb_pb_peaks_CTCF_distal <- sum(CTCF_distal$width)
nb_pb_peaks_CTCF_distal

nb_pb_peaks_CTCF_prox <- sum(CTCF_prox$width)
nb_pb_peaks_CTCF_prox

nb_pb_peaks_CTCF_total <- sum(CTCF_annot$width)
nb_pb_peaks_CTCF_total


nb_pb_peaks_CTCF_non_distal <- nb_pb_peaks_CTCF_total - nb_pb_peaks_CTCF_distal
nb_pb_peaks_CTCF_non_distal
nb_pb_peaks_CTCF_non_prox <- nb_pb_peaks_CTCF_total - nb_pb_peaks_CTCF_prox
nb_pb_peaks_CTCF_non_prox

# On exporte les fichiers contenant les peaks de TBP distaux et proximaux
write.table(x = CTCF_distal, file = "peaks_CTCF_distaux.tsv", sep = "\t")
write.table(x = CTCF_prox, file = "peaks_CTCF_proximaux.tsv", sep = "\t")


#### Pol2
Pol2_avant2000 <- subset(Pol2_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
Pol2_apres2000 <- subset(Pol2_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
Pol2_distal <- rbind(Pol2_avant2000, Pol2_apres2000) 
head(Pol2_distal)

Pol2_prox <- subset(Pol2_annot, distanceToTSS >= -400) # on garde les peaks entre -400 pb et +100 pb autour du TSS
Pol2_prox <- subset(Pol2_prox, distanceToTSS <= 100) 
head(Pol2_prox)

nb_pb_peaks_Pol2_distal <- sum(Pol2_distal$width)
nb_pb_peaks_Pol2_distal

nb_pb_peaks_Pol2_prox <- sum(Pol2_prox$width)
nb_pb_peaks_Pol2_prox

nb_pb_peaks_Pol2_total <- sum(Pol2_annot$width)
nb_pb_peaks_Pol2_total


nb_pb_peaks_Pol2_non_distal <- nb_pb_peaks_Pol2_total - nb_pb_peaks_Pol2_distal
nb_pb_peaks_Pol2_non_distal
nb_pb_peaks_Pol2_non_prox <- nb_pb_peaks_Pol2_total - nb_pb_peaks_Pol2_prox
nb_pb_peaks_Pol2_non_prox

# On exporte les fichiers contenant les peaks de Pol2 distaux et proximaux
write.table(x = Pol2_distal, file = "peaks_Pol2_distaux.tsv", sep = "\t")
write.table(x = Pol2_prox, file = "peaks_Pol2_proximaux.tsv", sep = "\t")


#### Chd8
Chd8_avant2000 <- subset(Chd8_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
Chd8_apres2000 <- subset(Chd8_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
Chd8_distal <- rbind(Chd8_avant2000, Chd8_apres2000) 
head(Chd8_distal)

Chd8_prox <- subset(Chd8_annot, distanceToTSS >= -400) # on garde les peaks entre -400 pb et +100 pb autour du TSS
Chd8_prox <- subset(Chd8_prox, distanceToTSS <= 100) 
head(Chd8_prox)

nb_pb_peaks_Chd8_distal <- sum(Chd8_distal$width)
nb_pb_peaks_Chd8_distal

nb_pb_peaks_Chd8_prox <- sum(Chd8_prox$width)
nb_pb_peaks_Chd8_prox

nb_pb_peaks_Chd8_total <- sum(Chd8_annot$width)
nb_pb_peaks_Chd8_total


nb_pb_peaks_Chd8_non_distal <- nb_pb_peaks_Chd8_total - nb_pb_peaks_Chd8_distal
nb_pb_peaks_Chd8_non_distal
nb_pb_peaks_Chd8_non_prox <- nb_pb_peaks_Chd8_total - nb_pb_peaks_Chd8_prox
nb_pb_peaks_Chd8_non_prox

# On exporte les fichiers contenant les peaks de Chd8 distaux et proximaux
write.table(x = Chd8_distal, file = "peaks_Chd8_distaux.tsv", sep = "\t")
write.table(x = Chd8_prox, file = "peaks_Chd8_proximaux.tsv", sep = "\t")


#### Chd8_Oct4_merge
Chd8_Oct4_merge_avant2000 <- subset(Chd8_Oct4_intersection_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
Chd8_Oct4_merge_apres2000 <- subset(Chd8_Oct4_intersection_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
Chd8_Oct4_merge_distal <- rbind(Chd8_Oct4_merge_avant2000, Chd8_Oct4_merge_apres2000) 
head(Chd8_Oct4_merge_distal)

Chd8_Oct4_merge_prox <- subset(Chd8_Oct4_intersection_annot, distanceToTSS >= -400) # on garde les peaks entre -400 pb et +100 pb autour du TSS
Chd8_Oct4_merge_prox <- subset(Chd8_Oct4_merge_prox, distanceToTSS <= 100) 
head(Chd8_Oct4_merge_prox)

nb_pb_peaks_Chd8_Oct4_merge_distal <- sum(Chd8_Oct4_merge_distal$width)
nb_pb_peaks_Chd8_Oct4_merge_distal

nb_pb_peaks_Chd8_Oct4_merge_prox <- sum(Chd8_Oct4_merge_prox$width)
nb_pb_peaks_Chd8_Oct4_merge_prox

nb_pb_peaks_Chd8_Oct4_merge_total <- sum(Chd8_Oct4_intersection_annot$width)
nb_pb_peaks_Chd8_Oct4_merge_total

nb_pb_peaks_Chd8_Oct4_non_distal <- nb_pb_peaks_Chd8_Oct4_merge_total - nb_pb_peaks_Chd8_Oct4_merge_distal
nb_pb_peaks_Chd8_Oct4_non_distal
nb_pb_peaks_Chd8_Oct4_non_prox <- nb_pb_peaks_Chd8_Oct4_merge_total - nb_pb_peaks_Chd8_Oct4_merge_prox
nb_pb_peaks_Chd8_Oct4_non_prox

# On exporte les fichiers contenant les peaks de Chd8 distaux et proximaux
write.table(x = Chd8_Oct4_merge_distal, file = "peaks_Chd8_Oct4_merge_distaux.tsv", sep = "\t")
write.table(x = Chd8_Oct4_merge_prox, file = "peaks_Chd8_Oct4_merge_proximaux.tsv", sep = "\t")


#### Oct4_CTCF_merge
Oct4_CTCF_merge_avant2000 <- subset(Oct4_CTCF_intersection_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
Oct4_CTCF_merge_apres2000 <- subset(Oct4_CTCF_intersection_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
Oct4_CTCF_merge_distal <- rbind(Oct4_CTCF_merge_avant2000, Oct4_CTCF_merge_apres2000) 
head(Oct4_CTCF_merge_distal)

Oct4_CTCF_merge_prox <- subset(Oct4_CTCF_intersection_annot, distanceToTSS >= -400) # on garde les peaks entre -400 pb et +100 pb autour du TSS
Oct4_CTCF_merge_prox <- subset(Oct4_CTCF_merge_prox, distanceToTSS <= 100) 
head(Oct4_CTCF_merge_prox)

nb_pb_peaks_Oct_CTCF_merge_distal <- sum(Oct4_CTCF_merge_distal$width)
nb_pb_peaks_Oct_CTCF_merge_distal

nb_pb_peaks_Oct_CTCF_merge_prox <- sum(Oct4_CTCF_merge_prox$width)
nb_pb_peaks_Oct_CTCF_merge_prox

nb_pb_peaks_Oct_CTCF_merge_total <- sum(Oct4_CTCF_intersection_annot$width)
nb_pb_peaks_Oct_CTCF_merge_total


nb_pb_peaks_Oct_CTCF_non_distal <- nb_pb_peaks_Oct_CTCF_merge_total - nb_pb_peaks_Oct_CTCF_merge_distal
nb_pb_peaks_Oct_CTCF_non_distal
nb_pb_peaks_Oct_CTCF_non_prox <- nb_pb_peaks_Oct_CTCF_merge_total - nb_pb_peaks_Oct_CTCF_merge_prox
nb_pb_peaks_Oct_CTCF_non_prox

# On exporte les fichiers contenant les peaks d'Oct4_CTCF_merge distaux et proximaux
write.table(x = Oct4_CTCF_merge_distal, file = "peaks_Oct4_CTCF_merge_distaux.tsv", sep = "\t")
write.table(x = Oct4_CTCF_merge_prox, file = "peaks_Oct4_CTCF_merge_proximaux.tsv", sep = "\t")



#### Génome de mm10
samplefiles_genome <- list.files("/home/mnocente/Bureau/", pattern= ".bed",full.names = T)
samplefiles_genome <- as.list(samplefiles_genome)
names(samplefiles_genome) <- c("genome_mm10")

peakAnnoList_genome <- lapply(samplefiles_genome, annotatePeak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)
print(peakAnnoList_genome)

genome_annot <- as.data.frame(peakAnnoList_genome[["genome_mm10"]]@anno)
sum(genome_annot$width)


############################################################################################################
#### Tables de contingence des differents facteurs et intersections ####
#### Test du Chi2 pour tester la significativité de l'enrichissement du facteur sur la region d'interet ####
############################################################################################################
nb_pb_genome_mm10 <- 2.72554e+09

nb_pb_region_distale_mm10 <- 2420292007
nb_pb_region_proximale_mm10_400 <- 38288993

nb_pb_genome_sans_distale <- 2.72554e+09 - 2420292007
nb_pb_genome_sans_proximale <- 2.72554e+09 - 38288993


#### Oct4
### Oct4 au niveau distal
nombre_pb_aux_peaks_Oct4_dist <- c(nb_pb_peaks_Oct_distal, nb_pb_peaks_Oct_non_distal)
nombre_pb_aux_peaks_Oct4_dist
nombre_pb_sur_genome_entier_dist <- c(nb_pb_region_distale_mm10, nb_pb_genome_sans_distale)
nombre_pb_sur_genome_entier_dist


contingence_table_Oct4_dist <- data.frame(nombre_pb_aux_peaks_Oct4_dist, nombre_pb_sur_genome_entier_dist, stringsAsFactors = FALSE)
contingence_table_Oct4_dist

colnames(contingence_table_Oct4_dist) <- c("nombre_pb_aux_peaks_Oct4", "nombre_pb_sur_genome_entier")
contingence_table_Oct4_dist

rownames(contingence_table_Oct4_dist) <- c("region_distale", "genome_complet_sans_region_distale")
contingence_table_Oct4_dist


test_Chi2_Oct4_distal <- chisq.test(x = contingence_table_Oct4_dist, y = NULL, correct = TRUE,
                                   p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                   simulate.p.value = FALSE, B = 2000)
test_Chi2_Oct4_distal


### Oct4 au niveau proximal : -400 : +100
nombre_pb_aux_peaks_Oct4_prox <- c(nb_pb_peaks_Oct_prox, nb_pb_peaks_Oct_non_prox)
nombre_pb_aux_peaks_Oct4_prox
nombre_pb_sur_genome_entier_prox <- c(nb_pb_region_proximale_mm10_400, nb_pb_genome_sans_proximale)
nombre_pb_sur_genome_entier_prox

contingence_table_Oct4_prox <- data.frame(nombre_pb_aux_peaks_Oct4_prox, nombre_pb_sur_genome_entier_prox, stringsAsFactors = FALSE)
contingence_table_Oct4_prox

colnames(contingence_table_Oct4_prox) <- c("nombre_pb_aux_peaks_Oct4", "nombre_pb_sur_genome_entier")
contingence_table_Oct4_prox

rownames(contingence_table_Oct4_prox) <- c("region_proximale", "genome_complet_sans_region_proximale")
contingence_table_Oct4_prox


testChi2_Oct4_prox <- chisq.test(x = contingence_table_Oct4_prox, y = NULL, correct = TRUE,
                                 p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                 simulate.p.value = FALSE, B = 2000)

testChi2_Oct4_prox


#### TBP_rep1
### TBP_rep1 au niveau distal
nombre_pb_aux_peaks_TBP_rep1_dist <- c(nb_pb_peaks_TBP_rep1_distal, nb_pb_peaks_TBP_rep1_non_distal)
nombre_pb_aux_peaks_TBP_rep1_dist
nombre_pb_sur_genome_entier_dist <- c(nb_pb_region_distale_mm10, nb_pb_genome_sans_distale)
nombre_pb_sur_genome_entier_dist

contingence_table_TBP_rep1_dist <- data.frame(nombre_pb_aux_peaks_TBP_rep1_dist, nombre_pb_sur_genome_entier_dist, stringsAsFactors = FALSE)
contingence_table_TBP_rep1_dist

colnames(contingence_table_TBP_rep1_dist) <- c("nombre_pb_aux_peaks_TBP_rep1", "nombre_pb_sur_genome_entier")
contingence_table_TBP_rep1_dist

rownames(contingence_table_TBP_rep1_dist) <- c("region_distale", "genome_complet_sans_region_distale")
contingence_table_TBP_rep1_dist


test_Chi2_TBP_rep1_distal <- chisq.test(x = contingence_table_TBP_rep1_dist, y = NULL, correct = TRUE,
                                    p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                    simulate.p.value = FALSE, B = 2000)
test_Chi2_TBP_rep1_distal


### TBP_rep1 au niveau proximal : -400 : +100
nombre_pb_aux_peaks_TBP_rep1_prox <- c(nb_pb_peaks_TBP_rep1_prox, nb_pb_peaks_TBP_rep1_non_prox)
nombre_pb_aux_peaks_TBP_rep1_prox
nombre_pb_sur_genome_entier_prox <- c(nb_pb_region_proximale_mm10_400, nb_pb_genome_sans_proximale)
nombre_pb_sur_genome_entier_prox

contingence_table_TBP_rep1_prox <- data.frame(nombre_pb_aux_peaks_TBP_rep1_prox, nombre_pb_sur_genome_entier_prox, stringsAsFactors = FALSE)
contingence_table_TBP_rep1_prox

colnames(contingence_table_TBP_rep1_prox) <- c("nombre_pb_aux_peaks_TBP_rep1", "nombre_pb_sur_genome_entier")
contingence_table_TBP_rep1_prox

rownames(contingence_table_TBP_rep1_prox) <- c("region_proximale", "genome_complet_sans_region_proximale")
contingence_table_TBP_rep1_prox


testChi2_TBP_rep1_prox <- chisq.test(x = contingence_table_TBP_rep1_prox, y = NULL, correct = TRUE,
                                 p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                 simulate.p.value = FALSE, B = 2000)

testChi2_TBP_rep1_prox


#### Pol2
### Pol2 au niveau distal
nombre_pb_aux_peaks_Pol2_dist <- c(nb_pb_peaks_Pol2_distal, nb_pb_peaks_Pol2_non_distal)
nombre_pb_aux_peaks_Pol2_dist
nombre_pb_sur_genome_entier_dist <- c(nb_pb_region_distale_mm10, nb_pb_genome_sans_distale)
nombre_pb_sur_genome_entier_dist

contingence_table_Pol2_dist <- data.frame(nombre_pb_aux_peaks_Pol2_dist, nombre_pb_sur_genome_entier_dist, stringsAsFactors = FALSE)
contingence_table_Pol2_dist

colnames(contingence_table_Pol2_dist) <- c("nombre_pb_aux_peaks_Pol2", "nombre_pb_sur_genome_entier")
contingence_table_Pol2_dist

rownames(contingence_table_Pol2_dist) <- c("region_distale", "genome_complet_sans_region_distale")
contingence_table_Pol2_dist


test_Chi2_Pol2_distal <- chisq.test(x = contingence_table_Pol2_dist, y = NULL, correct = TRUE,
                                   p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                   simulate.p.value = FALSE, B = 2000)
test_Chi2_Pol2_distal


### Pol2 au niveau proximal : -400 : +100
nombre_pb_aux_peaks_Pol2_prox <- c(nb_pb_peaks_Pol2_prox, nb_pb_peaks_Pol2_non_prox)
nombre_pb_aux_peaks_Pol2_prox
nombre_pb_sur_genome_entier_prox <- c(nb_pb_region_proximale_mm10_400, nb_pb_genome_sans_proximale)
nombre_pb_sur_genome_entier_prox

contingence_table_Pol2_prox <- data.frame(nombre_pb_aux_peaks_Pol2_prox, nombre_pb_sur_genome_entier_prox, stringsAsFactors = FALSE)
contingence_table_Pol2_prox

colnames(contingence_table_Pol2_prox) <- c("nombre_pb_aux_peaks_Pol2", "nombre_pb_sur_genome_entier")
contingence_table_Pol2_prox

rownames(contingence_table_Pol2_prox) <- c("region_proximale", "genome_complet_sans_region_proximale")
contingence_table_Pol2_prox


testChi2_Pol2_prox <- chisq.test(x = contingence_table_Pol2_prox, y = NULL, correct = TRUE,
                                p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                simulate.p.value = FALSE, B = 2000)

testChi2_Pol2_prox


#### CTCF
### CTCF au niveau distal
nombre_pb_aux_peaks_CTCF_dist <- c(nb_pb_peaks_CTCF_distal, nb_pb_peaks_CTCF_non_distal)
nombre_pb_aux_peaks_CTCF_dist
nombre_pb_sur_genome_entier_dist <- c(nb_pb_region_distale_mm10, nb_pb_genome_sans_distale)
nombre_pb_sur_genome_entier_dist

contingence_table_CTCF_dist <- data.frame(nombre_pb_aux_peaks_CTCF_dist, nombre_pb_sur_genome_entier_dist, stringsAsFactors = FALSE)
contingence_table_CTCF_dist

colnames(contingence_table_CTCF_dist) <- c("nombre_pb_aux_peaks_CTCF", "nombre_pb_sur_genome_entier")
contingence_table_CTCF_dist

rownames(contingence_table_CTCF_dist) <- c("region_distale", "genome_complet_sans_region_distale")
contingence_table_CTCF_dist


test_Chi2_CTCF_distal <- chisq.test(x = contingence_table_CTCF_dist, y = NULL, correct = TRUE,
                                    p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                    simulate.p.value = FALSE, B = 2000)
test_Chi2_CTCF_distal


### CTCF au niveau proximal : -400 : +100
nombre_pb_aux_peaks_CTCF_prox <- c(nb_pb_peaks_CTCF_prox, nb_pb_peaks_CTCF_non_prox)
nombre_pb_aux_peaks_CTCF_prox
nombre_pb_sur_genome_entier_prox <- c(nb_pb_region_proximale_mm10_400, nb_pb_genome_sans_proximale)
nombre_pb_sur_genome_entier_prox

contingence_table_CTCF_prox <- data.frame(nombre_pb_aux_peaks_CTCF_prox, nombre_pb_sur_genome_entier_prox, stringsAsFactors = FALSE)
contingence_table_CTCF_prox

colnames(contingence_table_CTCF_prox) <- c("nombre_pb_aux_peaks_CTCF", "nombre_pb_sur_genome_entier")
contingence_table_CTCF_prox

rownames(contingence_table_CTCF_prox) <- c("region_proximale", "genome_complet_sans_region_proximale")
contingence_table_CTCF_prox


testChi2_CTCF_prox <- chisq.test(x = contingence_table_CTCF_prox, y = NULL, correct = TRUE,
                                 p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                 simulate.p.value = FALSE, B = 2000)

testChi2_CTCF_prox


#### Chd8
### Chd8 au niveau distal
nombre_pb_aux_peaks_Chd8_dist <- c(nb_pb_peaks_Chd8_distal, nb_pb_peaks_Chd8_non_distal)
nombre_pb_aux_peaks_Chd8_dist
nombre_pb_sur_genome_entier_dist <- c(nb_pb_region_distale_mm10, nb_pb_genome_sans_distale)
nombre_pb_sur_genome_entier_dist

contingence_table_Chd8_dist <- data.frame(nombre_pb_aux_peaks_Chd8_dist, nombre_pb_sur_genome_entier_dist, stringsAsFactors = FALSE)
contingence_table_Chd8_dist

colnames(contingence_table_Chd8_dist) <- c("nombre_pb_aux_peaks_Chd8", "nombre_pb_sur_genome_entier")
contingence_table_Chd8_dist

rownames(contingence_table_Chd8_dist) <- c("region_distale", "genome_complet_sans_region_distale")
contingence_table_Chd8_dist


test_Chi2_Chd8_distal <- chisq.test(x = contingence_table_Chd8_dist, y = NULL, correct = TRUE,
                                    p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                    simulate.p.value = FALSE, B = 2000)
test_Chi2_Chd8_distal


### Chd8 au niveau proximal : -400 : +100
nombre_pb_aux_peaks_Chd8_prox <- c(nb_pb_peaks_Chd8_prox, nb_pb_peaks_Chd8_non_prox)
nombre_pb_aux_peaks_Chd8_prox
nombre_pb_sur_genome_entier_prox <- c(nb_pb_region_proximale_mm10_400, nb_pb_genome_sans_proximale)
nombre_pb_sur_genome_entier_prox

contingence_table_Chd8_prox <- data.frame(nombre_pb_aux_peaks_Chd8_prox, nombre_pb_sur_genome_entier_prox, stringsAsFactors = FALSE)
contingence_table_Chd8_prox

colnames(contingence_table_Chd8_prox) <- c("nombre_pb_aux_peaks_Chd8", "nombre_pb_sur_genome_entier")
contingence_table_Chd8_prox

rownames(contingence_table_Chd8_prox) <- c("region_proximale", "genome_complet_sans_region_proximale")
contingence_table_Chd8_prox


testChi2_Chd8_prox <- chisq.test(x = contingence_table_Chd8_prox, y = NULL, correct = TRUE,
                                 p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                 simulate.p.value = FALSE, B = 2000)

testChi2_Chd8_prox



#### Chd8_Oct4_merge
### Chd8_Oct4_merge au niveau distal
nombre_pb_aux_peaks_Chd8_Oct4_dist <- c(nb_pb_peaks_Chd8_Oct4_merge_distal, nb_pb_peaks_Chd8_Oct4_non_distal)
nombre_pb_aux_peaks_Chd8_Oct4_dist
nombre_pb_sur_genome_entier_dist <- c(nb_pb_region_distale_mm10, nb_pb_genome_sans_distale)
nombre_pb_sur_genome_entier_dist

contingence_table_Chd8_Oct4_dist <- data.frame(nombre_pb_aux_peaks_Chd8_Oct4_dist, nombre_pb_sur_genome_entier_dist, stringsAsFactors = FALSE)
contingence_table_Chd8_Oct4_dist

colnames(contingence_table_Chd8_Oct4_dist) <- c("nombre_pb_aux_peaks_Chd8_Oct4_inter", "nombre_pb_sur_genome_entier")
contingence_table_Chd8_Oct4_dist

rownames(contingence_table_Chd8_Oct4_dist) <- c("region_distale", "genome_complet_sans_region_distale")
contingence_table_Chd8_Oct4_dist


test_Chi2_Chd8_Oct4_inter_distal <- chisq.test(x = contingence_table_Chd8_Oct4_dist, y = NULL, correct = TRUE,
                                    p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                    simulate.p.value = FALSE, B = 2000)
test_Chi2_Chd8_Oct4_inter_distal




### Chd8_Oct4_merge au niveau proximal : -400 : +100
nombre_pb_aux_peaks_Chd8_Oct4_merge_prox <- c(nb_pb_peaks_Chd8_Oct4_merge_prox, nb_pb_peaks_Chd8_Oct4_non_prox)
nombre_pb_aux_peaks_Chd8_Oct4_merge_prox
nombre_pb_sur_genome_entier_prox <- c(nb_pb_region_proximale_mm10_400, nb_pb_genome_sans_proximale)
nombre_pb_sur_genome_entier_prox

contingence_table_Chd8_Oct4_inter_prox <- data.frame(nombre_pb_aux_peaks_Chd8_Oct4_merge_prox, nombre_pb_sur_genome_entier_prox, stringsAsFactors = FALSE)
contingence_table_Chd8_Oct4_inter_prox

colnames(contingence_table_Chd8_Oct4_inter_prox) <- c("nombre_pb_aux_peaks_Chd8_Oct4_merge", "nombre_pb_sur_genome_entier")
contingence_table_Chd8_Oct4_inter_prox

rownames(contingence_table_Chd8_Oct4_inter_prox) <- c("region_proximale", "genome_complet_sans_region_proximale")
contingence_table_Chd8_Oct4_inter_prox


testChi2_Chd8_Oct4_inter_prox <- chisq.test(x = contingence_table_Chd8_Oct4_inter_prox, y = NULL, correct = TRUE,
                                 p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                 simulate.p.value = FALSE, B = 2000)

testChi2_Chd8_Oct4_inter_prox


#### TBP_rep1_Pol2_Chd8
### TBP_rep1_Pol2_Chd8 au niveau distal
nombre_pb_aux_peaks_TBP_rep1_Pol2_Chd8_dist <- c(nb_pb_peaks_TBP_rep1_Pol2_Chd8_distal, nb_pb_peaks_TBP_rep1_Pol2_Chd8_non_distal)
nombre_pb_aux_peaks_TBP_rep1_Pol2_Chd8_dist
nombre_pb_sur_genome_entier_dist <- c(nb_pb_region_distale_mm10, nb_pb_genome_sans_distale)
nombre_pb_sur_genome_entier_dist

contingence_table_TBP_rep1_Pol2_Chd8_inter_dist <- data.frame(nombre_pb_aux_peaks_TBP_rep1_Pol2_Chd8_dist, nombre_pb_sur_genome_entier_dist, stringsAsFactors = FALSE)
contingence_table_TBP_rep1_Pol2_Chd8_inter_dist

colnames(contingence_table_TBP_rep1_Pol2_Chd8_inter_dist) <- c("nombre_pb_aux_peaks_TBP_rep1_Pol2_Chd8_inter", "nombre_pb_sur_genome_entier")
contingence_table_TBP_rep1_Pol2_Chd8_inter_dist

rownames(contingence_table_TBP_rep1_Pol2_Chd8_inter_dist) <- c("region_distale", "genome_complet_sans_region_distale")
contingence_table_TBP_rep1_Pol2_Chd8_inter_dist


test_Chi2_TBP_rep1_Pol2_Chd8_inter_distal <- chisq.test(x = contingence_table_TBP_rep1_Pol2_Chd8_inter_dist, y = NULL, correct = TRUE,
                                    p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                    simulate.p.value = FALSE, B = 2000)
test_Chi2_TBP_rep1_Pol2_Chd8_inter_distal


### TBP_rep1_Pol2_Chd8 au niveau proximal : -400 : +100
nombre_pb_aux_peaks_TBP_rep1_Pol2_Chd8_prox <- c(nb_pb_peaks_TBP_rep1_Pol2_Chd8_prox, nb_pb_peaks_TBP_rep1_Pol2_Chd8_non_prox)
nombre_pb_aux_peaks_TBP_rep1_Pol2_Chd8_prox
nombre_pb_sur_genome_entier_prox <- c(nb_pb_region_proximale_mm10_400, nb_pb_genome_sans_proximale)
nombre_pb_sur_genome_entier_prox

contingence_table_TBP_rep1_Pol2_Chd8_inter_prox <- data.frame(nombre_pb_aux_peaks_TBP_rep1_Pol2_Chd8_prox, nombre_pb_sur_genome_entier_prox, stringsAsFactors = FALSE)
contingence_table_TBP_rep1_Pol2_Chd8_inter_prox

colnames(contingence_table_TBP_rep1_Pol2_Chd8_inter_prox) <- c("nombre_pb_aux_peaks_TBP_rep1_Pol2_Chd8_inter", "nombre_pb_sur_genome_entier")
contingence_table_TBP_rep1_Pol2_Chd8_inter_prox

rownames(contingence_table_TBP_rep1_Pol2_Chd8_inter_prox) <- c("region_proximale", "genome_complet_sans_region_proximale")
contingence_table_TBP_rep1_Pol2_Chd8_inter_prox


testChi2_TBP_rep1_Pol2_Chd8_inter_prox <- chisq.test(x = contingence_table_TBP_rep1_Pol2_Chd8_inter_prox, y = NULL, correct = TRUE,
                                 p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                 simulate.p.value = FALSE, B = 2000)

testChi2_TBP_rep1_Pol2_Chd8_inter_prox



#### Oct4_CTCF_inter
### Oct4_CTCF_inter au niveau distal
nombre_pb_aux_peaks_Oct4_CTCF_merge_dist <- c(nb_pb_peaks_Oct_CTCF_merge_distal, nb_pb_peaks_Oct_CTCF_non_distal)
nombre_pb_aux_peaks_Oct4_CTCF_merge_dist
nombre_pb_sur_genome_entier_dist <- c(nb_pb_region_distale_mm10, nb_pb_genome_sans_distale)
nombre_pb_sur_genome_entier_dist

contingence_table_Oct4_CTCF_inter_dist <- data.frame(nombre_pb_aux_peaks_Oct4_CTCF_merge_dist, nombre_pb_sur_genome_entier_dist, stringsAsFactors = FALSE)
contingence_table_Oct4_CTCF_inter_dist

colnames(contingence_table_Oct4_CTCF_inter_dist) <- c("nombre_pb_aux_peaks_Oct4_CTCF_inter", "nombre_pb_sur_genome_entier")
contingence_table_Oct4_CTCF_inter_dist

rownames(contingence_table_Oct4_CTCF_inter_dist) <- c("region_distale", "genome_complet_sans_region_distale")
contingence_table_Oct4_CTCF_inter_dist


test_Chi2_Oct4_CTCF_inter_distal <- chisq.test(x = contingence_table_Oct4_CTCF_inter_dist, y = NULL, correct = TRUE,
                                    p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                    simulate.p.value = FALSE, B = 2000)
test_Chi2_Oct4_CTCF_inter_distal


### Oct4_CTCF_inter au niveau proximal : -400 : +100
nombre_pb_aux_peaks_Oct4_CTCF_merge_prox <- c(nb_pb_peaks_Oct_CTCF_merge_prox, nb_pb_peaks_Oct_CTCF_non_prox)
nombre_pb_aux_peaks_Oct4_CTCF_merge_prox
nombre_pb_sur_genome_entier_prox <- c(nb_pb_region_proximale_mm10_400, nb_pb_genome_sans_proximale)
nombre_pb_sur_genome_entier_prox

contingence_table_Oct4_CTCF_inter_prox <- data.frame(nombre_pb_aux_peaks_Oct4_CTCF_merge_prox, nombre_pb_sur_genome_entier_prox, stringsAsFactors = FALSE)
contingence_table_Oct4_CTCF_inter_prox

colnames(contingence_table_Oct4_CTCF_inter_prox) <- c("nombre_pb_aux_peaks_Oct4_CTCF_inter", "nombre_pb_sur_genome_entier")
contingence_table_Oct4_CTCF_inter_prox

rownames(contingence_table_Oct4_CTCF_inter_prox) <- c("region_proximale", "genome_complet_sans_region_proximale")
contingence_table_Oct4_CTCF_inter_prox


testChi2_Oct4_CTCF_inter_prox <- chisq.test(x = contingence_table_Oct4_CTCF_inter_prox, y = NULL, correct = TRUE,
                                 p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                 simulate.p.value = FALSE, B = 2000)

testChi2_Oct4_CTCF_inter_prox


