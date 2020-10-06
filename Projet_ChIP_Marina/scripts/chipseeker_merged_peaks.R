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
# les replicats de Pol2, TBP, CHd8 et Oct4 ont été mergés. On n'avait malheureusement pas de duplicat pour CTCF donc j'ai gardé tous les peaks.
samplefiles <- list.files("/home/mnocente/Bureau/Stage_projet_Marina/Projet_ChIP_Marina/results/macs2_peaks/merge/final", pattern= ".narrowPeak",full.names = T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("CTCF", "Chd8", "Oct4", "Pol2", "TBP")

print("les fichiers ont bien été chargés")
print(samplefiles)

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


##### Visualisation des annotations #####

## Barchart (multiple samples for comparison)
png(file = "Barchart_comparison_sample_annotation.png")
plotAnnoBar(peakAnnoList)
dev.off()

## Distribution of TF-binding loci relative to TSS
png(file = "Distribution_of_TF-binding_loci_relative_to_TSS.png")
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")
dev.off()


### Get annotation data frame
library(dplyr)
Oct4_annot <- as.data.frame(peakAnnoList[["Oct4"]]@anno)
TBP_annot <- as.data.frame(peakAnnoList[["TBP"]]@anno)
CTCF_annot <- as.data.frame(peakAnnoList[["CTCF"]]@anno)
Pol2_annot <- as.data.frame(peakAnnoList[["Pol2"]]@anno)
Chd8_annot <- as.data.frame(peakAnnoList[["Chd8"]]@anno)

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
entrezids_TBP <- TBP_annot$geneId %>% as.character() %>% unique()

ego_TBP <- enrichGO(gene = entrezids_TBP, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

## Output results from GO analysis to a table
#cluster_summary_TBP <- data.frame(ego_TBP)
#write.csv(cluster_summary_TBP, "results/clusterProfiler_TBP.csv")

## Dotplot visualization
png(file = "dotplotfunctional_enrichment_TBP.png", width = 750, height = 1000)
dotplot(ego_TBP, showCategory=20)
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

png(file = "KEGG_Pathway_Enrichment_Analysis.png", width = 1100, height = 1100)
dotplot(compKEGG, showCategory = 20, title = "KEGG_Pathway_Enrichment_Analysis")
dev.off()




#### Overlap of peaks and annotated genes ####
png(file = "Vennplot_allTF.png", width = 500, height = 500)
peakAnnoList_all <- c(peakAnnoList$Oct4, peakAnnoList$CTCF, peakAnnoList$Pol2, peakAnnoList$TBP, peakAnnoList$Chd8)
names(peakAnnoList_all) <- c("Oct4", "CTCF", "Pol2", "TBP", "Chd8")
peakAnnoList_all
genesOverlap= lapply(peakAnnoList_all, function(i) as.data.frame(i)$geneId)
vennplot(genesOverlap)
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



#### TBP
TBP_avant2000 <- subset(TBP_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
TBP_apres2000 <- subset(TBP_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
TBP_distal <- rbind(TBP_avant2000, TBP_apres2000) 
head(TBP_distal)

TBP_prox <- subset(TBP_annot, distanceToTSS >= -400) # on garde les peaks entre -400 pb et +100 pb autour du TSS
TBP_prox <- subset(TBP_prox, distanceToTSS <= 100) 
head(TBP_prox)

nb_pb_peaks_TBP_distal <- sum(TBP_distal$width)
nb_pb_peaks_TBP_distal

nb_pb_peaks_TBP_prox <- sum(TBP_prox$width)
nb_pb_peaks_TBP_prox

nb_pb_peaks_TBP_total <- sum(TBP_annot$width)
nb_pb_peaks_TBP_total


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


#### Génome de mm10
samplefiles_genome <- list.files("/home/mnocente/Bureau/", pattern= ".bed",full.names = T)
samplefiles_genome <- as.list(samplefiles_genome)
names(samplefiles_genome) <- c("genome_mm10")

peakAnnoList_genome <- lapply(samplefiles_genome, annotatePeak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)
print(peakAnnoList_genome)

genome_annot <- as.data.frame(peakAnnoList_genome[["genome_mm10"]]@anno)
sum(genome_annot$width)





####################################################################
#### Peaks communs Oct4_CTCF ####
####################################################################
##### Loading des data #####
# As input we need to provide the names of our BED files in a list format.
# les replicats de Oct4 ont été mergés puis avec CTCF. On regarde les peaks communs à CTCF et Oct4.
samplefiles_merge <- list.files("/home/mnocente/Bureau/Stage_projet_Marina/Projet_ChIP_Marina/results/macs2_peaks/merge/final/CTCF_Oct4", pattern= ".narrowPeak",full.names = T)
samplefiles_merge <- as.list(samplefiles_merge)
names(samplefiles_merge) <- c("CTCF_Oct4")

print(samplefiles_merge)

##### Assign annotation db #####

# We need to assign annotation databases generated from UCSC to a variable:

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
print("les annotations ont bien été chargés")

# ChIPseeker et sa fonction annotatePeak permettent d'annoter les pics avec le gène et la région génomique les plus proches de là où se trouve le pic.
# The annotatePeak function by default uses the TSS method, and provides parameters to specify a max distance cutoff.


##### Get annotations #####
peakAnnoList_merge <- lapply(samplefiles_merge, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)

print(peakAnnoList_merge) # annotation information is stored in the peakAnnoList_merge


##### Visualisation des annotations #####

## Barchart (multiple samples for comparison)
png(file = "Barchart_comparison_sample_annotation_CTCF_Oct4_commun.png")
plotAnnoBar(peakAnnoList_merge)
dev.off()

## Distribution of TF-binding loci relative to TSS
png(file = "Distribution_of_TF-binding_loci_relative_to_TSS_CTCF_Oct4_commun.png")
plotDistToTSS(peakAnnoList_merge, title="Distribution of transcription factor-binding loci \n relative to TSS")
dev.off()


### Get annotation data frame
library(dplyr)
CTCF_Oct4_annot <- as.data.frame(peakAnnoList_merge[["CTCF_Oct4"]]@anno)
head(CTCF_Oct4_annot$distanceToTSS)


### La liste de gènes des annotations Oct4 et CTCF commun en entree pour une analyse d'enrichissement GO.

## Run GO enrichment analysis 
entrezids_CTCF_Oct4 <- CTCF_Oct4_annot$geneId %>% as.character() %>% unique()

ego_CTCF_Oct4 <- enrichGO(gene = entrezids_CTCF_Oct4, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)


## Dotplot visualization
png(file = "dotplotfunctional_enrichment_CTCF_Oct4_commun.png", width = 750, height = 1000)
dotplot(ego_CTCF_Oct4, showCategory=50)
dev.off()


####################################################################
#### Tables de contingence et test de Fisher ####
####################################################################

### Creation de la table de contingence pour les éléments distaux du TSS:
## Oct4 rep1
Oct4_rep1_avant2000 <- subset(Oct4_rep1_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
Oct4_rep1_apres2000 <- subset(Oct4_rep1_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
Oct4_rep1_distal <- rbind(Oct4_rep1_avant2000, Oct4_rep1_apres2000) 
head(Oct4_rep1_distal)
min(Oct4_rep1_distal$distanceToTSS)
max(Oct4_rep1_distal$distanceToTSS)
dim(Oct4_rep1_distal) # il y a 30814 peaks distaux du TSS (avant -2000 et apres +2000)
dim(Oct4_rep1_annot) # il y a en tout 42848 peaks Oct4

nb_peaks_Oct4_rep1_distal <- nrow(Oct4_rep1_distal) 
nb_peaks_Oct4_rep1_distal

nb_peaks_Oct4_rep1_NON_distal <- nrow(Oct4_rep1_annot) - nrow(Oct4_rep1_distal)
nb_peaks_Oct4_rep1_NON_distal


count_distance_data <- function(df, dist) {
  # Recover distance information
  min_threshold <- df$distanceToTSS < -dist;
  max_threshold <- df$distanceToTSS >  dist;
  distal <- min_threshold | max_threshold;
  # Compute results
  nb_distal <- table(distal)["TRUE"];
  nb_non_distal <- table(distal)["FALSE"];
  total <- length(distal);
  # Build results object
  results <- c(nb_distal, nb_non_distal, total);
  names(results) <- c("nb_distal", "nb_non_distal", "total");
  return(results)
}

Nb_peaks_Oct4_rep1_fct <- count_distance_data(Oct4_rep1_annot, 2000)
print(Nb_peaks_Oct4_rep1_fct)





## Oct4 rep2
Oct4_rep2_avant2000 <- subset(Oct4_rep2_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
Oct4_rep2_apres2000 <- subset(Oct4_rep2_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
Oct4_rep2_distal <- rbind(Oct4_rep2_avant2000, Oct4_rep2_apres2000) 
head(Oct4_rep2_distal)
min(Oct4_rep2_distal$distanceToTSS)
max(Oct4_rep2_distal$distanceToTSS)
dim(Oct4_rep2_distal) # il y a 24358 peaks distaux du TSS (avant -2000 et apres +2000)
dim(Oct4_rep2_annot) # il y a en tout 34111 peaks Oct4

nb_peaks_Oct4_rep2_distal <- nrow(Oct4_rep2_distal) 
nb_peaks_Oct4_rep2_distal

nb_peaks_Oct4_rep2_NON_distal <- nrow(Oct4_rep2_annot) - nrow(Oct4_rep2_distal)
nb_peaks_Oct4_rep2_NON_distal

# On ajoute les réplicats:
nb_peaks_Oct4_distal <- nb_peaks_Oct4_rep1_distal + nb_peaks_Oct4_rep2_distal
nb_peaks_Oct4_distal
nb_peaks_Oct4_NON_distal <- nb_peaks_Oct4_rep1_NON_distal + nb_peaks_Oct4_rep2_NON_distal
nb_peaks_Oct4_NON_distal

## CTCF
CTCF_avant2000 <- subset(CTCF_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
CTCF_apres2000 <- subset(CTCF_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
CTCF_distal <- rbind(CTCF_avant2000, CTCF_apres2000) 
head(CTCF_distal)
min(CTCF_distal$distanceToTSS)
max(CTCF_distal$distanceToTSS)
dim(CTCF_distal) # il y a 43091 peaks distaux du TSS (avant -2000 et apres +2000)
dim(CTCF_annot) # il y a en tout 52760 peaks CTCF

nb_peaks_CTCF_distal <- nrow(CTCF_distal) 
nb_peaks_CTCF_distal

nb_peaks_CTCF_NON_distal <- nrow(CTCF_annot) - nrow(CTCF_distal)
nb_peaks_CTCF_NON_distal


## TBP rep1
TBP_rep1_avant2000 <- subset(TBP_rep1_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
TBP_rep1_apres2000 <- subset(TBP_rep1_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
TBP_rep1_distal <- rbind(TBP_rep1_avant2000, TBP_rep1_apres2000) 
head(TBP_rep1_distal)
min(TBP_rep1_distal$distanceToTSS)
max(TBP_rep1_distal$distanceToTSS)
dim(TBP_rep1_distal) # il y a 32520 peaks distaux du TSS (avant -2000 et apres +2000)
dim(TBP_rep1_annot) # il y a en tout 49448 peaks TBP

nb_peaks_TBP_rep1_distal <- nrow(TBP_rep1_distal) 
nb_peaks_TBP_rep1_distal

nb_peaks_TBP_rep1_NON_distal <- nrow(TBP_rep1_annot) - nrow(TBP_rep1_distal)
nb_peaks_TBP_rep1_NON_distal


## TBP rep2
TBP_rep2_avant2000 <- subset(TBP_rep2_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
TBP_rep2_apres2000 <- subset(TBP_rep2_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
TBP_rep2_distal <- rbind(TBP_rep2_avant2000, TBP_rep2_apres2000) 
head(TBP_rep2_distal)
min(TBP_rep2_distal$distanceToTSS)
max(TBP_rep2_distal$distanceToTSS)
dim(TBP_rep2_distal) # il y a 435 peaks distaux du TSS (avant -2000 et apres +2000)
dim(TBP_rep2_annot) # il y a en tout 743 peaks TBP

nb_peaks_TBP_rep2_distal <- nrow(TBP_rep2_distal) 
nb_peaks_TBP_rep2_distal

nb_peaks_TBP_rep2_NON_distal <- nrow(TBP_rep2_annot) - nrow(TBP_rep2_distal)
nb_peaks_TBP_rep2_NON_distal


# On ajoute les réplicats:
nb_peaks_TBP_distal <- nb_peaks_TBP_rep1_distal + nb_peaks_TBP_rep2_distal
nb_peaks_TBP_distal
nb_peaks_TBP_NON_distal <- nb_peaks_TBP_rep1_NON_distal + nb_peaks_TBP_rep2_NON_distal
nb_peaks_TBP_NON_distal


## Pol2 rep1
Pol2_rep1_avant2000 <- subset(Pol2_rep1_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
Pol2_rep1_apres2000 <- subset(Pol2_rep1_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
Pol2_rep1_distal <- rbind(Pol2_rep1_avant2000, Pol2_rep1_apres2000) 
head(Pol2_rep1_distal)
min(Pol2_rep1_distal$distanceToTSS)
max(Pol2_rep1_distal$distanceToTSS)
dim(Pol2_rep1_distal) # il y a 30055 peaks distaux du TSS (avant -2000 et apres +2000)
dim(Pol2_rep1_annot) # il y a en tout 48036 peaks Pol2

nb_peaks_Pol2_rep1_distal <- nrow(Pol2_rep1_distal) 
nb_peaks_Pol2_rep1_distal

nb_peaks_Pol2_rep1_NON_distal <- nrow(Pol2_rep1_annot) - nrow(Pol2_rep1_distal)
nb_peaks_Pol2_rep1_NON_distal


## Pol2 rep2
Pol2_rep2_avant2000 <- subset(Pol2_rep2_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
Pol2_rep2_apres2000 <- subset(Pol2_rep2_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
Pol2_rep2_distal <- rbind(Pol2_rep2_avant2000, Pol2_rep2_apres2000) 
head(Pol2_rep2_distal)
min(Pol2_rep2_distal$distanceToTSS)
max(Pol2_rep2_distal$distanceToTSS)
dim(Pol2_rep2_distal) # il y a 3383 peaks distaux du TSS (avant -2000 et apres +2000)
dim(Pol2_rep2_annot) # il y a en tout 11529 peaks Pol2

nb_peaks_Pol2_rep2_distal <- nrow(Pol2_rep2_distal) 
nb_peaks_Pol2_rep2_distal

nb_peaks_Pol2_rep2_NON_distal <- nrow(Pol2_rep2_annot) - nrow(Pol2_rep2_distal)
nb_peaks_Pol2_rep2_NON_distal


# On ajoute les réplicats:
nb_peaks_Pol2_distal <- nb_peaks_Pol2_rep1_distal + nb_peaks_Pol2_rep2_distal
nb_peaks_Pol2_distal
nb_peaks_Pol2_NON_distal <- nb_peaks_Pol2_rep1_NON_distal + nb_peaks_Pol2_rep2_NON_distal
nb_peaks_Pol2_NON_distal


## Chd8 rep1
Chd8_rep1_avant2000 <- subset(Chd8_rep1_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
Chd8_rep1_apres2000 <- subset(Chd8_rep1_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
Chd8_rep1_distal <- rbind(Chd8_rep1_avant2000, Chd8_rep1_apres2000) 
head(Chd8_rep1_distal)
min(Chd8_rep1_distal$distanceToTSS)
max(Chd8_rep1_distal$distanceToTSS)
dim(Chd8_rep1_distal) # il y a 8957 peaks distaux du TSS (avant -2000 et apres +2000)
dim(Chd8_rep1_annot) # il y a en tout 19942 peaks Chd8

nb_peaks_Chd8_rep1_distal <- nrow(Chd8_rep1_distal) 
nb_peaks_Chd8_rep1_distal

nb_peaks_Chd8_rep1_NON_distal <- nrow(Chd8_rep1_annot) - nrow(Chd8_rep1_distal)
nb_peaks_Chd8_rep1_NON_distal


## Chd8 rep2
Chd8_rep2_avant2000 <- subset(Chd8_rep2_annot, distanceToTSS <= -2000) # on garde les peaks avant -2000 pb du TSS
Chd8_rep2_apres2000 <- subset(Chd8_rep2_annot, distanceToTSS >= 2000) # on garde les peaks apres + 2000 pb du TSS
Chd8_rep2_distal <- rbind(Chd8_rep2_avant2000, Pol2_rep2_apres2000) 
head(Chd8_rep2_distal)
min(Chd8_rep2_distal$distanceToTSS)
max(Chd8_rep2_distal$distanceToTSS)
dim(Chd8_rep2_distal) # il y a 7081 peaks distaux du TSS (avant -2000 et apres +2000)
dim(Chd8_rep2_annot) # il y a en tout 21097 peaks Chd8

nb_peaks_Chd8_rep2_distal <- nrow(Chd8_rep2_distal) 
nb_peaks_Chd8_rep2_distal

nb_peaks_Chd8_rep2_NON_distal <- nrow(Chd8_rep2_annot) - nrow(Chd8_rep2_distal)
nb_peaks_Chd8_rep2_NON_distal


# On ajoute les réplicats:
nb_peaks_Chd8_distal <- nb_peaks_Chd8_rep1_distal + nb_peaks_Chd8_rep2_distal
nb_peaks_Chd8_distal
nb_peaks_Chd8_NON_distal <- nb_peaks_Chd8_rep1_NON_distal + nb_peaks_Chd8_rep2_NON_distal
nb_peaks_Chd8_NON_distal


## Table de contingence pour Oct4 pour les elements distaux:
peaks_distaux_autres <- nb_peaks_CTCF_distal + nb_peaks_Chd8_distal + nb_peaks_TBP_distal + nb_peaks_Pol2_distal
peaks_distaux_autres
peaks_distaux_Oct4_vs_autres <-c(nb_peaks_Oct4_distal, peaks_distaux_autres)
peaks_distaux_Oct4_vs_autres

nb_peaks_NON_distal_autre <- nb_peaks_CTCF_NON_distal + nb_peaks_Chd8_NON_distal + nb_peaks_TBP_NON_distal + nb_peaks_Pol2_NON_distal
nb_peaks_NON_distal_autre
peaks_NON_distaux_Oct4_vs_autres <-c(nb_peaks_Oct4_NON_distal, nb_peaks_NON_distal_autre)
peaks_NON_distaux_Oct4_vs_autres


contingence_table_Oct4 <- data.frame(peaks_distaux_Oct4_vs_autres, peaks_NON_distaux_Oct4_vs_autres, stringsAsFactors = FALSE)
contingence_table_Oct4

colnames(contingence_table_Oct4) <- c("Nombre_peaks_distaux", "Nombre_peaks_NON_distaux")
contingence_table_Oct4

rownames(contingence_table_Oct4) <- c("Oct4", "autre")
contingence_table_Oct4



## Table de contingence pour CTCF pour les elements distaux:
peaks_distaux_autres <- nb_peaks_Oct4_distal + nb_peaks_Chd8_distal + nb_peaks_TBP_distal + nb_peaks_Pol2_distal
peaks_distaux_autres
peaks_distaux_CTCF_vs_autres <-c(nb_peaks_CTCF_distal, peaks_distaux_autres)
peaks_distaux_CTCF_vs_autres

nb_peaks_NON_distal_autre <- nb_peaks_Oct4_NON_distal + nb_peaks_Chd8_NON_distal + nb_peaks_TBP_NON_distal + nb_peaks_Pol2_NON_distal
nb_peaks_NON_distal_autre
peaks_NON_distaux_CTCF_vs_autres <-c(nb_peaks_CTCF_NON_distal, nb_peaks_NON_distal_autre)
peaks_NON_distaux_CTCF_vs_autres


contingence_table_CTCF <- data.frame(peaks_distaux_CTCF_vs_autres, peaks_NON_distaux_CTCF_vs_autres, stringsAsFactors = FALSE)
contingence_table_CTCF

colnames(contingence_table_CTCF) <- c("Nombre_peaks_distaux", "Nombre_peaks_NON_distaux")
contingence_table_CTCF

rownames(contingence_table_CTCF) <- c("CTCF", "autre")
contingence_table_CTCF


## Table de contingence pour TBP pour les elements distaux:
peaks_distaux_autres <- nb_peaks_Oct4_distal + nb_peaks_Chd8_distal + nb_peaks_CTCF_distal + nb_peaks_Pol2_distal
peaks_distaux_autres
peaks_distaux_TBP_vs_autres <-c(nb_peaks_TBP_distal, peaks_distaux_autres)
peaks_distaux_TBP_vs_autres

nb_peaks_NON_distal_autre <- nb_peaks_Oct4_NON_distal + nb_peaks_Chd8_NON_distal + nb_peaks_CTCF_NON_distal + nb_peaks_Pol2_NON_distal
nb_peaks_NON_distal_autre
peaks_NON_distaux_TBP_vs_autres <-c(nb_peaks_TBP_NON_distal, nb_peaks_NON_distal_autre)
peaks_NON_distaux_TBP_vs_autres


contingence_table_TBP <- data.frame(peaks_distaux_TBP_vs_autres, peaks_NON_distaux_TBP_vs_autres, stringsAsFactors = FALSE)
contingence_table_TBP

colnames(contingence_table_TBP) <- c("Nombre_peaks_distaux", "Nombre_peaks_NON_distaux")
contingence_table_TBP

rownames(contingence_table_TBP) <- c("TBP", "autre")
contingence_table_TBP


## Table de contingence pour Pol2 pour les elements distaux:
peaks_distaux_autres <- nb_peaks_Oct4_distal + nb_peaks_Chd8_distal + nb_peaks_CTCF_distal + nb_peaks_TBP_distal
peaks_distaux_autres
peaks_distaux_Pol2_vs_autres <-c(nb_peaks_Pol2_distal, peaks_distaux_autres)
peaks_distaux_Pol2_vs_autres

nb_peaks_NON_distal_autre <- nb_peaks_Oct4_NON_distal + nb_peaks_Chd8_NON_distal + nb_peaks_CTCF_NON_distal + nb_peaks_TBP_NON_distal
nb_peaks_NON_distal_autre
peaks_NON_distaux_Pol2_vs_autres <-c(nb_peaks_Pol2_NON_distal, nb_peaks_NON_distal_autre)
peaks_NON_distaux_Pol2_vs_autres


contingence_table_Pol2 <- data.frame(peaks_distaux_Pol2_vs_autres, peaks_NON_distaux_Pol2_vs_autres, stringsAsFactors = FALSE)
contingence_table_Pol2

colnames(contingence_table_Pol2) <- c("Nombre_peaks_distaux", "Nombre_peaks_NON_distaux")
contingence_table_Pol2

rownames(contingence_table_Pol2) <- c("Pol2", "autre")
contingence_table_Pol2


## Table de contingence pour Chd8 pour les elements distaux:
peaks_distaux_autres <- nb_peaks_Oct4_distal + nb_peaks_Pol2_distal + nb_peaks_CTCF_distal + nb_peaks_TBP_distal
peaks_distaux_autres
peaks_distaux_Chd8_vs_autres <-c(nb_peaks_Chd8_distal, peaks_distaux_autres)
peaks_distaux_Chd8_vs_autres

nb_peaks_NON_distal_autre <- nb_peaks_Oct4_NON_distal + nb_peaks_Pol2_NON_distal + nb_peaks_CTCF_NON_distal + nb_peaks_TBP_NON_distal
nb_peaks_NON_distal_autre
peaks_NON_distaux_Chd8_vs_autres <-c(nb_peaks_Chd8_NON_distal, nb_peaks_NON_distal_autre)
peaks_NON_distaux_Chd8_vs_autres


contingence_table_Chd8 <- data.frame(peaks_distaux_Chd8_vs_autres, peaks_NON_distaux_Chd8_vs_autres, stringsAsFactors = FALSE)
contingence_table_Chd8

colnames(contingence_table_Chd8) <- c("Nombre_peaks_distaux", "Nombre_peaks_NON_distaux")
contingence_table_Chd8

rownames(contingence_table_Chd8) <- c("Chd8", "autre")
contingence_table_Chd8


## Test de fisher sur la table de contingence de Oct4
testF_Oct4 <- fisher.test(x=contingence_table_Oct4, # table de contingence
                          y = NULL, # a factor object; ignored if x is a matrix.
                          hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
                          control = list(), # a list with named components for low level algorithm control.
                          or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
                          alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
                          conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
                          conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
                          B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.



testF_CTCF <- fisher.test(x=contingence_table_CTCF, # table de contingence
                          y = NULL, # a factor object; ignored if x is a matrix.
                          hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
                          control = list(), # a list with named components for low level algorithm control.
                          or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
                          alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
                          conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
                          conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
                          B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.


testF_TBP <- fisher.test(x=contingence_table_TBP, # table de contingence
                         y = NULL, # a factor object; ignored if x is a matrix.
                         hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
                         control = list(), # a list with named components for low level algorithm control.
                         or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
                         alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
                         conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
                         conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
                         B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.


testF_Pol2 <- fisher.test(x=contingence_table_Pol2, # table de contingence
                          y = NULL, # a factor object; ignored if x is a matrix.
                          hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
                          control = list(), # a list with named components for low level algorithm control.
                          or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
                          alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
                          conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
                          conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
                          B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.


testF_Chd8 <- fisher.test(x=contingence_table_Chd8, # table de contingence
                          y = NULL, # a factor object; ignored if x is a matrix.
                          hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
                          control = list(), # a list with named components for low level algorithm control.
                          or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
                          alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
                          conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
                          conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
                          B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.


testF_Chd8
testF_Oct4
testF_TBP
testF_Pol2
testF_CTCF


#######################################################
# Région proximale : [-400 ; +100]
#######################################################


### Comptage du nombre de peaks pour tous les facteurs:

count_distance_data <- function(df, dist) {
  # Recover distance information
  min_threshold <- df$distanceToTSS >= -dist;
  max_threshold <- df$distanceToTSS <=  dist;
  distal <- min_threshold & max_threshold;
  # Compute results
  nb_distal <- table(distal)["TRUE"];
  nb_non_distal <- table(distal)["FALSE"];
  total <- length(distal);
  # Build results object
  results <- c(nb_distal, nb_non_distal, total);
  names(results) <- c("nb_distal", "nb_non_distal", "total");
  return(results)
}

### Oct4
Oct4_rep1_peaks_tot_prox <- count_distance_data(Oct4_rep1_annot, 300)
Oct4_rep1_peaks_tot_prox

Oct4_rep2_peaks_tot_prox <- count_distance_data(Oct4_rep2_annot, 300)
Oct4_rep2_peaks_tot_prox

nb_peaks_Oct4_prox <- 8772 + 7217
nb_peaks_Oct4_NON_prox <- 34076 + 26894

### CTCF
CTCF_peaks_tot_prox <- count_distance_data(CTCF_annot, 300)
CTCF_peaks_tot_prox

nb_peaks_CTCF_prox <- 4320
nb_peaks_CTCF_NON_prox <- 48440

### TBP
TBP_rep1_peaks_tot_prox <- count_distance_data(TBP_rep1_annot, 300)
TBP_rep1_peaks_tot_prox

TBP_rep2_peaks_tot_prox <- count_distance_data(TBP_rep2_annot, 300)
TBP_rep2_peaks_tot_prox

nb_peaks_TBP_prox <- 12560 + 258
nb_peaks_TBP_NON_prox <- 36888 + 485


### Pol2
Pol2_rep1_peaks_tot_prox <- count_distance_data(Pol2_rep1_annot, 300)
Pol2_rep1_peaks_tot_prox

Pol2_rep2_peaks_tot_prox <- count_distance_data(Pol2_rep2_annot, 300)
Pol2_rep2_peaks_tot_prox

nb_peaks_Pol2_prox <- 13962 + 7380
nb_peaks_Pol2_NON_prox <- 34074 + 4149


### Chd8
Chd8_rep1_peaks_tot_prox <- count_distance_data(Chd8_rep1_annot, 300)
Chd8_rep1_peaks_tot_prox

Chd8_rep2_peaks_tot_prox <- count_distance_data(Chd8_rep2_annot, 300)
Chd8_rep2_peaks_tot_prox

nb_peaks_Chd8_prox <- 9844 + 8877
nb_peaks_Chd8_NON_prox <- 10098 + 12220


### Tables de contingence:

## Table de contingence pour Oct4 pour les elements proximaux:
peaks_prox_autres <- nb_peaks_CTCF_prox + nb_peaks_Chd8_prox + nb_peaks_TBP_prox + nb_peaks_Pol2_prox
peaks_prox_autres
peaks_prox_Oct4_vs_autres <-c(nb_peaks_Oct4_prox, peaks_prox_autres)
peaks_prox_Oct4_vs_autres

nb_peaks_NON_prox_autre <- nb_peaks_CTCF_NON_prox + nb_peaks_Chd8_NON_prox + nb_peaks_TBP_NON_prox + nb_peaks_Pol2_NON_prox
nb_peaks_NON_prox_autre
peaks_NON_prox_Oct4_vs_autres <-c(nb_peaks_Oct4_NON_prox, nb_peaks_NON_prox_autre)
peaks_NON_prox_Oct4_vs_autres


contingence_table_Oct4_prox <- data.frame(peaks_prox_Oct4_vs_autres, peaks_NON_prox_Oct4_vs_autres, stringsAsFactors = FALSE)
contingence_table_Oct4_prox

colnames(contingence_table_Oct4_prox) <- c("Nombre_peaks_distaux", "Nombre_peaks_NON_distaux")
rownames(contingence_table_Oct4_prox) <- c("Oct4", "autre")
contingence_table_Oct4_prox



## Table de contingence pour CTCF pour les elements distaux:
peaks_prox_autres <- nb_peaks_Oct4_prox + nb_peaks_Chd8_prox + nb_peaks_TBP_prox + nb_peaks_Pol2_prox
peaks_prox_autres
peaks_prox_CTCF_vs_autres <-c(nb_peaks_CTCF_prox, peaks_prox_autres)
peaks_prox_CTCF_vs_autres

nb_peaks_NON_prox_autre <- nb_peaks_Oct4_NON_prox + nb_peaks_Chd8_NON_prox + nb_peaks_TBP_NON_prox + nb_peaks_Pol2_NON_prox
nb_peaks_NON_prox_autre
peaks_NON_prox_CTCF_vs_autres <-c(nb_peaks_CTCF_NON_prox, nb_peaks_NON_prox_autre)
peaks_NON_prox_CTCF_vs_autres


contingence_table_CTCF_prox <- data.frame(peaks_prox_CTCF_vs_autres, peaks_NON_prox_CTCF_vs_autres, stringsAsFactors = FALSE)
contingence_table_CTCF_prox

colnames(contingence_table_CTCF_prox) <- c("Nombre_peaks_distaux", "Nombre_peaks_NON_distaux")
rownames(contingence_table_CTCF_prox) <- c("CTCF", "autre")
contingence_table_CTCF_prox


## Table de contingence pour TBP pour les elements distaux:
peaks_prox_autres <- nb_peaks_Oct4_prox + nb_peaks_Chd8_prox + nb_peaks_CTCF_prox + nb_peaks_Pol2_prox
peaks_prox_autres
peaks_prox_TBP_vs_autres <-c(nb_peaks_TBP_prox, peaks_prox_autres)
peaks_prox_TBP_vs_autres

nb_peaks_NON_prox_autre <- nb_peaks_Oct4_NON_prox + nb_peaks_Chd8_NON_prox + nb_peaks_CTCF_NON_prox + nb_peaks_Pol2_NON_prox
nb_peaks_NON_prox_autre
peaks_NON_prox_TBP_vs_autres <-c(nb_peaks_TBP_NON_prox, nb_peaks_NON_prox_autre)
peaks_NON_prox_TBP_vs_autres


contingence_table_TBP_prox <- data.frame(peaks_prox_TBP_vs_autres, peaks_NON_prox_TBP_vs_autres, stringsAsFactors = FALSE)
contingence_table_TBP_prox

colnames(contingence_table_TBP_prox) <- c("Nombre_peaks_distaux", "Nombre_peaks_NON_distaux")
rownames(contingence_table_TBP_prox) <- c("TBP", "autre")
contingence_table_TBP_prox


## Table de contingence pour Pol2 pour les elements distaux:
peaks_prox_autres <- nb_peaks_Oct4_prox + nb_peaks_Chd8_prox + nb_peaks_CTCF_prox + nb_peaks_TBP_prox
peaks_prox_autres
peaks_prox_Pol2_vs_autres <-c(nb_peaks_Pol2_prox, peaks_prox_autres)
peaks_prox_Pol2_vs_autres

nb_peaks_NON_prox_autre <- nb_peaks_Oct4_NON_prox + nb_peaks_Chd8_NON_prox + nb_peaks_CTCF_NON_prox + nb_peaks_TBP_NON_prox
nb_peaks_NON_prox_autre
peaks_NON_prox_Pol2_vs_autres <-c(nb_peaks_Pol2_NON_prox, nb_peaks_NON_prox_autre)
peaks_NON_prox_Pol2_vs_autres


contingence_table_Pol2_prox <- data.frame(peaks_prox_Pol2_vs_autres, peaks_NON_prox_Pol2_vs_autres, stringsAsFactors = FALSE)
contingence_table_Pol2_prox

colnames(contingence_table_Pol2_prox) <- c("Nombre_peaks_distaux", "Nombre_peaks_NON_distaux")
rownames(contingence_table_Pol2_prox) <- c("Pol2", "autre")
contingence_table_Pol2_prox


## Table de contingence pour Chd8 pour les elements distaux:
peaks_prox_autres <- nb_peaks_Oct4_prox + nb_peaks_Pol2_prox + nb_peaks_CTCF_prox + nb_peaks_TBP_prox
peaks_prox_autres
peaks_prox_Chd8_vs_autres <-c(nb_peaks_Chd8_prox, peaks_prox_autres)
peaks_prox_Chd8_vs_autres

nb_peaks_NON_prox_autre <- nb_peaks_Oct4_NON_prox + nb_peaks_Pol2_NON_prox + nb_peaks_CTCF_NON_prox + nb_peaks_TBP_NON_prox
nb_peaks_NON_prox_autre
peaks_NON_prox_Chd8_vs_autres <-c(nb_peaks_Chd8_NON_prox, nb_peaks_NON_prox_autre)
peaks_NON_prox_Chd8_vs_autres


contingence_table_Chd8_prox <- data.frame(peaks_prox_Chd8_vs_autres, peaks_NON_prox_Chd8_vs_autres, stringsAsFactors = FALSE)
contingence_table_Chd8_prox

colnames(contingence_table_Chd8_prox) <- c("Nombre_peaks_distaux", "Nombre_peaks_NON_distaux")
rownames(contingence_table_Chd8_prox) <- c("Chd8", "autre")
contingence_table_Chd8_prox


## Test de fisher sur la table de contingence de Oct4
testF_Oct4_prox <- fisher.test(x=contingence_table_Oct4_prox, # table de contingence
                          y = NULL, # a factor object; ignored if x is a matrix.
                          hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
                          control = list(), # a list with named components for low level algorithm control.
                          or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
                          alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
                          conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
                          conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
                          B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.



testF_CTCF_prox <- fisher.test(x=contingence_table_CTCF_prox, # table de contingence
                          y = NULL, # a factor object; ignored if x is a matrix.
                          hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
                          control = list(), # a list with named components for low level algorithm control.
                          or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
                          alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
                          conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
                          conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
                          B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.


testF_TBP_prox <- fisher.test(x=contingence_table_TBP_prox, # table de contingence
                         y = NULL, # a factor object; ignored if x is a matrix.
                         hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
                         control = list(), # a list with named components for low level algorithm control.
                         or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
                         alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
                         conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
                         conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
                         B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.


testF_Pol2_prox <- fisher.test(x=contingence_table_Pol2_prox, # table de contingence
                          y = NULL, # a factor object; ignored if x is a matrix.
                          hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
                          control = list(), # a list with named components for low level algorithm control.
                          or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
                          alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
                          conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
                          conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
                          B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.


testF_Chd8_prox <- fisher.test(x=contingence_table_Chd8_prox, # table de contingence
                          y = NULL, # a factor object; ignored if x is a matrix.
                          hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
                          control = list(), # a list with named components for low level algorithm control.
                          or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
                          alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
                          conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
                          conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
                          B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.


testF_Chd8_prox
testF_Oct4_prox
testF_TBP_prox
testF_Pol2_prox
testF_CTCF_prox







##### Writing annotations to file #####

##### Pour Oct4 ######
## Get unique entrez gene Ids (Id uniques)
entrezids <- unique(Oct4_annot$geneId)

## Get mm10 entrez to gene symbol mappings (obtenir les symboles de gene pour mm10)
entrez2gene <- grcm38 %>% filter(entrez %in% entrezids) %>% dplyr::select(entrez, symbol)

## Match to each annotation dataframe (matcher les dataframe d'annotation)
m <- match(Oct4_annot$geneId, entrez2gene$entrez)
Oct4_annot_modif<- cbind(Oct4_annot[,1:14], geneSymbol=entrez2gene$symbol[m], Oct4_annot[,15:16])

write.table(Oct4_annot, file="./Oct4_annotation.txt", sep="\t", quote=F, row.names=F)


##### Pour CTCF ######
# Get unique entrez gene Ids
entrezids <- unique(CTCF_annot$geneId)

## Get mm10 entrez to gene symbol mappings #### A modifier pour la souris ####
entrez2gene <- grcm38 %>% filter(entrez %in% entrezids) %>% dplyr::select(entrez, symbol)

## Match to each annotation dataframe
m <- match(CTCF_annot$geneId, entrez2gene$entrez)
CTCF_annot_modif <- cbind(CTCF_annot[,1:14], geneSymbol=entrez2gene$symbol[m], CTCF_annot[,15:16])

write.table(CTCF_annot_modif, file="./CTCF_annotation.txt", sep="\t", quote=F, row.names=F)


