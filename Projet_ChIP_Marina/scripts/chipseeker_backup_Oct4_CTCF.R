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

samplefiles <- list.files("results/macs2_peaks/just_chr_{sample_wildcard}.narrowPeak", full.names = T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("CTCF", "Oct4") #### Mettre les noms dans le bon ordre !!!! ####

print("les fichiers ont bien été chargés")
head(samplefiles)

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

## Pie chart of genomic region annotation
png(file = "plot_camembert_annotation_Oct4.png")
plotAnnoPie(peakAnnoList[["Oct4"]])
dev.off()

#plotAnnoPie(peakAnnoList[["TBP"]])
#plotAnnoPie(peakAnnoList[["Pol2"]])

png(file = "plot_camembert_annotation_CTCF.png")
plotAnnoPie(peakAnnoList[["CTCF"]])
dev.off()

#plotAnnoPie(peakAnnoList[["Chd8"]])

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
CTCF_annot <- as.data.frame(peakAnnoList[["CTCF"]]@anno)

# dans ce dataframe, la derniere colonne correspond à la distance au TSS.
# Dans ces dataframes, on doit voir des colonnes correspondant au fichier BED d'entrée et des colonnes supplémentaires contenant le ou les gènes les plus proches, la distance entre le pic et le TSS du gène le plus proche, la région génomique du pic et d'autres informations. 
# Ordre de priorite : Promoter, 5’ UTR, 3’ UTR, Exon, Intron, Downstream (defined as the downstream of gene end), Intergenic.
# Il manque juste les symboles de gènes répertoriés dans le tableau.


### Creation de la table de contingence entre 0 et 1 kb:
## Oct4
Oct4_zero_un_kb <- subset(Oct4_annot, 0 <= distanceToTSS) # on garde les peaks a plus de 0kb du TSS
Oct4_zero_un_kb <- subset(Oct4_zero_un_kb, distanceToTSS <=1000) # on garde les peaks entre 0 et 1 kb du TSS
head(Oct4_zero_un_kb)
min(Oct4_zero_un_kb$distanceToTSS)
max(Oct4_zero_un_kb$distanceToTSS)
dim(Oct4_zero_un_kb) # il y a 8477 peaks compris entre 0 et 1 kb
dim(Oct4_annot) # il y a en tout 42848 peaks Oct4

peaks_Oct4_zero_un_kb <- nrow(Oct4_zero_un_kb) # nombre de peaks de Oct4 localise entre 0 et 1 kb
peaks_Oct4_zero_un_kb

peaks_Oct4_PAS_zero_un_kb <- nrow(Oct4_annot) - nrow(Oct4_zero_un_kb) # nombre de peaks de Oct4 qui ne sont PAS localise entre 0 et 1 kb
peaks_Oct4_PAS_zero_un_kb


## CTCF
CTCF_zero_un_kb <- subset(CTCF_annot, 0 <= distanceToTSS) # on garde les peaks a plus de 0kb du TSS
CTCF_zero_un_kb <- subset(CTCF_zero_un_kb, distanceToTSS <=1000) # on garde les peaks entre 0 et 1 kb du TSS
head(CTCF_zero_un_kb)
min(CTCF_zero_un_kb$distanceToTSS)
max(CTCF_zero_un_kb$distanceToTSS)
dim(CTCF_zero_un_kb) # il y a 4550 peaks compris entre 0 et 1 kb
dim(CTCF_annot) # il y a en tout 52760 peaks CTCF

peaks_CTCF_zero_un_kb <- nrow(CTCF_zero_un_kb) # nombre de peaks de CTCF localise entre 0 et 1 kb
peaks_CTCF_zero_un_kb

peaks_CTCF_PAS_zero_un_kb <- nrow(CTCF_annot) - nrow(CTCF_zero_un_kb) # nombre de peaks de Oct4 qui ne sont PAS localise entre 0 et 1 kb
peaks_CTCF_PAS_zero_un_kb



## Table de contingence pour Oct4 pour 0 à 1 kb:
peaks_zero_un_kb_Oct4_vs_autre <- c(peaks_Oct4_zero_un_kb, peaks_CTCF_zero_un_kb)
peaks_PAS_zero_un_kb_Oct4_vs_autre <- c(peaks_Oct4_PAS_zero_un_kb, peaks_CTCF_PAS_zero_un_kb)

contingence_table_Oct4 <- data.frame(peaks_zero_un_kb_Oct4_vs_autre, peaks_PAS_zero_un_kb_Oct4_vs_autre, stringsAsFactors = FALSE)
contingence_table_Oct4

colnames(contingence_table_Oct4) <- c("Nombre_peaks_entre_0-1_kb", "Nombre_peaks_PAS_entre_0-1_kb")
contingence_table_Oct4

rownames(contingence_table_Oct4) <- c("Oct4", "autre")
contingence_table_Oct4

## Table de contingence pour CTCF pour 0 à 1 kb:
peaks_zero_un_kb_CTCF_vs_autre <- c(peaks_CTCF_zero_un_kb, peaks_Oct4_zero_un_kb)
peaks_PAS_zero_un_kb_CTCF_vs_autre <- c(peaks_CTCF_PAS_zero_un_kb, peaks_Oct4_PAS_zero_un_kb)

contingence_table_CTCF <- data.frame(peaks_zero_un_kb_CTCF_vs_autre, peaks_PAS_zero_un_kb_CTCF_vs_autre, stringsAsFactors = FALSE)
contingence_table_CTCF

colnames(contingence_table_CTCF) <- c("Nombre_peaks_entre_0-1_kb", "Nombre_peaks_PAS_entre_0-1_kb")
contingence_table_CTCF

rownames(contingence_table_CTCF) <- c("CTCF", "autre")
contingence_table_CTCF



## Test de fisher sur la table de contingence de Oct4
fisher.test(x=contingence_table_Oct4, # table de contingence
            y = NULL, # a factor object; ignored if x is a matrix.
            hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
            control = list(), # a list with named components for low level algorithm control.
            or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
            alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
            conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
            conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
            B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.


#data:  contingence_table_Oct4
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval: 2.513752 2.716552
#sample estimates:
#odds ratio : 2.613199 

fisher.test(x=contingence_table_CTCF, # table de contingence
            y = NULL, # a factor object; ignored if x is a matrix.
            hybridPars = c(expect = 5, percent = 80, Emin = 1), # a numeric vector of length 3, by default describing “Cochran's conditions” for the validity of the chisquare approximation, see ‘Details’.
            control = list(), # a list with named components for low level algorithm control.
            or = 1, # the hypothesized odds ratio. Only used in the 2×2 case.
            alternative = "two.sided",  # indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2×2 case.
            conf.int = TRUE, # logical indicating if a confidence interval for the odds ratio in a 2×2 table should be computed (and returned).
            conf.level = 0.95, # confidence level for the returned confidence interval. Only used in the 2×2 case and if conf.int = TRUE.
            B = 2000) #an integer specifying the number of replicates used in the Monte Carlo test.


#data:  contingence_table_CTCF
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval: 0.3681137 0.3978116
#sample estimates:
#odds ratio : 0.3826728 





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



##### Functional enrichment analysis #####

#### Single sample analysis ####

### La liste de gènes des annotations Oct4 en entree pour une analyse d'enrichissement GO.

## Run GO enrichment analysis 
entrezids_Oct4 <- Oct4_annot$geneId %>% as.character() %>% unique()

ego_Oct4 <- enrichGO(gene = entrezids_Oct4, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_Oct4 <- data.frame(ego_Oct4)
write.csv(cluster_summary_Oct4, "results/clusterProfiler_Oct4.csv")

## Dotplot visualization
png(file = "dotplotfunctional_enrichment_Oct4.png")
dotplot(ego_Oct4, showCategory=50)
dev.off()



### La liste de gènes des annotations CTCF en entree pour une analyse d'enrichissement GO.

## Run GO enrichment analysis 
entrezids_CTCF <- CTCF_annot$geneId %>% as.character() %>% unique()

ego_CTCF <- enrichGO(gene = entrezids_CTCF, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_CTCF <- data.frame(ego_CTCF)
write.csv(cluster_summary_CTCF, "results/clusterProfiler_CTCF.csv")

## Dotplot visualization
png(file = "dotplotfunctional_enrichment_CTCF.png")
dotplot(ego_CTCF, showCategory=50)
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

png(file = "KEGG_Pathway_Enrichment_Analysis.png")
dotplot(compKEGG, showCategory = 20, title = "KEGG_Pathway_Enrichment_Analysis")
dev.off()











