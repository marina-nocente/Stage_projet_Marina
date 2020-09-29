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


##### Loading des data #####
# As input we need to provide the names of our BED files in a list format.

samplefiles <- list.files("/home/mnocente/Bureau/Stage_projet_Marina/Projet_ChIP_Marina/results/macs2_peaks", pattern= ".narrowPeak", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("CTCF", "Oct4") #### Mettre les noms dans le bon ordre !!!! ####


##### Assign annotation db #####

# We need to assign annotation databases generated from UCSC to a variable:

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


# ChIPseeker et sa fonction annotatePeak permettent d'annoter les pics avec le gène et la région génomique les plus proches de là où se trouve le pic.
# The annotatePeak function by default uses the TSS method, and provides parameters to specify a max distance cutoff.


##### Get annotations #####
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)

peakAnnoList # annotation information is stored in the peakAnnoList

##### Visualisation des annotations #####

# Pie chart of genomic region annotation
plotAnnoPie(peakAnnoList[["Oct4"]])
#plotAnnoPie(peakAnnoList[["TBP"]])
#plotAnnoPie(peakAnnoList[["Pol2"]])
plotAnnoPie(peakAnnoList[["CTCF"]])
#plotAnnoPie(peakAnnoList[["Chd8"]])

# Barchart (multiple samples for comparison)
plotAnnoBar(peakAnnoList)

# Distribution of TF-binding loci relative to TSS
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")


##### Writing annotations to file #####

### Get annotation data frame
Oct4_annot <- as.data.frame(peakAnnoList[["Oct4"]]@anno)
CTCF_annot <- as.data.frame(peakAnnoList[["CTCF"]]@anno)

# Dans ces dataframes, on doit voir des colonnes correspondant au fichier BED d'entrée et des colonnes supplémentaires contenant le ou les gènes les plus proches, la distance entre le pic et le TSS du gène le plus proche, la région génomique du pic et d'autres informations. 
# Ordre de priorite : Promoter, 5’ UTR, 3’ UTR, Exon, Intron, Downstream (defined as the downstream of gene end), Intergenic.
# Il manque juste les symboles de gènes répertoriés dans le tableau.

# Get unique entrez gene Ids
entrezids <- unique(Oct4_annot$geneId)

# Get mm10 entrez to gene symbol mappings #### A modifier pour la souris ####
entrez2gene <- grch37 %>% 
  filter(entrez %in% entrezids) %>% 
  dplyr::select(entrez, symbol)

# Match to each annotation dataframe
m <- match(Oct4_annot$geneId, entrez2gene$entrez)
Oct4_annot <- cbind(Oct4_annot[,1:14], geneSymbol=entrez2gene$symbol[m], Oct4_annot[,15:16])

write.table(Oct4_annot, file="results/Oct4_annotation.txt", sep="\t", quote=F, row.names=F)


##### Functional enrichment analysis #####

### Single sample analysis ###
# La liste de gènes des annotations Oct4 en entree pour une analyse d'enrichissement GO.

# Run GO enrichment analysis 
entrezids <- Oct4_annot$geneId %>% 
  as.character() %>% 
  unique()


ego <- enrichGO(gene = entrezids, 
                keyType = "ENTREZID", 
                OrgDb = org.Mm.eg.db, 
                ont = "BP", # biological process
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "results/clusterProfiler_Oct4.csv")

# Dotplot visualization
dotplot(ego, showCategory=50)


### Multiple samples ###

# Notre ensemble de données se compose de peaks de plusieurs facteurs différents, il serait donc utile de comparer les résultats d'enrichissement fonctionnel.

# Create a list with genes from each sample
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

# Run KEGG analysis
compKEGG <- compareCluster(geneCluster = genes, 
                           fun = "enrichKEGG",
                           organism = "mouse",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")












