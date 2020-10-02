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
samplefiles <- list.files("/shared/home/mnocente/macs2", pattern= ".narrowPeak",full.names = T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("Chd8_rep1", "Chd8_rep2", "CTCF", "Oct4_rep1", "Oct4_rep2", "Pol2_rep1", "Pol2_rep2", "TBP_rep1", "TBP_rep2")

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
Oct4_rep1_annot <- as.data.frame(peakAnnoList[["Oct4_rep1"]]@anno)
Oct4_rep2_annot <- as.data.frame(peakAnnoList[["Oct4_rep2"]]@anno)

TBP_rep1_annot <- as.data.frame(peakAnnoList[["TBP_rep1"]]@anno)
TBP_rep2_annot <- as.data.frame(peakAnnoList[["TBP_rep2"]]@anno)

CTCF_annot <- as.data.frame(peakAnnoList[["CTCF"]]@anno)

Pol2_rep1_annot <- as.data.frame(peakAnnoList[["Pol2_rep1"]]@anno)
Pol2_rep2_annot <- as.data.frame(peakAnnoList[["Pol2_rep2"]]@anno)

Chd8_rep1_annot <- as.data.frame(peakAnnoList[["Chd8_rep1"]]@anno)
Chd8_rep2_annot <- as.data.frame(peakAnnoList[["Chd8_rep2"]]@anno)

# dans ce dataframe, la derniere colonne correspond à la distance au TSS.
# Dans ces dataframes, on doit voir des colonnes correspondant au fichier BED d'entrée et des colonnes supplémentaires contenant le ou les gènes les plus proches, la distance entre le pic et le TSS du gène le plus proche, la région génomique du pic et d'autres informations. 
# Ordre de priorite : Promoter, 5’ UTR, 3’ UTR, Exon, Intron, Downstream (defined as the downstream of gene end), Intergenic.
# Il manque juste les symboles de gènes répertoriés dans le tableau.


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

a <- data.frame(c("Oct4_rep1_annot", "Oct4_rep2_annot"))
a$distanceToTSS <- c(100, 24)
a
b <- count_distance_data(Oct4_rep1_annot, 1000)
print(b)





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


### Creation de la table de contingence pour les éléments proxiamux du TSS:
## Oct4 rep1
Oct4_rep1_avant400 <- subset(Oct4_rep1_annot, distanceToTSS >= -400) # on garde les peaks avant -400 pb du TSS
Oct4_rep1_avant100 <- subset(Oct4_rep1_avant400, distanceToTSS <= 100) # on garde les peaks apres + 100 pb du TSS







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











