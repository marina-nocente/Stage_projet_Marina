# Stage_projet_Marina : contenu

Ce projet GitHub contient un répertoire "Projet_ChIP_Marina" qui contient un répertoire "scripts", un répertoire "fichiers_test", un fichier "Snakefile" et son fichier de configuration "config.yaml", ainsi que les fichiers .rst qui vont avec le Snakefile (pour le rapport généré avec snakemake).


Le répertoire "script" contient le script "fonction_parler_cluster" qui permet d'appeler un script avec des variables et des paramètres sur un cluster Slurm ou Torque. 

Le répertoire "script" contient également lui-même 11 répertoires (un pour chaque outil et étape): 
- "fastQC" (observation ds qualité des reads)
- "fastqScreen" (observation d'éventuelles contaminations)
- "fastp" (trimming des reads)
- "bowtie2" (alignement des reads sur le génome de la souris mm10)
- "samtools_flagstat" (statistiques sur la qualité de l'alignement)
- "samtools_sort" (les reads sont triés par coordonnées de chromosome)
- "Picard" (marquages des reads dupliqués)
- "samtools_view_sort_fichiers_nettoyes" (élimination ds reads dupliqués)
- "samtools_index_bam" (indexation des bam pour la visualisation dans IGV)


Le projet GitHub contient également le fichier "feuille_route.ods" qui répertorie toute les informations connues de mes échantillons, les différentes étapes de l'analyse à faire et leur avancement. La feuille de route a été mise à jour avec le nom des échantillons de la plateforme de séquencage et le nom avec lequel je les ai renommés.


# Creation d'un Snakefile utilisable avec snakemake:

## Installation de Snakemake via le fichier yaml de Thibault (bonnes versions, outils etc...)

```{bash, eval=FALSE}
conda env create -f ./snakemake_env.yaml
```

## Quelques mots sur snakemake et un snakefile:

Un workflow Snakemake est défini en spécifiant des règles dans un Snakefile. Les règles décomposent le workflow en petites étapes (par exemple, une règle par outil) en spécifiant comment créer des ensembles de fichiers de sortie à partir d'ensembles de fichiers d'entrée. Snakemake détermine automatiquement les dépendances entre les règles en faisant correspondre les noms de fichiers.
Le langage Snakemake est basé sur du Python. Une règle Snakemake a un nom (par ex fastp) et un certain nombre de directives (par ex: input, output et shell). Les directives d'entrée et de sortie sont suivies de listes de fichiers qui devraient être utilisés ou créés par la règle. On peut spécifier les paramètres de l'outil utilisé dans la règle (directement ou en utilisant un fichier de config).
On utilise des wildcards pour avoir des noms génériques de nos fichiers.

Nous allons utiliser des wrappers. Les wrappers snakemake sont une collection de wrappers réutilisables qui permettent d'utiliser rapidement les règles et workflows de snakemake. Les wrappers gèrent de facon autonome les environnements conda, les échantillons, les variables etc...


## Utilisation des wrappers

J'utilise les wrappers dont Thibault a rendu disponible le code.
[GitHub de Thibault pour acceder a ses wrappers](https://github.com/tdayris-perso/snakemake-wrappers)
[Ex de la doc du wrapper de Picard markduplicates](file:///home/mnocente/snakemake-wrappers/docs/_build/html/wrappers/picard/markduplicates.html)


```{bash, eval=FALSE}
git clone https://github.com/tdayris/snakemake-wrappers.git
git checkout Unofficial
cd snakemake-wrappers/docs/_build/html
firefox index.html

```

# Les differentes etapes de mon Snakefile:

## Etape 1 : vérification de la qualité des reads séquencés

L'outil fastQC permet de visualiser la qualite et le contenu des reads apres sequencage.
Il prend en entrée un fichier .fastq.gz et il fournit pour chaque read (R1 et R2) un fichier .zip et un fichier .html.

Dans le fichier de resultat de fastQC, on peut ainsi acceder a differentes rubriques comme :
- les statistiques basiques (nombre total de reads, longueur des reads, pourcentage en GC...)
- la qualite de la sequence pour chaque base (on a un score de qualite pour chaque base, c'est a dire le risque de se tromper)
- le contenu en A/T/G/C le long du read.
- le contenu en GC sur la sequence totale du read. Un contenu en GC "normal" pour la souris est entre 40 et 50%.
- le contenu en base N des reads (ce sont les bases que le sequenceur n'a pas reussi a identifier).
- le niveau de duplicat.
- les séquences surreprésentées.
- la contamination par les adaptateurs.

Les differentes rubriques du rapport fastQC aident a prendre des decisions pour le trimming des fichiers fastq.gz.

[doc fastQC](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf)

## Etape 2 : trimming des reads

Les fichiers input sont les fichiers fastq.gz fournis par la plateforme de sequencage et dont la qualite a ete evaluee par fastQC.
L'outil fastp permet de nettoyer les reads des adaptateurs et si besoin de couper une partie des reads de mauvaise qualite.

Voici quelques caracteristiques de l'outil fastp (d'apres sa documentation):
-	suit la qualite des donnees avant et apres trimming (courbes de qualité, contenu de base, KMER, Q20 / Q30, rapport GC, duplication, contenu en adaptateur ...).
-	filtre les mauvaises reads (qualité trop faible, trop courte ou trop de N ...)
-	coupe des bases de faible qualité en 5 'et 3' des reads en évaluant la qualité moyenne d'une fenêtre coulissante.
-	elimine automatiquement les adaptateurs.
-	fournit un rapport des resultats au format JSON ou HTML (comme fastQC).

On obtient ainsi des reads nettoyes au format fastq.gz et un rapport de ce que fastp a fait au format .html ou .json.

[manuel options fastp](https://manpages.debian.org/testing/fastp/fastp.1.en.html)
[github fastp](https://github.com/OpenGene/fastp#features)


## Etape 3 : Alignement des reads sur le génome de référence de la souris

### Indexation du génome mm10

J'ai pris le FASTA du genome de la souris sur ENSEMBL: "contenu Genome sequence, primary assembly (GRCm38) mm10, region PRI. Nucleotide sequence of the GRCm38 primary genome assemby (chromosomes and scaffolds (suffisant pour ce que l'on veut nous, pas besoin du reste)). The sequence region names are in the same as in the GTF/GFF file".

J'ai utilisé bowtie2 pour indexé le génome.

### Mapping des reads séquencés et nettoyés sur le génome

Les reads nettoyes par fastp sont mappés sur le genome indexe (par bowtie2) de la souris (version mm10).
Bowtie utilise donc entree des fichiers .fastq.gz et un fichier fasta indexé.

D'apres sa documentation, Bowtie 2 est un outil ultrarapide et efficace en mémoire pour aligner les reads de séquençage sur de longues séquences de référence (comme les genomes de mammifere).
Bowtie 2 peut faire des alignements avec des gaps, des alignements locaux et paired-end.
Par défaut, Bowtie 2 effectue un alignement de read "de bout en bout" (= end-to-end). C'est-à-dire qu'il recherche les alignements impliquant toutes les bases du read.
Lorsque l'option –local est spécifiée, Bowtie 2 effectue un alignement local du read. Dans ce mode, Bowtie 2 peut "découper" ou "couper" certaines bases du read à l'une ou aux deux extrémités de l'alignement si cela maximise le score d'alignement.

Bowtie2 a tourne avec ses parametres par defaut (alignement end-to-end) et j'ai precise le parametre --align-paired-reads pour lui indiquer que les reads sont paires.
Bowtie2 mappe les reads paires sur le genome indexe de la souris mm10.

Le resultat de bowtie2 est des fichiers au format .bam.

[doc bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)


## Etape 4 : Verification de la qualité de l'alignement

samtools flagstat prend en entree des fichiers .bam issus de l'alignement.
Samtools flagstat permet d'acceder a differentes informations sur le mapping de nos echantillons (on recherche des information sur les flag dans un fichier bam).

On peut ainsi acceder a ces informations :
- total : nombre total de reads (somme des reads forward, reverse, single et supplumentaires).
- supplementary : nombre de fois ou un read chimerique a mappu a plusieurs endroits.
- mapped : nombre de reads qui ont mappe sur les contigs et pourcentage mappes sur le total.
- paired in sequence : somme des reads 1 et 2
- properly paired : nombre de reads forward et reverse qui ont mappes sur les même contigs orientes l’un vers l’autre.
- with itself and mate mapped : nombre de reads apparies où les deux sont mappes.
- singletons : nombre de reads apparies ou l’un est mappe et pas l’autre.

[doc samtools_flagstat](http://www.htslib.org/doc/samtools-flagstat.html)


## Etape 5 : Tri des fichiers bam:

samtools sort prend en entree des fichiers .bam.

samtools sort permet de trier dans un fichier .bam les reads en fonction de leur coordonnees.
Ainsi, si deux reads ont les meme coordonnees ils seront ecrits a la suite.
L'utilisation de certains outils, comme Picard, nécessite que les reads dans les fichiers .bam soient tries par leur coordonnees.
On obtient donc en sortie des fichiers tries au format .bam.

[doc samtools_sort](http://www.htslib.org/doc/samtools-sort.html)



## Etape 6 : Marquage et élimination des reads dupliqués:

Picard Markduplicates prend en entree des fichiers .bam tries (les reads sont tries en fonction de leurs coordonnees).

Picard Markduplicates va reperer les duplicats dans les fichiers bam et va les marquer.
Grace a l'option "REMOVE_DUPLICATES=true", les reads dupliqués ne seront pas ecrits dans les fichiers .bam de sortie, ils sont directement éliminés.
Picard fournit egalement un fichier "metrics" qui indique le nombre de reads dupliques.

[doc Picard markDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037225972-MarkDuplicates-Picard-)
[gitHub Picard](https://broadinstitute.github.io/picard/)


## Etape 7 : Vérification du fichier de sortie de Picard

On vérifie à l'aide de samtools flagstat que les reads dupliqués n'ont pas été écrits dans le bam de sortie.
On veut uniquement des reads qui sont uniques et "bien mappés" (flag = 99, 147, 83, 163).

[doc samtools_flagstat](http://www.htslib.org/doc/samtools-flagstat.html)


## Etape 8 : Création d'un bam indexé pour la visualisation sous IGV des reads mappés dédupliqués.

samtools index prend en entrée des fichiers .bam où les reads ont été triés et dédupliqués.

samtools index permet d'indexer des fichiers .bam et de fournir des fichiers .bam.bai afin de visualiser les reads alignes, tries et dedupliques sous IGV.

[doc samtools_index](http://www.htslib.org/doc/samtools-index.html)


## Etape 9 : Peak-calling

Les reads nettoyés et dédupliqués au format .bam sont utilisés en input par MACS2 pour faire le peak-calling.

Le peak-calling est une méthode qui permet d'identifier les zones du génome qui ont été enrichies avec les reads alignés suite à une expérience de ChIP-seq. Pour faire cela, on utilise l'outil MACS2.

J'ai utilisé les options suivantes :
- NAME_peaks.xls : a table with information about called peaks (.xlsx) = une table d'infos sur les peaks.

- NAME_treat_pileup.bdg : pileup signals from treatment sample = empiler les signaux de l'échantillon de traitement (.bedGrapH, option : –bdg or -B)

- NAME_peaks.narrowPeak : contains the peak locations, peak summit, p-value and q-value = contient les emplacements des pics, le sommet des pics, les p-value et q-value (BED 6+4, option : if not set –broad)

- NAME_summits.bed : peak summits locations for every peak = emplacements des sommets pour chaque pic (BED, option : if not set –broad)

Le fichier "peaks.narrowPeak"  est un format BED 6 + 4 (6eres colonnes d'un fichier BED standard (chromosome, coordonnees start, coordonnees end, name, score, brin) avec 4 champs supplementaires (signalValue (Measurement of overall enrichment for the region), pValue, qValue et peak).

[doc macs2 callpeak](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html)
[autre doc macs2 callpeak](https://github.com/macs3-project/MACS)


## Etape 10 : Annotation des peaks

Les peaks sont annotés avec ChIPseeker qui prend en entrée des peaks et des fichiers d'annotation (ex: TxDb.Mmusculus.UCSC.mm10.knownGene ou org.Mm.eg.db). ChIPseeker est un package R (il faut écrire un script R que l'on pourra inclure dans le Snakefile).
On va pouvoir obtenir des listes de peaks avec leur annotation, des graphes pour voir au niveau de quelles régions tombent les peaks annotés et on va pourvoir faire des analyses d'enrichissement fonctionnel.
Les analyses d'enrichissement fonctionnel pourront se faire sur un seul facteur ou en intégrant mes différents facteurs.
On pourra utiliser les bases de données GO (gènes qui ont des choses en communs) ou KEG (réaction et voies de signalisation liés entre eux).

[github de chipseeker](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/12_annotation_functional_analysis.html)
[autre doc de chipseeker](http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#abstract)


# Lancement du snakefile pour test :

## Permet de tester le Snakefile
```{bash, eval=FALSE}

conda activate snakemake
snakemake --configfile config.yaml --snakefile Snakefile -j5 --printshellcmd --use-conda --reason --wrapper-prefix https://raw.githubusercontent.com/tdayris/snakemake-wrappers/ --dryrun

```

## Permet d'afficher une image du pipeline
```{bash, eval=FALSE}
snakemake --configfile config.yaml --snakefile Snakefile -j1 --printshellcmd --use-conda --dryrun --reason --rulegraph | dot -T png > test.png

eog test.png

```


# Modification du Snakefile pour le lancer sur le cluster de l'IGR:

Pour chacune de mes règles du Snakefile, il faut définir la mémoire, le temps d'exécution et le nombre de threads.

## Fonction lambda

En python (existe aussi dans d'autres languages), une fonction lambda est une très petite fonction, qui fait quelque chose d'assez simple et qui est peu utilisée.
Cette fonction lambda n'a pas de nom, mais elle a des arguments.

Ex :
lambda wildcare attempt : 1024 * 10
avec :
lamda : fonction lambda
wildcare : argument 1 (si on veut faire un traitement sur un fichier en particulier par ex)
attempt : argument 2 (si une regle échoue, on peut la relancer en réservant un peu plus de resources)
: = retourne
1024 * 10 : valeur de l'argument attempt



## Le paramètre "resources" dans un Snakefile:
Quand on veut lancer un Snakefile sur le cluster via snakemake, il est essentiel de lui préciser la mémoire et le temps d'exécution dont chaque outil a besoin dans le paramètre "resources" avec mem_mb et time_min.
Pour cela, on utilise des fonctions lambda qui précise la mémoire et le temps d'exécution pour chaque rule. 
Par exemple, pour la règle de Bowtie2:
resources:
        mem_mb = lambda wildcards, attempt: min((1024 * 8) + (2 * attempt), 10 * 1024),
        time_min = lambda wildcards, attempt: 45 + (15 * attempt)

Bowtie2 a besoin de 8/10gb et 45/60min par rapport à la taille de mes fichiers.
min((1024 * 8) + (2 * attempt), 10 * 1024), mémoire nécessaire minimale d'env 8 Gb + ajout de 2 Gb à chaque fois qu'il retente et la mémoire max est 10 Gb.
Pour le temps minimal, on met 45 min + on ajoute 15 min à chaque nouvel essai si le précédent a échoué.


## Le paramètre "threads" dans un Snakefile:
Quand on veut lancer un Snakefile sur le cluster via snakemake, il est également essentiel de lui préciser le nombre de threads pour chaque rule. Ainsi, il ne faut pas hésiter à mettre 20 threads pour Bowtie2 ou fastp. On met un seul thread pour fastQC car il prend les fichiers un par un.
Pour samtools_sort, on met 10 threads car on s'attend à avoir une mémoire de 10 Gb et donc on donne 1 Gb par thread.
Quand l'outil n'est pas multi-threadés (à priori comme Picard et MACS2), il faut quand même préciser threads =1 pour le cluster.



# Lancement d'un Snakefile sur le cluster:

## Lancement d'un Snakefile sur le cluster pour un test (sans exécution):
```{bash, eval=FALSE}
snakemake --profile slurm -n
```

## Lancement d'un Snakefile sur le cluster avec exécution:
```{bash, eval=FALSE}
snakemake --profile slurm

nohup snakemake --profile slurm --rerun-incomplete > nohup.ok.log 2> nohup.errors.log &
nohup snakemake --profile slurm --rerun-incomplete > nohup.complete.log 2>&1 &
```

[option --profile de Snakemake](https://snakemake.readthedocs.io/en/stable/executing/cli.html?highlight=--profile)

l'option --profile slurm de snakemake permet de lancer le Snakefile grace à un fichier de configuration "slurm" dans lequel est écrit tous les environnements et options que l'on veut pour snakemake.
Ce fichier est ici : "~/.config/snakemake/slurm/config.yaml"
Il contient :

restart-times: 3 # snakemake va refaire jusqu'à 3 essais si il y a eu un échec
jobscript: "slurm-jobscript.sh"
cluster: "slurm-submit.py" 
cluster-status: "slurm-status.py"
max-jobs-per-second: 1 # maximum 1 job par seconde (pour qu'on ait le temps de la voir apparaitre à l'écran)
max-status-checks-per-second: 10
local-cores: 1
jobs: 10
keep-going: true
reason: true 
printshellcmds: true
jobname: "{name}.{jobid}.snakejob.sh"
conda-prefix: /home/m_nocente@intra.igr.fr/snakemake/conda # utilise mon conda
wrapper-prefix: https://raw.githubusercontent.com/tdayris/snakemake-wrappers/ # utilise les wrappers de Thibault
use-conda: true # utilise conda
~                 

Il a été crée par Thibault et provient de son GitHub :
[GitHUb Thibault pour slurm](https://github.com/tdayris-perso/slurm)


nohup permet de détacher la commande du parent (on veut que le processus continu comme un orphelin).
& (esperluette) permet de reprendre la main, c'est-à-dire que la commande est mise en arrière plan.
--rerun-incomplete permet de dire à snakemake que l'on veut qu'il recommence le job qu'il a commencé et pas fini (car on l'a arrêté).



## Donner des droits:
Suite à un problème de compte (je n'ai pas le droit de soumettre des jobs sur le cluster pour l'instant), j'ai donné à Thibault (et à tout le monde) des droits sur mon home afin qu'il puisse lancer mon Snakefile pour ne pas perdre de temps.

```{bash, eval=FALSE}
chmod -R 770 $PWD
```

## Protection de mes fichiers output:
Afin d'éviter qu'une personne travaillant sur le cluster puisse modifier mes fichiers output, je les ai protéger en écrivant protected() autour de chaque output.
Par exemple : 
trimmed=protected(["results/fastp/cleaned_filtered_{sample_wildcard}_1.fastq.gz", "results/fastp/cleaned_filtered_{sample_wildcard}_2.fastq.gz"])

[option --protected de Snakemake](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?highlight=protected#protected-and-temporary-files)


## Repertoire de travail sur le cluster

Il faut travailler dans :
/mnt/beegfs/userdata/m_nocente/

Réinstallation de conda, puis de snakemake et je relance mon Snakefile.




# On merge les réplicats des fichiers de peaks avec HOMER mergePeaks


srun /mnt/beegfs/userdata/m_diop/for_homer/bin/mergePeaks -d 100 just_chr_ChIPseq_Chd8_rep1_peaks.narrowPeak just_chr_ChIPseq_Chd8_rep2_peaks.narrowPeak -venn "merge_Chd8_peaks_venn" -prefix merge

## On supprime la ligne commencant par # des fichiers de peaks produits par HOMER mergePeaks
grep "#" -v merge_just_chr_ChIPseq_Chd8_rep1_peaks.narrowPeak_just_chr_ChIPseq_Chd8_rep2_peaks.narrowPeak > merge_modif_Chd8_peaks.narrowPeak

## On garde que les colonnes "chr", "start", "end".
cut -f 2-4 merge_modif_Chd8_peaks.narrowPeak > merge_final_Chd8_peaks.narrowPeak

## Mieux ecrit:
grep "#" -v merge_just_chr_ChIPseq_Chd8_rep1_peaks.narrowPeak_just_chr_ChIPseq_Chd8_rep2_peaks.narrowPeak | cut -f 2-4 > merge_final_Chd8_peaks.narrowPeak



# On compte le nombre de base sous les peaks pour chaque facteur
awk '{SUM += $3-$2} END {print SUM}' merge_final_TBP_peaks.narrowPeak 

Le nombre de bases sous TOUS les peaks pour chaque facteur:
TBP : 73495 pb
Pol2 : 942432 pb
Oct4 : 3469030 pb
Chd8 : 1423805 pb
CTCF : 19756243 pb



# Chercher les peaks communs à CTCF et Oct4:
srun /mnt/beegfs/userdata/m_diop/for_homer/bin/mergePeaks -d 100 merge_just_chr_ChIPseq_Oct4_rep1_peaks.narrowPeak_just_chr_ChIPseq_Oct4_rep2_peaks.narrowPeak just_chr_ChIPseq_CTCF_peaks.narrowPeak -venn "merge_Oct4_CTCF_peaks_venn" -prefix merge

On obtient 3 listes :
- une liste contenant que les peaks communs à CTCF et aux 2 réplicats d'Oct4 : merge_merge_just_chr_ChIPseq_Oct4_rep1_peaks.narrowPeak_just_chr_ChIPseq_Oct4_rep2_peaks.narrowPeak_just_chr_ChIPseq_CTCF_peaks.narrowPeak
- une liste ne contenant que les peaks présents dans CTCF : merge_just_chr_ChIPseq_CTCF_peaks.narrowPeak
- une liste ne contenant que les peaks présents dans Oct4 : merge_merge_just_chr_ChIPseq_Oct4_rep1_peaks.narrowPeak_just_chr_ChIPseq_Oct4_rep2_peaks.narrowPeak



# Récupérer le fichier GTF de mm10, garder uniquement les "start_codon", puis les colonnes chr, start, end puis soustraire -2000 à la colonne start et ajouter + 2000 à la colonne end.

grep "start_codon" hgTables.gtf > only_start_codon_hgTables.gtf

cut -f 1,4,5 only_start_codon_hgTables.gtf > good_columns_only_start_codon_hgTables.gtf

awk '{print $1"\t"$2-2000"\t"$3+2000}' good_columns_only_start_codon_hgTables.gtf > coordonnees_good_columns_only_start_codon_hgTables.bed

Je vérifie sur IGV en important mon fichier .bed final et le gtf du début (comparable à ce qu'on a dans RefSeq avec des transcrits en plus).


# Annotation avec HOMER et stats

srun --mem=10G /mnt/beegfs/userdata/m_diop/for_homer/bin/annotatePeaks.pl merge_just_chr_ChIPseq_TBP_rep1_peaks.narrowPeak_just_chr_ChIPseq_TBP_rep2_peaks.narrowPeak mm10 -annStats stats_TBP_annot_HOMER > annotated_peaks_HOMER_TBP_merge.txt

srun --mem=10G /mnt/beegfs/userdata/m_diop/for_homer/bin/annotatePeaks.pl merge_just_chr_ChIPseq_Chd8_rep1_peaks.narrowPeak_just_chr_ChIPseq_Chd8_rep2_peaks.narrowPeak mm10 -annStats stats_Chd8_annot_HOMER > ./annotated_peaks_HOMER/annotated_peaks_HOMER_Chd8_merge.txt



Dans le fichier de statistiques, on a accès à différentes stats pour les différents elements annotes: 
- log2Ratio : enrichissement (si + enrichissement positif, si - enrichissement negatif)
- logP : log des p-values (plus c'est négatif et plus c'est significatif)


# Détection de motifs 
srun --mem=10G --cpus-per-task=6 /mnt/beegfs/userdata/m_diop/for_homer/bin/findMotifsGenome.pl merge_just_chr_ChIPseq_TBP_rep1_peaks.narrowPeak_just_chr_ChIPseq_TBP_rep2_peaks.narrowPeak mm10 ./motifs -size given -preparsedDir /mnt/beegfs/userdata/m_nocente/Stage_projet_Marina/Projet_ChIP_Marina/results/macs2/motif/preparsed/ -p 6
