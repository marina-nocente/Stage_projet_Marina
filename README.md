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


# Explication du Snakefile:

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


## Etape 5 : Tri des fichiers bam:

samtools sort prend en entree des fichiers .bam.

samtools sort permet de trier dans un fichier .bam les reads en fonction de leur coordonnees.
Ainsi, si deux reads ont les meme coordonnees ils seront ecrits a la suite.
L'utilisation de certains outils, comme Picard, nécessite que les reads dans les fichiers .bam soient tries par leur coordonnees.
On obtient donc en sortie des fichiers tries au format .bam.


## Etape 6 : Marquage et élimination des reads dupliqués:

Picard Markduplicates prend en entree des fichiers .bam tries (les reads sont tries en fonction de leurs coordonnees).

Picard Markduplicates va reperer les duplicats dans les fichiers bam et va les marquer.
Grace a l'option "REMOVE_DUPLICATES=true", les reads dupliqués ne seront pas ecrits dans les fichiers .bam de sortie, ils sont directement éliminés.
Picard fournit egalement un fichier "metrics" qui indique le nombre de reads dupliques.

__https://gatk.broadinstitute.org/hc/en-us/articles/360037225972-MarkDuplicates-Picard-


## Etape 7 : Vérification du fichier de sortie de Picard

On vérifie à l'aide de samtools flagstat que les reads dupliqués n'ont pas été écrits dans le bam de sortie.
On veut uniquement des reads qui sont uniques et "bien mappés" (flag = 99, 147, 83, 163).


## Etape 8 : Création d'un bam indexé pour la visualisation sous IGV des reads mappés dédupliqués.

samtools index prend en entrée des fichiers .bam où les reads ont été triés et dédupliqués.

samtools index permet d'indexer des fichiers .bam et de fournir des fichiers .bam.bai afin de visualiser les reads alignes, tries et dedupliques sous IGV.


## Etape 9 : Peak-calling

Les reads nettoyés et dédupliqués au format .bam sont utilisés en input par MACS2 pour faire le peak-calling.

Le peak-calling est une méthode qui permet d'identifier les zones du génome qui ont été enrichies avec les reads alignés suite à une expérience de ChIP-seq. Pour faire cela, on utilise l'outil MACS2.

J'ai utilisé les options suivantes :
- NAME_peaks.xls : a table with information about called peaks (.xlsx) = une table d'infos sur les peaks.

- NAME_treat_pileup.bdg : pileup signals from treatment sample = empiler les signaux de l'échantillon de traitement (.bedGrapH, option : –bdg or -B)

- NAME_peaks.narrowPeak : contains the peak locations, peak summit, p-value and q-value = contient les emplacements des pics, le sommet des pics, les p-value et q-value (BED 6+4, option : if not set –broad)

- NAME_summits.bed : peak summits locations for every peak = emplacements des sommets pour chaque pic (BED, option : if not set –broad)

Le fichier "peaks.narrowPeak"  est un format BED 6 + 4 (6eres colonnes d'un fichier BED standard (chromosome, coordonnees start, coordonnees end, name, score, brin) avec 4 champs supplementaires (signalValue (Measurement of overall enrichment for the region), pValue, qValue et peak).



# Lancement du snakefile pour test :

## Permet de tester le Snakefile
```{bash, eval=FALSE}

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
attempt : argument 2 (si une regle échoue, on peut la relancer en reservant un peu plus de resources)
: = retourne
1024 * 10 : valeur de l'argument attempt



## le paramètre "resource" dans un Snakefile:








