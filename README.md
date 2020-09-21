# Stage_projet_Marina

Ce projet GitHub contient un répertoire "Projet_ChIP_Marina" qui contient un répertoire "scripts", un répertoire "fichiers_test" et un fichier "Snakefile" et son fichier de configuration "config.yaml".


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


Le projet GitHub contient également le fichier "feuille_route.ods" qui répertorie toute les informations connues de mes échantillons, les différentes étapes de l'analyse à faire et leur avancement.

Un répertoire "raw_data" contenant les fichiers fastq.gz des 9 échantillons a été créé.
La feuille de route a été mise à jour avec le nom des échantillons de la plateforme de séquencage et le nom avec lequel je les ai renommés.
