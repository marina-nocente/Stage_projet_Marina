Contenu:
=========
Ceci est le resultat du mapping de bowtie2 sur {{snakemake.input}}

Fichiers input:
===============
Les reads nettoyes par fastp sont mappes sur le genome indexe (par bowtie2) de la souris (version mm10).
Bowtie utilise donc entree des fichiers .fastq.gz et un fichier fasta indexe.

Outil utilise et resultats:
===========================
D'apres sa documentation, Bowtie 2 est un outil ultrarapide et efficace en mémoire pour aligner les reads de séquençage sur de longues séquences de référence (comme les genomes de mammifere).
Bowtie 2 peut faire des alignements avec des gaps, des alignements locaux et paired-end.
Par défaut, Bowtie 2 effectue un alignement de read "de bout en bout" (= end-to-end). C'est-à-dire qu'il recherche les alignements impliquant toutes les bases du read.
Lorsque l'option –local est spécifiée, Bowtie 2 effectue un alignement local du read. Dans ce mode, Bowtie 2 peut "découper" ou "couper" certaines bases du read à l'une ou aux deux extrémités de l'alignement si cela maximise le score d'alignement.

Bowtie2 a tourne avec ses parametres par defaut (alignement end-to-end) et j'ai precise le parametre --align-paired-reads pour lui indiquer que les reads sont paires.
Bowtie2 mappe les reads paires sur le genome indexe de la souris mm10.
J'ai pris le FASTA du genome de la souris sur GENCODE: "contenu Genome sequence, primary assembly (GRCm38) mm10, region PRI. Nucleotide sequence of the GRCm38 primary genome assemby (chromosomes and scaffolds (suffisant pour ce que l'on veut nous, pas besoin du reste)). The sequence region names are in the same as in the GTF/GFF file".

Le resultat de bowtie2 est des fichiers au format .bam.
