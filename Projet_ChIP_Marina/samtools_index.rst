Contenu:
=========
Ceci est le resultat de samtools index sur {{snakemake.input}}


Fichiers input:
===============
samtools index prend en entree des fichiers .bam ou les reads ont ete tries et dedupliques.


Outil utilise et resultats:
===========================
samtools index permet d'indexer des fichiers .bam et de fournir des fichiers .bam.bai afin de visualiser les reads alignes, tries et dedupliques sous IGV.
