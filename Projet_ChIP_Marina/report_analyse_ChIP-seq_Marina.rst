fastQC sur mes échantillons apres sequencage:
=============================================
On a utilise fastQC pour verifier la qualite des echantillons apres sequençage.


fastp sur mes echantillons fastq.gz:
====================================
On utilise fastp pour l'etape de trimming : on elimine les adaptateurs des fichiers fastq.gz et on les nettoie.


bowtie2 sur les echantillons fastq.gz nettoyes:
===============================================
On utilise bowtie2 pour mapper les reads nettoyes sur le genome indexe de la souris (mm10)


samtools_flagstat sur mes bam:
==============================
On utilise bowtie2 pour verifier la qualite de l'alignement.


samtools_sort sur mes bam:
==========================
On utilise samtools_sort pour trier les reads en fonction de leurs coordonnees.


Picard markduplicates sur les bam tries:
========================================
On utilise Picard markduplicates pour marquer les reads dupliques et les eliminer. Le nouveau bam trie ne conteient pas de reads dupliques.


samtools_flagstat sur mes bam tries et dedupliques:
===================================================
On verifie que nos bam tries et dedupliques n'ont bien pas de reads dupliques


samtools_index sur mes bam tries et dedupliques:
================================================
On utilise samtools_index pour indexer mes bam tries et dedupliques


MACS2 callpeak bam tries et dedupliques:
========================================
On utilise MACS2 callpeak pour generer des peaks a partir de mes bam.
