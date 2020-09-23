Contenu:
=========
Ceci est le resultat de multiQC sur {{snakemake.input}}


Fichiers input:
===============
multiQC prend en entree des fichiers de fastQC (.zip), de fastp (.fastq), de samtools_flagstat (.bam.flagstat) et de picard (.metrics.txt).


Outil utilise et resultats:
===========================

multiQC permet de faire un synthese des resultats de differents outils et de presenter une visualisation agreable et plus facilement comparable.

Dans mon rapport multiQC, on peut voir les resultats de :
- fastQC
- fastp
- samtools_flagstat apres l'alignement avec bowtie2
- picard markduplicates avec les metriques
- samtools_flagstat apres picard pour verifier que les reads dupliques ont tous bien ete supprimes.
