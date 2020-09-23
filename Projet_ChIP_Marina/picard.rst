Contenu:
=========
Ceci est le resultat de picard Markduplicates sur {{snakemake.input}}


Fichiers input:
===============
picard Markduplicates prend en entree des fichiers .bam tries (les reads sont tries en fonction de leurs coordonnees).


Outil utilise et resultats:
===========================
Picard Markduplicates va reperer les duplicats dans les fichiers bam et va les marquer.
Grace a l'option "REMOVE_DUPLICATES=true", les reads dupliques ne seront pas ecrits dans les fichiers .bam de sortie.
Picard fournit egalement un fichier "metrics" qui indique le nombre de reads dupliques.

__https://gatk.broadinstitute.org/hc/en-us/articles/360037225972-MarkDuplicates-Picard-
