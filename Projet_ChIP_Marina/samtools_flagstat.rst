Contenu:
=========
Ceci est le resultat de samtools_flagstat apres le mapping et apres picard sur {{snakemake.input}}


Fichiers input:
===============
samtools flagstat prend en entree des fichiers .bam.


Outil utilise et resultats:
===========================

Samtools flagstat permet d'acceder a differentes informations sur le mapping de nos echantillons (on recherche des information sur les flag dans un fichier bam).

On peut ainsi acceder a ces informations :
- total : nombre total de reads (somme des reads forward, reverse, single et supplumentaires).
- supplementary : nombre de fois ou un read chimerique a mappu a plusieurs endroits.
- mapped : nombre de reads qui ont mappe sur les contigs et pourcentage mappes sur le total.
- paired in sequence : somme des reads 1 et 2
- properly paired : nombre de reads forward et reverse qui ont mappes sur les même contigs orientes l’un vers l’autre.
- with itself and mate mapped : nombre de reads apparies où les deux sont mappes.
- singletons : nombre de reads apparies ou l’un est mappe et pas l’autre.


J'utilise samtools flagstat :
- une premiere fois pour verifier la qualite de l'alignement et
- une seconde fois pour verifier que les duplicats ont bien ete supprimes apres Picard markduplicates.
