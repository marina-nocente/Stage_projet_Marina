Contenu:
=========
Ceci est le resultat de samtools sort sur {{snakemake.input}}


Fichiers input:
===============
samtools sort prend en entree des fichiers .bam.


Outil utilise et resultats:
===========================
samtools sort permet de trier dans un fichier .bam les reads en fonction de leur coordonnees.
Ainsi, si deux reads ont les meme coordonnees ils seront ecrits a la suite.
L'utilisation de certains outils, comme Picard, nécessite que les reads dans les fichiers .bam soient tries par leur coordonnees.
On obtient donc des fichiers tries au format .bam.
