Contenu:
=========
Ceci est le resultat du controle de fastp sur {{snakemake.input}}


Fichiers input:
===============
Les fichiers input sont les fichiers fastq.gz fournis par la plateforme de sequencage et dont la qualite a ete evaluee par fastQC.


Outil utilise et resultats:
===========================
L'outil fastp permet de nettoyer les reads des adaptateurs et si besoin de couper une partie des reads de mauvaise qualite.

Voici quelques caracteristiques de l'outil fastp (d'apres sa documentation):
-	suit la qualite des donnees avant et apres trimming (courbes de qualité, contenu de base, KMER, Q20 / Q30, rapport GC, duplication, contenu en adaptateur ...).
-	filtre les mauvaises reads (qualité trop faible, trop courte ou trop de N ...)
-	coupe des bases de faible qualité en 5 'et 3' des reads en évaluant la qualité moyenne d'une fenêtre coulissante.
-	elimine automatiquement les adaptateurs.
-	fournit un rapport des resultats au format JSON ou HTML (comme fastQC).

On obtient ainsi des reads nettoyes au format fastq.gz et un rapport de ce que fastp a fait au format .html ou .json.
