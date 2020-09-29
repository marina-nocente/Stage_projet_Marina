Contenu:
=========
Cette partie est le resultat du controle de fastQC sur {{snakemake.input}}.

Fichiers input:
===============
Les fichiers input sont les fichiers fastq.gz fournis par la plateforme de sequencage.

Outil utilise et resultats:
===========================
L'outil fastQC permet de visualiser la qualite et le contenu des reads apres sequencage.
Il fourni pour chaque read (R1 et R2) un fichier .zip et un fichier .html.

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
