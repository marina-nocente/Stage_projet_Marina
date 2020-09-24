Contenu:
=========
Ceci est le resultat du peak-calling de MACS2 sur {{snakemake.input}}

Fichiers input:
===============
Les reads nettoyes et dedepliques au format .bam sont utilises en input par MACS2 pour faire le peak-calling.

Outil utilise et resultats:
===========================
Le peak-calling est une methode qui permet d'identifier les zones du genome qui ont ete enrichiees avec les reads alignes suite a une experience de ChIP-seq.
Pour faire cela, on utilise l'outil MACS2.

J'ai utilise les options suivantes :
- NAME_peaks.xls : a table with information about called peaks (.xlsx) = une table d'infos sur les peaks.

- NAME_treat_pileup.bdg : pileup signals from treatment sample = empiler les signaux de l'échantillon de traitement (.bedGrapH, option : –bdg or -B)

- NAME_peaks.narrowPeak : contains the peak locations, peak summit, p-value and q-value = contient les emplacements des pics, le sommet des pics, les p-value et q-value (BED 6+4, option : if not set –broad)

- NAME_summits.bed : peak summits locations for every peak = emplacements des sommets pour chaque pic (BED, option : if not set –broad)

Le fichier "peaks.narrowPeak"  est un format BED 6 + 4 (6eres colonnes d'un fichier BED standard (chromosome, coordonnees start, coordonnees end, name, score, brin) avec 4 champs supplementaires (signalValue (Measurement of overall enrichment for the region), pValue, qValue et peak)
