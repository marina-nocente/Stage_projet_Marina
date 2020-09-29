#!/bin/bash

## Je veux concatener dans le meme fichier tous les fichiers peaks.narrowPeak (provenant des ChIP Oct4, TBP, CTCF, Pol2 et Chd8), puis je veux trier les peaks par coordonnees, puis je veux garder uniquement les coordonnees uniques.


## les fichiers peaks.narrowPeak sont des BED 6 + 4 (6eres colonnes d'un fichier BED standard (chromosome, coordonnees start, coordonnees end, name, score, brin) avec 4 champs supplementaires (signalValue (Measurement of overall enrichment for the region), pValue, qValue et peak).


cat ./results/macs2/*_peaks.narrowPeak | sort | unique | cut -c1-3 > "./results/bed_peaks_all_tries/tri_unique_all_peaks.narrowPeak"



## On s'interesse ensuite a l'outil bedtools pour pouvoir croiser les coordonnees de nos peaks et l'endroit ou le peak est present dans le ChIP

bedtools intersect -a "./results/bed_peaks_all_tries/tri_unique_all_peaks.narrowPeak" \
                  -b
