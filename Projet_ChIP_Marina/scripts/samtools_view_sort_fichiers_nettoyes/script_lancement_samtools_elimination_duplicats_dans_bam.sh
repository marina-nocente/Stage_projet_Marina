#!/bin/bash

set -euop pipefail
# si il y a une erreur, ou qu'une variable n'est pas definie, ou que le chemin n'est pas bien defini, ou si un pipe a un probleme il va s'arreter tout de suite

# Appelle de ma fonction pour les parametres du cluster et lance un script
source ./scripts/fonction_parler_cluster.sh

## Repertoire qui contiendra les resultats de fastQC:
outputDir="${PWD}/results/samtools_view_sort_fichiers_nettoyes"

## Cluster utilise:
cluster="Slurm"



## Pour mon fichier dans ma liste de fichiers: faire:
for bamFile in ${PWD}/results/picard/Marked_duplicate_*; do

  # Afficher le fichier traite:
  echo "${bamFile}";

  if [ "${cluster}" == "Torque" ]
  then
    # chemin vers mon script
    parler_cluster_Torq ${PWD}/scripts/samtools_view_sort_fichiers_nettoyes/samtools_elimination_duplicat_dans_bam_conda_marina.qsub \
    "bamFile='${bamFile}',outputDir='${outputDir}'" \ # mes variables
    "${outputDir}" # mon repertoire de sortie pour la sortie standard et erreur


  elif [ "${cluster}" == "Slurm" ]
  then
    # chemin vers mon script
    parler_cluster_Slurm ${PWD}/scripts/samtools_view_sort_fichiers_nettoyes/samtools_elimination_duplicat_dans_bam_conda_marina.qsub \
    "bamFile='${bamFile}',outputDir='${outputDir}'" # mes variables


  else [ "${cluster}" == "" ]
    echo "Preciser le nom du cluster utilise"

  fi

done
