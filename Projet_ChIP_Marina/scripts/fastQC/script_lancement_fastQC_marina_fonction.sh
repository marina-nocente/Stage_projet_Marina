#!/bin/bash

set -euop pipefail
# si il y a une erreur, ou qu'une variable n'est pas definie, ou que le chemin n'est pas bien defini, ou si un pipe a un probleme il va s'arreter tout de suite

# Appelle de ma fonction pour les parametres du cluster et lance un script
source /Users/mn242062/Desktop/Stage_projet_Marina/Projet_ChIP_Marina/scripts/fonction_parler_cluster.sh

## Repertoire qui contiendra les resultats de fastQC:
outputDir="/home/mnocente/Bureau/Projet_ChIP_Marina/scripts/samtools_flagstat/results"

## Cluster utilise:
cluster="Slurm"

## Chemin de mes donnees:
# data=(/store/EQUIPES/REMOD/200120_novogene/raw_data/C57_*.fastq.gz)


## Pour mon fichier dans ma liste de fichiers: faire:
for fastqFile in /Users/mn242062/Desktop/Stage_projet_Marina/Projet_ChIP_Marina/fichiers_test/*.fastq.gz; do

  # Afficher le fichier traite:
  echo "${fastqFile}";

  if [ "${cluster}" == "Torque" ]
  then
    # chemin vers mon script
    parler_cluster_Torq /home/mnocente/Bureau/Projet_ChIP_Marina/scripts/fastQC/fastQC_conda_marina.qsub \
    "fastqFile='${fastqFile}',outputDir='${outputDir}'" \ # mes variables
    "${outputDir}" # mon repertoire de sortie pour la sortie standard et erreur


  elif [ "${cluster}" == "Slurm" ]
  then
    # chemin vers mon script
    parler_cluster_Slurm /home/mnocente/Bureau/Projet_ChIP_Marina/scripts/fastQC/fastQC_conda_marina.qsub \
    "fastqFile='${fastqFile}',outputDir='${outputDir}'" # mes variables


  else [ "${cluster}" == "" ]
    echo "Preciser le nom du cluster utilise"

  fi

done
