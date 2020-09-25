#!/bin/bash

set -euop pipefail
# si il y a une erreur, ou qu'une variable n'est pas definie, ou que le chemin n'est pas bien defini, ou si un pipe a un probleme il va s'arreter tout de suite

# Appelle de ma fonction pour les parametres du cluster et lance un script
source ./scripts/fonction_parler_cluster.sh # il faut etre dans le repertoire "Projet_ChIP_Marina" pour que cela fonctionne

## Repertoire qui contiendra les resultats de fastQC:
outputDir="${PWD}/results/fastQC"

## Cluster utilise:
cluster=""

## Chemin de mes donnees:
# data=(/store/EQUIPES/REMOD/200120_novogene/raw_data/C57_*.fastq.gz)


## Pour mon fichier dans ma liste de fichiers: faire:
for fastqFile in ${PWD}/raw_data/*.fastq.gz; do

  # Afficher le fichier traite:
  echo "${fastqFile}";

  if [ "${cluster}" == "Torque" ]
  then
    # chemin vers mon script
    parler_cluster_Torq ${PWD}/scripts/fastQC/fastQC_conda_marina.qsub \
    "fastqFile='${fastqFile}',outputDir='${outputDir}'" \ # mes variables
    "${outputDir}" # mon repertoire de sortie pour la sortie standard et erreur


  elif [ "${cluster}" == "Slurm" ]
  then
    # chemin vers mon script
    parler_cluster_Slurm ${PWD}/scripts/fastQC/fastQC_conda_marina.qsub \
    "fastqFile='${fastqFile}',outputDir='${outputDir}'" # mes variables


  else [ "${cluster}" == "" ]
    echo "Preciser le nom du cluster utilise"
    bash ${PWD}/scripts/fastQC/fastQC_conda_marina.qsub \
    "fastqFile='${fastqFile}',outputDir='${outputDir}'" # mes variables
  fi

done
