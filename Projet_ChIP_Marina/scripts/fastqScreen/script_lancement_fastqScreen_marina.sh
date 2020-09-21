#!/bin/bash

set -euop pipefail
# si il y a une erreur, ou qu'une variable n'est pas definie, ou que le chemin n'est pas bien defini, ou si un pipe a un probleme il va s'arreter tout de suite

# Appelle de ma fonction pour les parametres du cluster et lance un script
source ./scripts/fonction_parler_cluster.sh

### Definition de variables:
outputDir="${PWD}/results/fastqScreen"
fastqScreenConf="/store/EQUIPES/REMOD/scripts/fastq_screen.conf" ### a modifier quand je serais sur le cluster ###
cluster="Slurm"

## Pour mon fichier dans ma liste de fichiers: faire:
for fastqFile in ${PWD}/raw_data/*.fastq.gz; do

  # Afficher le nom du fichier en cours de traitement:
  echo "${fastqFile}";

  if [ "${cluster}" == "Torque" ]
  then
    # chemin vers mon script
    parler_cluster_Torq ${PWD}/scripts/fastScreen/fastqScreen_marina_conda.qsub \
    "fastqFile='${fastqFile}',outputDir='${outputDir}'" \ # mes variables
    "${outputDir}" # mon repertoire de sortie pour la sortie standard et erreur


  elif [ "${cluster}" == "Slurm" ]
  then
    # chemin vers mon script
    parler_cluster_Slurm ${PWD}/scripts/fastScreen/fastqScreen_marina_conda.qsub \
    "fastqFile='${fastqFile}',outputDir='${outputDir}'" # mes variables


  else [ "${cluster}" == "" ]
    echo "Preciser le nom du cluster utilise"

  fi

done

## Parametres de qsub (pour lancer une ligne de commande sur le cluster I2BC):
# mettre le chemin du script de fastqScreen
# -v : variable_list
# outputDir : variable qui contient le chemin ou les nouveaux fichiers seront ecrits
# fastqScreenConf : renvoie au fichier de configuration de fastqScreen
