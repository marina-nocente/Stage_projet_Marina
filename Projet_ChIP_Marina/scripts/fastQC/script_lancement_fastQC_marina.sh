#!/bin/bash

set -euop pipefail
# si il y a une erreur, ou qu'une variable n'est pas definie, ou que le chemin n'est pas bien defini, ou si un pipe a un probleme il va s'arreter tout de suite

##PBS -q remod


## Repertoire qui contiendra les resultats de fastQC:
outputDir="/store/EQUIPES/REMOD/200120_novogene/Marina/RNA-seq/Fastqc"

## Cluster utilise:
cluster="Slurm"

## Chemin de mes donnees:
# data=(/store/EQUIPES/REMOD/200120_novogene/raw_data/C57_*.fastq.gz)

## Pour mon fichier dans ma liste de fichiers: faire:
for fastqFile in /store/EQUIPES/REMOD/200120_novogene/raw_data/C57_*.fastq.gz; do

  # Afficher le fichier traite:
  echo "${fastqFile}";

  if [ "${cluster}" == "Torq" ]
  then
    # Lancer le script fastQC:
    qsub /store/EQUIPES/REMOD/200120_novogene/Marina/fastQC_conda_marina.qsub -v fastqFile="${fastqFile}",outputDir="${outputDir}"

  elif [ "${cluster}" == "Slurm" ]
  then
    sbatch /store/EQUIPES/REMOD/200120_novogene/Marina/fastQC_conda_marina.qsub --export fastqFile="${fastqFile}",outputDir="${outputDir}"

  else [ "${cluster}" == "" ]
    echo "Preciser le nom du cluster utilise"

  fi

done


## Parametres de qsub (pour lancer une ligne de commande sur le cluster I2BC):
# mettre le chemin du script de fastp
# -v : variable_list
# --export : donne toutes les variables utilisees dans SLURM.
# outputDir : variable qui contient le chemin ou les nouveaux fichiers seront ecrits
