#!/bin/bash

set -euop pipefail
# si il y a une erreur, ou qu'une variable n'est pas definie, ou que le chemin n'est pas bien defini, ou si un pipe a un probleme il va s'arreter tout de suite


source /home/mnocente/Bureau/fonction_parler_cluster.sh

## Definition du repertoire de sortie pour les resultats de fastp:
outputDir="/home/mnocente/Bureau/Projet_ChIP_Marina/scripts/fastQC/results"

## Definition du cluster utilise
cluster="Slurm"

## Pour faire tourner le script sur tous les fichiers reads 1:
for fastqFile1 in /home/mnocente/Bureau/Projet_ChIP_Marina/scripts/fastQC/test/*.fastq.gz ; do

	# Afficher le fichier read1 utilise
	echo "${fastqFile1}";

	# Afficher le nom du fichier en enlevant ce qui est avant le dernier /
	echo "${fastqFile1##*/}";

	# Creation d'une variable fastqFile2 qui contient les reads2 (on a remplace dans fastqFile1 "1.fastq.gz" par "2.fastq.gz")
	fastqFile2="${fastqFile1/1.fastq.gz/2.fastq.gz}";

	# Afficher le nom du fichier en enlevant ce qui est avant le dernier /
	echo "${fastqFile2##*/}";

# Creation de la variable fastqFile:
fastqFile="${fastqFile1/_R1.fastq.gz/}" # on garde juste le nom de "base" du fichier (on enleve read1 et l'extension)


  if [ "${cluster}" == "Torq" ]
  then
    # chemin vers mon script
    parler_cluster_Torq /home/mnocente/Bureau/Projet_ChIP_Marina/scripts/fastp/fastp_conda_marina.qsub \
    "fastqFile='${fastqFile}',fastqFile1='${fastqFile1}',fastqFile2='${fastqFile2}',outputDir='${outputDir}'" \ # mes variables
    "${outputDir}" # mon repertoire de sortie pour la sortie standard et erreur


  elif [ "${cluster}" == "Slurm" ]
  then
    # chemin vers mon script
    parler_cluster_Slurm /home/mnocente/Bureau/Projet_ChIP_Marina/scripts/fastp/fastp_conda_marina.qsub \
    "fastqFile='${fastqFile}',fastqFile1='${fastqFile1}',fastqFile2='${fastqFile2}',outputDir='${outputDir}'" # mes variables

  else [ "${cluster}" == "" ]
    echo "Preciser le nom du cluster utilise"

  fi

# -v : elargie la liste des variable qui sont exportees vers la tache.
# variable_list est une liste de chaines de caracteres separees par des virgules.
# Chaque élément est de la forme : variable ou variable=value. Ces variables et leurs valeurs seront passes a la tache.
# http://www.delafond.org/traducmanfr/man/man1/qsub.1.html

	# Afficher une ligne de separation entre les fichiers
	echo "______________"

done
