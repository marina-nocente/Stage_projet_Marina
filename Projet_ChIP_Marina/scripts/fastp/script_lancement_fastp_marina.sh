#!/bin/bash

set -euop pipefail
# si il y a une erreur, ou qu'une variable n'est pas definie, ou que le chemin n'est pas bien defini, ou si un pipe a un probleme il va s'arreter tout de suite

#PBS -q remod
#PBS -l ncpus=8

## Explication des parametres pour le cluster (PBS)
# -q : notre noeud sur lequel on veut travailler
# -l ncpus : nombre de cpus

source fonction_parler_cluster

## Definition du repertoire de sortie pour les resultats de fastp:
outputDir="/store/EQUIPES/REMOD/200120_novogene/Marina/RNA-seq/fastp/"

## Definition du cluster utilise
cluster="Slurm"

## Pour faire tourner le script sur tous les fichiers reads 1:
for fastqFile1 in /store/EQUIPES/REMOD/200120_novogene/raw_data/C57_*1.fastq.gz ; do

	# Afficher le fichier read1 utilise
	echo "${fastqFile1}";

	# Afficher le nom du fichier en enlevant ce qui est avant le dernier /
	echo "${fastqFile1##*/}";

	# Creation d'une variable fastqFile2 qui contient les reads2 (on a remplace dans fastqFile1 "1.fastq.gz" par "2.fastq.gz")
	fastqFile2="${fastqFile1/1.fastq.gz/2.fastq.gz}";

	# Afficher le nom du fichier en enlevant ce qui est avant le dernier /
	echo "${fastqFile2##*/}";

  if [ "${cluster}" == "Torq" ]
  then
    # On lance le script contena le programme et les parametres fastp:
	  qsub /store/EQUIPES/REMOD/200120_novogene/Marina/fastp_conda_marina.qsub \ # chemin vers mon script
	  -o "${outputDir}" \ # le repertoire de sortie ou fastp va ecrire les nouveaux fichiers
	  -e "${outputDir}" \ # le chemin pour le flux d'erreur standard
	  -v "fastqFile1='${fastqFile1}',fastqFile2='${fastqFile2}',outputDir='${outputDir}'"

  elif [ ${cluster} == "Slurm" ]
  then
    sbatch /store/EQUIPES/REMOD/200120_novogene/Marina/fastp_conda_marina.qsub \ # chemin vers mon script
	  -o "${outputDir}" \ # le repertoire de sortie ou fastp va ecrire les nouveaux fichiers
	  -e "${outputDir}" \ # le chemin pour le flux erreur standard
	  --export "fastqFile1='${fastqFile1}',fastqFile2='${fastqFile2}',outputDir='${outputDir}'"

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
