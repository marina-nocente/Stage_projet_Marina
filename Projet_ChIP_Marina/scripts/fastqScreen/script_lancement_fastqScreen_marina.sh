#!/bin/bash

set -euop pipefail
# si il y a une erreur, ou qu'une variable n'est pas definie, ou que le chemin n'est pas bien defini, ou si un pipe a un probleme il va s'arreter tout de suite

### Explication des parametres pour le cluster (PBS)
## notre noeud
#PBS -q remod

## nombre de cpu par job
#PBS -l ncpus=8

### Definition de variables:
outputDir="/store/EQUIPES/REMOD/200120_novogene/Marina/RNA-seq/fastqScreen"
fastqScreenConf="/store/EQUIPES/REMOD/scripts/fastq_screen.conf"
cluster="Slurm"

## Pour mon fichier dans ma liste de fichiers: faire:
for fastqFile in /store/EQUIPES/REMOD/200120_novogene/Marina/RNA-seq/fastq_sans_rRNA/cleaned_filtered_C57_*.fastq; do

  # Afficher le nom du fichier en cours de traitement:
  echo "${fastqFile}";

  if [ "${cluster}" == "Torq" ]
  then
    # Lancer le script fastScreen:
    qsub /store/EQUIPES/REMOD/200120_novogene/Marina/fastqScreen_marina_conda.qsub \
    -o "${outputDir}" \ # le repertoire de sortie ou fastp va ecrire les nouveaux fichiers
    -e  "${outputDir}" \ # le chemin pour le flux d'erreur standard
    -v "fastqFile='${fastqFile}',outputDir='${outputDir}',fastqScreenConf='${fastqScreenConf}'"

  elif [ "${cluster}" == "Slurm" ]
  then
    sbatch /store/EQUIPES/REMOD/200120_novogene/Marina/fastqScreen_marina_conda.qsub \
    -o "${outputDir}" \ # le repertoire de sortie ou fastp va ecrire les nouveaux fichiers
    -e  "${outputDir}" \ # le chemin pour le flux d erreur standard
    --export "fastqFile='${fastqFile}',outputDir='${outputDir}',fastqScreenConf='${fastqScreenConf}'"

  else [ "${cluster}" == "" ]
    echo "Preciser le nom du cluster utilise"

  fi

done

## Parametres de qsub (pour lancer une ligne de commande sur le cluster I2BC):
# mettre le chemin du script de fastqScreen
# -v : variable_list
# outputDir : variable qui contient le chemin ou les nouveaux fichiers seront ecrits
# fastqScreenConf : renvoie au fichier de configuration de fastqScreen
