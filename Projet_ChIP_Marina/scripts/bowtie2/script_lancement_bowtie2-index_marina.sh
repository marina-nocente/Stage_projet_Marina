#!/bin/bash

set -euop pipefail
# si il y a une erreur, ou qu'une variable n'est pas definie, ou que le chemin n'est pas bien defini, ou si un pipe a un probleme il va s'arreter tout de suite


#PBS -q remod
#PBS -l ncpus=8

## Definition du repertoire de sortie pour les resultats de bowtie2:
outputDir="/store/EQUIPES/REMOD/200120_novogene/Marina/RNA-seq/mapping_bowtie2/"

## Definition du cluster utilise
cluster="Slurm"


## Pour mon fichier read1 : faire:
for fastqFile1 in /store/EQUIPES/REMOD/200120_novogene/Marina/RNA-seq/fastq_sans_rRNA/cleaned_filtered_C57_*paired1.fastq ; do
	echo "${fastqFile1}";
	echo "${fastqFile1##*/}";
	fastqFile2="${fastqFile1/paired1.fastq/paired2.fastq}";
	echo "${fastqFile2##*/}";

  if [ "${cluster}" == "Torq" ]
  then
    # Lancer le script fastQC:
    qsub /store/EQUIPES/REMOD/200120_novogene/Marina/bowtie2_index_conda_marina.qsub -o ${outputDir} -e ${outputDir} -v "fastqFile1='${fastqFile1}'","fastqFile2='${fastqFile2}'","outputDir='${outputDir}'"
	  echo "______________";

  elif [ ${cluster} == "Slurm" ]
  then
    sbatch /store/EQUIPES/REMOD/200120_novogene/Marina/bowtie2_index_conda_marina.qsub -o ${outputDir} -e ${outputDir} --export "fastqFile1='${fastqFile1}'","fastqFile2='${fastqFile2}'","outputDir='${outputDir}'"
	  echo "______________";

  else [ "${cluster}" == "" ]
    echo "Preciser le nom du cluster utilise"

  fi

done
