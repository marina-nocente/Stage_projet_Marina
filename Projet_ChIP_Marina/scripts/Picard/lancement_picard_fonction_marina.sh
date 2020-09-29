#!/bin/bash

set -euop pipefail
# si il y a une erreur, ou qu'une variable n'est pas definie, ou que le chemin n'est pas bien defini, ou si un pipe a un probleme il va s'arreter tout de suite

# Appelle de ma fonction pour les parametres du cluster et lance un script
source /Users/mn242062/Desktop/Stage_projet_Marina/Projet_ChIP_Marina/scripts/fonction_parler_cluster.sh

### Definition de variables:
outputDir="${PWD}/results/picard"

cluster="Slurm"

## Pour mon fichier dans ma liste de fichiers: faire:
for bamFile in ${PWD}/results/samtools_sort_tri/*.bam ; do
	echo "${bamFile}";


	if [ "${cluster}" == "Torque" ]
	then
		# chemin vers mon script
		parler_cluster_Torq ${PWD}/scripts/Picard/picard_conda_marina.qsub \
		"bamFile='${bamFile}',outputDir='${outputDir}'" \ # mes variables
		"${outputDir}" # mon repertoire de sortie pour la sortie standard et erreur


	elif [ "${cluster}" == "Slurm" ]
	then
		# chemin vers mon script
		parler_cluster_Slurm ${PWD}/scripts/Picard/picard_conda_marina.qsub \
		"bamFile='${bamFile}',outputDir='${outputDir}'" # mes variables


	else [ "${cluster}" == "" ]
		echo "Preciser le nom du cluster utilise"

	fi

done
