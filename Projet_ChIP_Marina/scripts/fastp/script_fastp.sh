#!/bin/bash

set -euop pipefail
# si il y a une erreur, ou qu'une variable n'est pas definie, ou que le chemin n'est pas bien defini, ou si un pipe a un probleme il va s'arreter tout de suite


source ./scripts/fonction_parler_cluster.sh

# J'ai une liste de 2 échantillons.
liste="1 2"

#Pour chaque element de ma liste, j'ajoute R1 puis R2 dans le nom.
for  nom_ech in $liste
do
  echo "nom_ech =" $nom_ech

R1="${nom_ech}_R1.fastq.gz"
echo $R1

R2="${nom_ech}_R2.fastq.gz"
echo $R2

prefix="${PWD}/raw_data"
#echo $prefix

#J'obtiens alors 2 paires de fichiers avec leur chemin
R1_path="${prefix}/$R1"
#echo $R1_path

R2_path="${prefix}/$R2"
#echo $R2_path


# Je verifie que mes fichiers existent:

if [ -f  ${R1_path} ]
then
    echo "Mes fichiers existent"
    parler_cluster_Slurm ${PWD}/scripts/fastp/fastp_conda_marina.qsub "fastqFileR1='${R1_path},fastqFileR2='${R2_path}'"

else
    echo "Mes fichiers n'existent pas"
fi

if [ -f  ${R2_path} ]
then
    echo "Mes fichiers existent"
    #parler_cluster_Slurm ${PWD}/scripts/fastp/fastp_conda_marina.qsub "fastqFileR1='${R1_path},fastqFileR2='${R2_path}'"
else
    echo "Mes fichiers n'existent pas"
fi


done


#Je peux ensuite lancer mon fastp.
