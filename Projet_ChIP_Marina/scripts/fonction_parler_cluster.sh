#!/bin/bash

### Fonction pour traiter un script et des variables toujours de la meme facon sur le cluster ###

function parler_cluster_Torq() {
  local script="${1}"
  local variable="${2}"
  local repertoire_sortie="${3}"

  echo "mon script est : ${script}, les variables sonts : ${variable} et mon repertoire de sortie est : ${repertoire_sortie}"
  echo qsub ${script} -v ${variable} -q remod -l ncpus=8 -o ${repertoire_sortie} -e ${repertoire_sortie}
}

## Exemple pour lancer les fonctions
#parler_cluster_Torq script1 var8moi /Bureau/dossier

## Explication des parametres pour le cluster I2BC (PBS)
# -q : notre noeud sur lequel on veut travailler
# -l ncpus : nombre de cpus (un peu comme coeur)
# -v : liste des variables
# -o : repertoire de sortie pour la sortie standard
# -e : repertoire de sortie pour la sortie erreur


function parler_cluster_Slurm() {
  local script="${1}"
  local variable="${2}"

  echo "mon script est : ${script} et les variables sont : ${variable}"
  echo sbatch ${script} --export ${variable} --cpus-per-task 8 --output "${script##*/}.%j.out" --error "${script##*/}.%j.err"
}


## Exemple pour lancer les fonctions
#parler_cluster_Slurm script_slurm mes_var

## Explication des parametres pour le cluster SLURM
# --export : contient mes variables
# -o : le nom du fichier de la sortie standard
#	-e : le nom du fichier de la sortie erreur
