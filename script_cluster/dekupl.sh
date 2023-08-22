#!/bin/bash

##Configuration SLURM ########
##on définit le nom du job 
#SBATCH --job-name=dekuplN

##définir la partition utilisée
#SBATCH -p normal

##définir le nombre de  à utiliser 
#SBATCH -c 9

###SBATCH --output=OutDEkuplN.txt

#SBATCH --mem 54
##SBATCH --array=1-19%7 ##index table between 0 to 19 with 4 parallele jobs

##Creer le répertoire perso 

cd /scratch
mkdir -p ${USER}_${SLURM_JOB_ID}_DEKUPL
cd ${USER}_${SLURM_JOB_ID}_DEKUPL
mkdir -p output
cd output

#PARAMETERS
image_singularity=$1
config_file=$2
fasta=$3

echo "COPYING DATA" 

mkdir -p input/
cd input
scp ${image_singularity} .
scp ${config_file}

mkdir -p data/
cd data

scp ${fasta} .
cd ../

#scp nas:/home/boubred/data/image/dekupl-run.img image/
#scp nas:/home/boubred/data/image/config.json .


##Charger les logiciel à utiliser 
module load system/singularity/3.6.0

#scp nas3:/data3/projects/integration/rnaseq_riz/output/RC/RC*cutadapt*.fq.gz .
#scp nas3:/data3/projects/integration/rnaseq_riz/output/RS/RS*cutadapt*fq.gz . 



fasta = $(basename ${fasta})
image_singularity = $(basename ${image_singularity})
config_file = $(basename ${config_file})

path_output = $(realpath ${fasta})
path_output = $(dirname ${path_output})

singularity run input/${image_singularity} --configfile input/${config_file} -j 9 --resources ram=20000 -p
echo -e "\n#----------------------------------------------#\n"
echo -e "\n################# DE-KUPL ################\n"
echo -e "\n#----------------------------------------------#\n"

rm -r input/
rm -r data/

tree
#TRANSFERT DES RESULTATS
scp -r output ${path_output}

cd ../
rm -rf /scratch/DEKUPL_boubred
echo -e "\n#----------------------------------------------#\n"
echo -e "\n################# Suppr /scratch ################\n"
echo -e "\n#----------------------------------------------#\n"
