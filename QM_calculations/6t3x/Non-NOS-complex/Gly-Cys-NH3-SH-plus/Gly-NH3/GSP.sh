#!/bin/bash
#SBATCH -p elephant 
#SBATCH -N 1   # node count
#SBATCH --mem=8G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12

WORK_DIR=`pwd`
echo $WORK_DIR > tmp
sed s/home/scratch/g < tmp > tmp2
SCR_DIR=`cat tmp2`
rm tmp tmp2
mkdir -p $SCR_DIR

echo $HOSTNAME
echo $WORK_DIR
echo $SCR_DIR
head -n 3 /proc/meminfo

source /etc/profile.d/modules.sh
. /etc/profile
module load g16_a03

g16 Gly-NH3.gjf

