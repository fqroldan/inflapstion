#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=7
#SBATCH --mem=60GB
#SBATCH --job-name=cddp_%a
#SBATCH --time=5-10:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=froldan@nyu.edu
#SBATCH --output=cddp_%A_%a.out
#SBATCH --error=cddp_%A_%a.err

x=$SLURM_ARRAY_TASK_ID

module purge
module load julia/1.3.0
module load rclone/1.38
echo $x
echo $N

mkdir -p run$x/Graphs/tests/
cd run$x/

git clone https://github.com/fqroldan/inflapstion.git
cd inflapstion

export JULIA_NUM_THREADS=7
julia comp_stats.jl $x $N
