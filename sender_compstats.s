#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name=cddp_sender
#SBATCH --time=00:05:00
#SBATCH --mail-type=END
#SBATCH --mail-user=froldan@nyu.edu
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err


# Load modules
module purge
module load rclone/1.38

# rm cddp_compstats -rf
mkdir -p cddp_compstats/
cp -n ./ct_1.jld cddp_compstats/
cp -n ./ct_1_temp.jld cddp_compstats/
cp -n ./ct_opt.jld cddp_compstats/
cd cddp_compstats

N=20
rclone copy remote_dropbox:NYU/InflAPStion/Codes/compstats_cddp.s ./
echo start

touch compstats.txt
sbatch --export=N=$N --array=1-$N compstats_cddp.s

echo all jobs sent