#!/bin/bash
# Job name:
#SBATCH --job-name=new_KLd_10
#
# Partition:
#SBATCH --partition=CMT
#
# Processors:
#SBATCH --ntasks=32
#SBATCH --nodes=2
#
# Wall clock limit:
#SBATCH --time=7-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxd190113@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o new_KLd_10.o%j
#------------------------------------------------------

echo "Starting job new_KLd_10"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=new_KLd_10.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job new_KLd_10 complete"
