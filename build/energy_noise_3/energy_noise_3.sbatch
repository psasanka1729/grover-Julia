#!/bin/bash
# Job name:
#SBATCH --job-name=energy_noise_3
#
# Partition:
#SBATCH --partition=normal
#
# Processors:
#SBATCH --ntasks=16
#SBATCH --nodes=1
#
# Wall clock limit:
#SBATCH --time=96:00:00
#
# Email at end:
#SBATCH --mail-user=sxd190113@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o energy_noise_3.o%j
#------------------------------------------------------

echo "Starting job energy_noise_3"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=energy_noise_3.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job energy_noise_3 complete"
