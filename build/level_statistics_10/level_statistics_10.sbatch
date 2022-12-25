#!/bin/bash
# Job name:
#SBATCH --job-name=level_statistics_10
#
# Partition:
#SBATCH --partition=normal
#
# Processors:
#SBATCH --ntasks=256
#SBATCH --nodes=16
#
# Wall clock limit:
#SBATCH --time=4-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxd190113@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o level_statistics_10.o%j
#------------------------------------------------------

echo "Starting job level_statistics_10"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=level_statistics_10.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job level_statistics_10 complete"
