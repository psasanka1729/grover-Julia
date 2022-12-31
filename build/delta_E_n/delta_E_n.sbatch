#!/bin/bash
# Job name:
#SBATCH --job-name=delta_E_n
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
#SBATCH -o delta_E_n.o%j
#------------------------------------------------------

echo "Starting job delta_E_n"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=delta_E_n.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job delta_E_n complete"
