#!/bin/bash
# Job name:
#SBATCH --job-name=special_states_1
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
#SBATCH -o special_states_1.o%j
#------------------------------------------------------

echo "Starting job special_states_1"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=special_states_1.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job special_states_1 complete"
