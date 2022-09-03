#!/bin/bash
# Job name:
#SBATCH --job-name=run2_03092022
#
# Partition:
#SBATCH --partition=CMT
#
# Processors:
#SBATCH --ntasks=16
#SBATCH --nodes=1
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
# Email at end:
#SBATCH --mail-user=sxd190113@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run2_03092022.o%j
#------------------------------------------------------

echo "Starting job run2_03092022"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run2_03092022.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job run2_03092022 complete"
