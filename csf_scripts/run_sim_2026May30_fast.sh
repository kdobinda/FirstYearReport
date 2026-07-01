#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# run_sim_2026May30_fast.sh  --  methods 1, 2, 3, 5
#
# Submit for each p:
#   P_VAL=200 RUN_ID=run_May30 sbatch run_sim_2026May30_fast.sh
#   P_VAL=300 RUN_ID=run_May30 sbatch run_sim_2026May30_fast.sh
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH --partition=multicore
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=1-100
#SBATCH --job-name=mixsim_fast
#SBATCH --output=logs/fast_%A_%a.out
#SBATCH --error=logs/fast_%A_%a.err

P_VAL=${P_VAL:-200}
RUN_ID=${RUN_ID:-"default_run"}
MGROUP="fast"

start_time=$SECONDS

SCRATCH_DIR="$HOME/scratch/mixed_sim_runs/$RUN_ID"
mkdir -p "$SCRATCH_DIR"
cd "$SCRATCH_DIR"

module load apps/gcc/R/4.4.2
module load python/3.13.1
export JAVA_HOME="/mnt/iusers01/bk01-icvs/z55517kd/mixed_sim/jdk-21.0.3+9"
export PATH="$JAVA_HOME/bin:$PATH"

# Keep BLAS single-threaded so mclapply workers don't fight each other
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "====================================================="
echo "Array task : $SLURM_ARRAY_TASK_ID"
echo "p          : $P_VAL"
echo "Method grp : $MGROUP"
echo "Run ID     : $RUN_ID"
echo "CPUs       : $SLURM_CPUS_PER_TASK"
echo "Node       : $SLURM_NODELIST"
echo "====================================================="

Rscript $HOME/mixed_sim/simulation_run_2026May30_csf.R \
    $SLURM_ARRAY_TASK_ID \
    $P_VAL \
    $MGROUP \
    $RUN_ID

end_time=$SECONDS
echo "Finished task $SLURM_ARRAY_TASK_ID in $((end_time - start_time))s at: $(date)"
