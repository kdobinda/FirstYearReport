#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# run_sim_2026May30.sh
#
# Usage examples:
#
#   # Fast methods (1,2,3,5) for p=100:
#   P_VAL=100 MGROUP=fast RUN_ID=run_May30 sbatch run_sim_2026May30.sh
#
#   # DAGSLAM only for p=100 (separate array so it doesn't block fast jobs):
#   P_VAL=100 MGROUP=slow RUN_ID=run_May30 sbatch run_sim_2026May30.sh
#
#   # Fast methods for p=200:
#   P_VAL=200 MGROUP=fast RUN_ID=run_May30 sbatch run_sim_2026May30.sh
#
# The fast and slow arrays can run simultaneously -- submit both independently.
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --partition=multicore
#SBATCH --ntasks=1
#SBATCH --array=1-100
#SBATCH --job-name=mixed_sim
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err

# ── Resource allocation differs by method group ───────────────────────────────
# Controlled via MGROUP env var at submission time.
# fast: parallelises methods 1/2/5 across 8 CPUs; short wall-time
# slow: DAGSLAM is single-threaded Python; request fewer CPUs but more time

if [[ "${MGROUP:-fast}" == "slow" ]]; then
    # DAGSLAM: single-threaded, long wall-time, modest RAM
    #SBATCH --cpus-per-task=2
    #SBATCH --time=12:00:00
    #SBATCH --mem=16G
else
    # Fast methods: parallelised mclapply over 8 CPUs
    #SBATCH --cpus-per-task=8
    #SBATCH --time=04:00:00
    #SBATCH --mem=32G
fi

P_VAL=${P_VAL:-100}
MGROUP=${MGROUP:-fast}
RUN_ID=${RUN_ID:-"default_run"}

start_time=$SECONDS

# ── Scratch workspace ─────────────────────────────────────────────────────────
SCRATCH_DIR="$HOME/scratch/mixed_sim_runs/$RUN_ID"
mkdir -p "$SCRATCH_DIR"
cd "$SCRATCH_DIR"

# ── Environment ───────────────────────────────────────────────────────────────
module load apps/gcc/R/4.4.2
module load python/3.13.1
export JAVA_HOME="/mnt/iusers01/bk01-icvs/z55517kd/mixed_sim/jdk-21.0.3+9"
export PATH="$JAVA_HOME/bin:$PATH"

# Limit OpenBLAS / MKL threads so they don't fight with mclapply workers
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
duration=$((end_time - start_time))
echo "Finished task $SLURM_ARRAY_TASK_ID in ${duration}s at: $(date)"
