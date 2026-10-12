#!/bin/bash

# Parameters
SETUP="timing"

if [[ $SETUP == "timing" ]]; then
  # Survival model sets (names in MODEL_SETS of the R script)
  MODEL_TYPE_LIST=("grf")
  # List of calibration sample sizes
  N_CAL_LIST=(100 200 500 1000 2000 5000 10000 20000)
  # List of cohort sizes
  M_LIST=(100 1000 10000)
  # Sequence of batches
  BATCH_LIST=$(seq 1 1)
  # Memory
  MEMO=20G

fi

# Slurm parameters: one core, so timings are comparable
TIME=00-00:10:00
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1123 --partition=main"

LOGS="logs/"$SETUP
OUT_DIR="results/"$SETUP
mkdir -p $LOGS $OUT_DIR

for BATCH in $BATCH_LIST; do
  for M in "${M_LIST[@]}"; do
    for N_CAL in "${N_CAL_LIST[@]}"; do
      for MODEL_TYPE in "${MODEL_TYPE_LIST[@]}"; do

        JOBN="surv_${MODEL_TYPE}_cal${N_CAL}_m${M}_batch${BATCH}.txt"
        OUT_FILE=$OUT_DIR"/"$JOBN

        if [[ -f $OUT_FILE ]]; then
          echo "Skipping $JOBN (already complete)"
          continue
        fi

        SCRIPT="./experiment_timing.sh $SETUP $MODEL_TYPE $N_CAL $M $BATCH"
        OUTF=$LOGS"/"$JOBN".out"
        ERRF=$LOGS"/"$JOBN".err"

        ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" $SCRIPT"
        echo $ORD
        $ORD

      done
    done
  done
done
