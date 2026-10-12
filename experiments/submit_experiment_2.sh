#!/bin/bash

# Parameters
SETUP="v0h"

# Defaults (overridden per setup below)
CENS_TILT_LIST=(0)   # strength of the covariate-dependent censoring tilt
                     # (log hazard ratio per SD of X32); only used by "v0t"
CENS_SHIFT_LIST=(0)  # log-hazard shift raising the overall censoring rate; "v0t" only

if [[ $SETUP == "v0t" ]]; then
  # Whether to use real data
  REAL_LIST=(0)
  # Generative model types
  GEN_MODEL_TYPE_LIST=("grf")
  # Survival model types (use same model type for censoring model)
  MODEL_TYPE_LIST=("grf") # "cox" "grf" "xgb"
  # List of training sample sizes for censoring model
  N_TRAIN_LIST=(2000)
  # List of calibration sample sizes
  N_CAL_LIST=(1000)
  # List of test sample sizes
  N_TEST_LIST=(100 1000)
  # List of time points
  TIME_LIST=(3)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 5)
  # Tilt strengths (0 = untilted control point) and shifts
  CENS_TILT_LIST=(0)
  CENS_SHIFT_LIST=(0 0.5 1 1.5 2 2.5 3 3.5)
  # Memory
  MEMO=20G

elif [[ $SETUP == "v0h" ]]; then
  # Whether to use real data
  REAL_LIST=(0)
  # Generative model types
  GEN_MODEL_TYPE_LIST=("grf")
  # Survival model types (use same model type for censoring model)
  MODEL_TYPE_LIST=("grf") # "cox" "grf" "xgb"
  # List of training sample sizes for censoring model
  N_TRAIN_LIST=(2000)
  # List of calibration sample sizes
  N_CAL_LIST=(1000)
  # List of test sample sizes
  N_TEST_LIST=(100 1000)
  # List of time points
  TIME_LIST=(3)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 5)
  # Tilt strengths (0 = untilted control point) and shifts
  CENS_TILT_LIST=(-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3)
  CENS_SHIFT_LIST=(0)
  # Memory
  MEMO=20G

fi

# Slurm parameters
TIME=00-01:00:00                    # Time required (20 m)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1123 --partition=main"

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS
mkdir -p $LOGS"/"$SETUP

OUT_DIR="results"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/"$SETUP

# Loop over configurations
for BATCH in $BATCH_LIST; do
  for CENS_SHIFT in "${CENS_SHIFT_LIST[@]}"; do
    for CENS_TILT in "${CENS_TILT_LIST[@]}"; do
      for REAL in "${REAL_LIST[@]}"; do
        for N_TEST in "${N_TEST_LIST[@]}"; do
          for N_CAL in "${N_CAL_LIST[@]}"; do
            for N_TRAIN in "${N_TRAIN_LIST[@]}"; do
              for GEN_MODEL_TYPE in "${GEN_MODEL_TYPE_LIST[@]}"; do
                for MODEL_TYPE in "${MODEL_TYPE_LIST[@]}"; do
                  for TIME in "${TIME_LIST[@]}"; do

                    # File-name tag for tilt and shift; must match `mix_tag` in the R script
                    MIX_TAG="_tilt${CENS_TILT}_shift${CENS_SHIFT}"
                    
                    # Generate a unique and interpretable file name based on the input parameters
                    JOBN="real_${REAL}_gen_${GEN_MODEL_TYPE}_surv_${MODEL_TYPE}_train${N_TRAIN}${MIX_TAG}_cal${N_CAL}_test${N_TEST}_time${TIME}_batch${BATCH}.txt"
                    OUT_FILE=$OUT_DIR"/"$SETUP"/"$JOBN
                    COMPLETE=0

                    if [[ -f $OUT_FILE ]]; then
                      COMPLETE=1
                    fi

                    if [[ $COMPLETE -eq 0 ]]; then
                      # R script to be run with command line arguments (10th = tilt, 11th = shift)
                      SCRIPT="./experiment_2.sh $SETUP $REAL $GEN_MODEL_TYPE $MODEL_TYPE $N_TRAIN $N_CAL $N_TEST $TIME $BATCH $CENS_TILT $CENS_SHIFT"

                      # Define job name for this configuration
                      OUTF=$LOGS"/"$SETUP"/"$JOBN".out"
                      ERRF=$LOGS"/"$SETUP"/"$JOBN".err"

                      # Assemble slurm order for this job
                      ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" $SCRIPT"

                      # Print order
                      echo $ORD
                      # Submit order
                      $ORD
                      # Run command now
                      #./$SCRIPT

                    fi

                  done
                done
              done
            done
          done
        done
      done
    done
  done
done
