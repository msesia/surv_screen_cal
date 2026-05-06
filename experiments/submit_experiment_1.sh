#!/bin/bash

# Parameters
SETUP="v0"

if [[ $SETUP == "v0" ]]; then
  # Whether to use real data
  REAL_LIST=(0)
  # Generative model types
  GEN_MODEL_TYPE_LIST=("grf")
  # Survival model types (use same model type for censoring model)
  MODEL_TYPE_LIST=("grf") # "cox" "grf2")
#  MODEL_TYPE_LIST=("grf")
  # List of training sample sizes for censoring model
  N_TRAIN_LIST=(5000)
  # List of calibration sample sizes
  N_CAL_LIST=(100)
  # List of test sample sizes
  N_TEST_LIST=(100)
  # List of time points
  TIME_LIST=(3)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 1)
  # Memory
  MEMO=5G

elif [[ $SETUP == "v1" ]]; then
  # Whether to use real data
  REAL_LIST=(1)
  # Generative model types
  GEN_MODEL_TYPE_LIST=("grf")
  # Survival model types (use same model type for censoring model)
  MODEL_TYPE_LIST=("grf" "cox" "xgb")
  # List of training sample sizes for censoring model
  N_TRAIN_LIST=(5000)
  # List of calibration sample sizes
  N_CAL_LIST=(100 1000)
  # List of test sample sizes
  N_TEST_LIST=(100)
  # List of time points
  TIME_LIST=(2 3 6 9)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 5)
  # Memory
  MEMO=5G

fi

# Slurm parameters
TIME=00-00:20:00                    # Time required (20 m)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --partition=main"

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS
mkdir -p $LOGS"/"$SETUP

OUT_DIR="results"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/"$SETUP

# Loop over configurations
for BATCH in $BATCH_LIST; do
  for REAL in "${REAL_LIST[@]}"; do
    for N_TEST in "${N_TEST_LIST[@]}"; do
      for N_CAL in "${N_CAL_LIST[@]}"; do
        for N_TRAIN in "${N_TRAIN_LIST[@]}"; do
          for GEN_MODEL_TYPE in "${GEN_MODEL_TYPE_LIST[@]}"; do
            for MODEL_TYPE in "${MODEL_TYPE_LIST[@]}"; do
              for TIME in "${TIME_LIST[@]}"; do
                
                # Generate a unique and interpretable file name based on the input parameters
                JOBN="real_${REAL}_gen_${GEN_MODEL_TYPE}_surv_${MODEL_TYPE}_train${N_TRAIN}_cal${N_CAL}_test${N_TEST}_time${TIME}_batch${BATCH}.txt"
                OUT_FILE=$OUT_DIR"/"$SETUP"/"$JOBN
                #ls $OUT_FILE
                COMPLETE=0

                if [[ -f $OUT_FILE ]]; then
                  COMPLETE=1
                fi

                if [[ $COMPLETE -eq 0 ]]; then
                  # R script to be run with command line arguments
                  SCRIPT="./experiment_1.sh $SETUP $REAL $GEN_MODEL_TYPE $MODEL_TYPE $N_TRAIN $N_CAL $N_TEST $TIME $BATCH"

                  # Define job name for this configuration
                  OUTF=$LOGS"/"$JOBN".out"
                  ERRF=$LOGS"/"$JOBN".err"

                  # Assemble slurm order for this job
                  ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" $SCRIPT"

                  # Print order
                  echo $ORD
                  # Submit order
                  #$ORD
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
