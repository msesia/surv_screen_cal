#!/bin/bash

# Parameters
SETUP="real_imputed"

# Defaults (overridden per setup below)
BATCH_LIST=$(seq 1 1)   # batches (random seeds); only used when USE_BATCH=1
USE_BATCH=0             # 1 = pass batch to the script and add _batch<b> to file names

if [[ $SETUP == "real_bands_1" ]]; then
  CENS_MODEL_TYPE_LIST=("grf")
  MODEL_TYPE_LIST=("grf")
  N_TRAIN_LIST=(5000)
  N_CAL_LIST=(10000)
  M_LIST=(100 1000)
  TIME_LIST=(3 6 9)
  SCRIPT="analysis_bands.sh"
  MEMO=20G

elif [[ $SETUP == "real_bands_2" ]]; then
  CENS_MODEL_TYPE_LIST=("grf")
  MODEL_TYPE_LIST=("grf" "cox" "rsf" "xgb" "mix")
  N_TRAIN_LIST=(5000)
  N_CAL_LIST=(10000)
  M_LIST=(100)
  TIME_LIST=(6)
  SCRIPT="analysis_bands.sh"
  MEMO=20G

elif [[ $SETUP == "real_selection_1" ]]; then
  CENS_MODEL_TYPE_LIST=("grf")
  MODEL_TYPE_LIST=("grf" "cox" "rsf" "xgb")
  N_TRAIN_LIST=(5000)
  N_CAL_LIST=(10000)
  M_LIST=(100 1000)
  TIME_LIST=(3 6 9)
  SCRIPT="analysis_selection.sh"
  BATCH_LIST=$(seq 0 0)
  USE_BATCH=1
  MEMO=20G

elif [[ $SETUP == "real_selection_2" ]]; then
  CENS_MODEL_TYPE_LIST=("grf")
  MODEL_TYPE_LIST=("mix")
  N_TRAIN_LIST=(5000)
  N_CAL_LIST=(10000)
  M_LIST=(100 1000)
  TIME_LIST=(3 6 9)
  SCRIPT="analysis_selection.sh"
  BATCH_LIST=$(seq 0 0)
  USE_BATCH=1
  MEMO=20G

elif [[ $SETUP == "real_selection_3" ]]; then
  CENS_MODEL_TYPE_LIST=("grf" "cox" "rsf" "xgb")
  MODEL_TYPE_LIST=("grf" "xgb")
  N_TRAIN_LIST=(5000)
  N_CAL_LIST=(10000)
  M_LIST=(100 1000)
  TIME_LIST=(3 6 9)
  SCRIPT="analysis_selection.sh"
  BATCH_LIST=$(seq 0 0)
  USE_BATCH=1
  MEMO=20G

elif [[ $SETUP == "real_selection_1_cal" ]]; then
  CENS_MODEL_TYPE_LIST=("grf" "cox" "rsf" "xgb")
  MODEL_TYPE_LIST=("grf")
  N_TRAIN_LIST=(5000)
  N_CAL_LIST=(100 200 500 1000 2000 5000 10000)
  M_LIST=(1000)
  TIME_LIST=(3)
  SCRIPT="analysis_selection.sh"
  BATCH_LIST=$(seq 1 1)
  USE_BATCH=1
  MEMO=20G

elif [[ $SETUP == "real_selection_1b_cal" ]]; then
  CENS_MODEL_TYPE_LIST=("grf")
  MODEL_TYPE_LIST=("grf")
  N_TRAIN_LIST=(5000)
  N_CAL_LIST=(100 200 500 1000 2000 5000 10000)
  M_LIST=(100 1000)
  TIME_LIST=(3 6 9)
  SCRIPT="analysis_selection.sh"
  BATCH_LIST=$(seq 1 5)
  USE_BATCH=1
  MEMO=20G


elif [[ $SETUP == "real_varsel_1" ]]; then
  CENS_MODEL_TYPE_LIST=("grf")
  MODEL_TYPE_LIST=("grf" "cox" "rsf" "xgb")
  N_TRAIN_LIST=(5000)
  N_CAL_LIST=(10000)
  M_LIST=(100 1000)
  TIME_LIST=(3 6 9)
  BATCH_LIST=$(seq 1 20)
  USE_BATCH=1
  SCRIPT="analysis_varsel.sh"
  MEMO=20G

elif [[ $SETUP == "real_selection_varsel_1" ]]; then
  CENS_MODEL_TYPE_LIST=("grf")
  MODEL_TYPE_LIST=("grf" "cox" "rsf" "xgb")
  N_TRAIN_LIST=(5000)
  N_CAL_LIST=(10000)
  M_LIST=(100 1000)
  TIME_LIST=(3 6 9)
  SCRIPT="analysis_selection_varsel.sh"
  MEMO=20G

fi

# Slurm parameters
TIME=00-04:00:00                    # Time required (4 h)
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
  for M in "${M_LIST[@]}"; do
    for N_CAL in "${N_CAL_LIST[@]}"; do
      for N_TRAIN in "${N_TRAIN_LIST[@]}"; do
        for CENS_MODEL_TYPE in "${CENS_MODEL_TYPE_LIST[@]}"; do
          for MODEL_TYPE in "${MODEL_TYPE_LIST[@]}"; do
            for TIME_T0 in "${TIME_LIST[@]}"; do

              # Generate a unique and interpretable file name based on the input parameters
              if [[ $USE_BATCH -eq 1 ]]; then
                JOBN="cens_${CENS_MODEL_TYPE}_surv_${MODEL_TYPE}_train${N_TRAIN}_cal${N_CAL}_m${M}_time${TIME_T0}_batch${BATCH}.txt"
                ARGS="$SETUP $CENS_MODEL_TYPE $MODEL_TYPE $N_TRAIN $N_CAL $M $TIME_T0 $BATCH"
              else
                JOBN="cens_${CENS_MODEL_TYPE}_surv_${MODEL_TYPE}_train${N_TRAIN}_cal${N_CAL}_m${M}_time${TIME_T0}.txt"
                ARGS="$SETUP $CENS_MODEL_TYPE $MODEL_TYPE $N_TRAIN $N_CAL $M $TIME_T0"
              fi
              OUT_FILE=$OUT_DIR"/"$SETUP"/"$JOBN
              COMPLETE=0

              if [[ -f $OUT_FILE ]]; then
                COMPLETE=1
              fi

              if [[ $COMPLETE -eq 0 ]]; then
                # R script to be run with command line arguments
                COMMAND="./$SCRIPT $ARGS"

                # Define job name for this configuration
                OUTF=$LOGS"/"$SETUP"/"$JOBN".out"
                ERRF=$LOGS"/"$SETUP"/"$JOBN".err"

                # Assemble slurm order for this job
                ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" $COMMAND"

                # Print order
                echo $ORD
                # Submit order
                $ORD
                # Run command now
                #./$COMMAND

              fi

            done
          done
        done
      done
    done
  done
done
