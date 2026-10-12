#!/bin/bash

# Sensitivity to missing-data handling: imputed data vs. complete cases
SETUP="real_imputed"

CENS_MODEL_TYPE_LIST=("grf")
MODEL_TYPE_LIST=("grf" "cox" "rsf" "xgb")
N_TRAIN_LIST=(5000)
N_CAL_LIST=(1000 10000)
M_LIST=(100 1000)
TIME_LIST=(3 6 9)
IMPUTED_LIST=(1 0)        # 1 = imputed data, 0 = complete cases
SCRIPT="analysis_imputed.sh"
MEMO=20G

# Slurm parameters
TIME=00-04:00:00
ORDP="sbatch --mem=$MEMO --nodes=1 --ntasks=1 --cpus-per-task=1 --time=$TIME --account=sesia_1123 --partition=main"

LOGS="logs/$SETUP"
OUT_DIR="results/$SETUP"
mkdir -p $LOGS $OUT_DIR

for IMPUTED in "${IMPUTED_LIST[@]}"; do
  for M in "${M_LIST[@]}"; do
    for N_CAL in "${N_CAL_LIST[@]}"; do
      for N_TRAIN in "${N_TRAIN_LIST[@]}"; do
        for CENS_MODEL_TYPE in "${CENS_MODEL_TYPE_LIST[@]}"; do
          for MODEL_TYPE in "${MODEL_TYPE_LIST[@]}"; do
            for TIME_T0 in "${TIME_LIST[@]}"; do

              # Must match output_file in the R script
              JOBN="cens_${CENS_MODEL_TYPE}_surv_${MODEL_TYPE}_train${N_TRAIN}_cal${N_CAL}_m${M}_time${TIME_T0}_imp${IMPUTED}.txt"
              if [[ -f $OUT_DIR/$JOBN ]]; then
                continue
              fi

              ARGS="$SETUP $CENS_MODEL_TYPE $MODEL_TYPE $N_TRAIN $N_CAL $M $TIME_T0 $IMPUTED"
              ORD="$ORDP -J $JOBN -o $LOGS/$JOBN.out -e $LOGS/$JOBN.err ./$SCRIPT $ARGS"
              echo $ORD
              $ORD

            done
          done
        done
      done
    done
  done
done
