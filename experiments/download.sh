DATA=$1

mkdir -p results_hpc

rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/merck_project_exploration/code_exploration/analysis/results/* results_hpc/

#rsync -auv ../../data/sample_data_full.rds  sesia@discovery.usc.edu:/home1/sesia/Workspace/screening_calibration_private/data/
