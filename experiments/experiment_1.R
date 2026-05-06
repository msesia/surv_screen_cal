# Load required libraries
library(tidyverse)

## Flag to determine if input should be parsed from command line
parse_input <- TRUE

source("../R/load_methods.R")

######################
## Input parameters ##
######################

if(parse_input) {
    ## Reading command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    ## Checking if the correct number of arguments is provided
    if (length(args) < 9) {
        stop("Insufficient arguments provided. Expected 9 arguments.")
    }
    ## Assigning command line arguments to variables
    setup <- args[1]
    real_data <- as.integer(args[2])
    gen_model_type <- args[3]
    surv_model_type <- args[4]
    num_samples_train <- as.integer(args[5])
    num_samples_cal <- as.integer(args[6])
    num_samples_test <- as.integer(args[7])
    screening_time <- as.numeric(args[8])
    batch <- as.integer(args[9])
} else {
    setup <- "v0"
    real_data <- 0
    gen_model_type <- "grf"
    surv_model_type <- "grf"
    num_samples_train <- 5000
    num_samples_cal <- 100
    num_samples_test <- 100
    screening_time <- 6
    batch <- 1
}

## DO NOT CHANGE THIS
screening_crit <- "low risk"

## Significance level
alpha <- 0.05

## Censoring model type
cens_model_type <- "grf"

## Do not use weights (default: TRUE, use FALSE only for debugging)
use_weights <- TRUE

## Number of boostrap samples
B_boot <- 500

## Number of repetitions
batch_size <- 20

####################
## Prepare output ##
####################

## Store important parameters including model types
header <- tibble(real_data = real_data,
                 gen_model_type = gen_model_type,
                 surv_model_type = surv_model_type,
                 cens_model_type = cens_model_type,
                 n_train = num_samples_train,
                 n_cal = num_samples_cal,
                 n_test = num_samples_test,
                 alpha = alpha,
                 screening_time = screening_time,
                 batch = batch)

## Generate a unique and interpretable file name based on the input parameters
output_file <- paste0("results/", setup, "/",
                      "real_", real_data,
                      "_gen_", gen_model_type,
                      "_surv_", surv_model_type,
                      "_train", num_samples_train,
                      "_cal", num_samples_cal,
                      "_test", num_samples_test,
                      "_time", screening_time,
                      "_batch", batch, ".txt")

## Print the output file name to verify
cat("Output file name:", output_file, "\n")


###############################################################
## Load the raw data and initialize semi-synthetic generator ##
###############################################################

load_raw_data <- function() {
    ## NOTE: data set is not provided
    combined <- readRDS("../data/fhrd_data_full.rds")
    combined <- as_tibble(do.call("rbind", combined))
    combined <- combined %>%
        mutate(time=month, status=event) %>%
        select(-month, -event) %>%
        select(time, status, everything())
    colnames(combined) <- c("time", "status", paste("X", seq(1,ncol(combined)-2), sep=""))
    return(combined)
}
data.raw <- load_raw_data()

## Instantiate generator (trains models once)
if(real_data) {
    data.generator <- RealDataGenerator$new(data = data.raw)
} else {
    data.generator <- SemiSyntheticDataGenerator$new(data = data.raw, surv_model_type = gen_model_type, cens_model_type = gen_model_type)
}

###################################################
## Instantiate the survival and censoring models ##
###################################################

surv_model <- init_surv_model(surv_model_type)
surv_model_large <- init_surv_model(surv_model_type)

# Instantiate censoring model based on the specified type
cens_base_model <- init_censoring_model(cens_model_type)
cens_model <- CensoringModel$new(model = cens_base_model)



analyze_data <- function(data.train, data.cal, data.test, surv_model, cens_model, data.test.oracle) {
    ## Fit the survival model on the training data
    surv_model$fit(Surv(time, status) ~ ., data = data.train)

    ## Fit the survival model on all training and calibration data
    data.supervised <- rbind(data.train, data.cal)
    surv_model_large$fit(Surv(time, status) ~ ., data = data.supervised)

    ## Fit the censoring model on training data
    cens_model$fit(data = data.train)

    ## Compute conformity scores for calibration data
    scores.cal <- as.numeric(surv_model$predict(data.cal, screening_time)$predictions)

    ## Compute conformity scores for test data
    scores.test <- as.numeric(surv_model$predict(data.test, screening_time)$predictions)

    ## Evaluate survival probabilities for estimated censoring distribution, for ET-IPCW weights
    prob_O.cal <- compute_censoring_survival(data.cal, cens_model, data.cal$time, data.cal$status, method = "et")

    lambda_seq <- seq(0.5, 1, by=0.1)

    df_eval <- lapply(lambda_seq, function(lambda) {
        ## Select test patients
        sel_test <- which(scores.test >= lambda)
        ## Evaluate performance
        eval <- evaluate_yield_and_event_rate(data.test.oracle, selected_idx = sel_test, t0 = screening_time)
        ## Return one-row tibble
        tibble::tibble(lambda = lambda, num_sel = eval$num_selected,
                       yield = eval$yield, ppv = 1 - eval$event_rate,
                       )
    })   
    df_eval <- bind_rows(df_eval)

    ## Compute calibration bands
    O.cal <- as.integer(data.cal$status)
    E.cal <- as.integer(data.cal$time <= screening_time)
    E.cal[O.cal == 0L] <- NA_integer_

    ## Boostrap (pointwise)
    band_lambda_boot <- purrr::map_dfr(lambda_seq, function(lambda) {
        S <- as.integer(scores.cal >= lambda)
        estim <- bootstrap_pointwise(
            S = S, E = E.cal, O = O.cal, prob_O = prob_O.cal,
            B = B_boot, m = num_samples_test, alpha = alpha, seed = 1,
            ipcw_method="ht"
        )
        make_lambda_row(lambda, S, estim, E, O, prob_O, num_samples_test)
    })
    estim_boot_sim <- bootstrap_simultaneous(scores = scores.cal, E = E.cal, O = O.cal, prob_O = prob_O.cal,
                                             B = B_boot, m = num_samples_test, alpha = alpha, seed = 1,
                                             ipcw_method="ht", lambda_seq=lambda_seq)    

    df_bands <- rbind(
        convert_band_to_long_targets(estim = estim_boot_sim, band_method = "Bootstrap-sim", m = num_samples_test),
        convert_pointwise_to_long_targets(band_df = band_lambda_boot, band_method = "Bootstrap-point", m = num_samples_test)
    )
    
    df_joined <- dplyr::left_join(df_eval, df_bands, by = "lambda")
    df_joined    
}


#######################################
## Define function to run experiment ##
#######################################

run_experiment <- function(random.state) {
    set.seed(random.state)

    ## Generate training, calibration, and test data
    data.synthetic.oracle <- data.generator$sample(shuffle=TRUE, return_oracle=TRUE)$oracle
    splits <- split_data_n(data.synthetic.oracle, n_train = num_samples_train, n_cal = num_samples_cal, n_test = num_samples_test)
    data.train.oracle <- splits$data.train
    data.cal.oracle <- splits$data.cal
    data.test.oracle <- splits$data.test

    ## Remove true event and censoring times from the data (right-censoring)
    data.train <- data.train.oracle |> select(-event_time, -censoring_time)
    data.cal <- data.cal.oracle |> select(-event_time, -censoring_time)
    data.test <- data.test.oracle |> select(-event_time, -censoring_time)

    ## Run analysis
    results <- analyze_data(data.train, data.cal, data.test, surv_model, cens_model, data.test.oracle)

    return(results)
}


## Function to run multiple experiments and gather results
## Args:
##   batch_size: Number of repetitions for each experimental setting
## Returns:
##   A tibble containing the combined results of all experiments
run_multiple_experiments <- function(batch_size) {
    results_df <- data.frame()  # Initialize an empty data frame to store cumulative results

    # Print a progress bar header
    cat("Running experiments\n")
    pb <- txtProgressBar(min = 0, max = batch_size, style = 3)  # Initialize progress bar

    # Loop over each repetition
    for (i in 1:batch_size) {
        random.state <- batch*1000 + i
        res <- run_experiment(random.state)  # Run experiment and get the result
        
        ## Combine the results with experiment metadata
        result_df <- tibble(Seed = random.state) |> cbind(header) |> cbind(res)

        # Add the result to the cumulative data frame
        results_df <- rbind(results_df, result_df) %>% as_tibble()

        # Write the cumulative results to the CSV file
        write.csv(results_df, output_file, row.names = FALSE)

        setTxtProgressBar(pb, i)  # Update progress bar
    }

    close(pb)  # Close the progress bar

    return(results_df)  # Return the cumulative results data frame
}

#####################
## Run experiments ##
#####################

## Run the experiments with specified parameters
results <- run_multiple_experiments(batch_size)

