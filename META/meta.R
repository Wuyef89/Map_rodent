# Load required packages
library(metafor)
library(car) 
library(openxlsx) 
library(tidyverse) 

# Set global paths and parameters
base_path <- "/data/home/test01/meta"
input_file <- file.path(base_path, "Supplementary Table 2.xlsx")
output_path <- file.path(base_path, "results")
set.seed(123)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Define global control object for rma.mv function
control <- list(
  optimizer = "optim",    # Optimizer
  optmethod = "L-BFGS-B", # Optimization method
  maxiter = 10000,        # Maximum iterations
  reltol = 1e-8,          # Relative convergence tolerance
  stepadj = 0.5           # Step adjustment factor
)

# Create output directory
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# --- Data loading and preprocessing ---
df <- read.xlsx(input_file, 1)
message("Original number of rows: ", nrow(df))

# Remove unnecessary columns
dat <- df[, setdiff(names(df), c("Original_habitat", "Original_rodents",
                                 "Chinese_rodent_name", "Specific_mixed_sample",
                                 "Testing_method", "Specific_testing_method",
                                 "Original_pathogens", "Standardization"))]

# Filter data: healthy status, non-missing NDVI and Time_interval, non-unknown habitat, non-unclassified samples
# And total Events for Pathogen_species >= 5
data <- dat %>%
  group_by(Pathogen_species) %>%
  filter(sum(Events) >= 5) %>%
  ungroup() %>%
  filter(Health_status == "Healthy",
         !is.na(NDVI), !is.na(Time_interval),
         Habitat_type != "Unknown habitat",
         !Sample_source %in% c("Unclassified sample"))
message("Number of rows after filtering: ", nrow(data))

# Define a function to set factor levels
define_factor_levels <- function(data, column_name, levels_order) {
  if (column_name %in% names(data)) {
    data[[column_name]] <- factor(data[[column_name]], levels = levels_order)
  } else {
    warning(paste("Column", column_name, "does not exist in the data for factor level definition."))
  }
  return(data)
}

# Define factor level order for each column
levels_list <- list(
  Habitat_type = c('Artificial habitat', 'Natural habitat', 'Mixed habitat'),
  Sample_source = c('Alimentary or intestinal sample', 'Blood or serum sample',
                    'Body surface sample', 'Brain sample','Diaphragm sample',
                    'Ear sample','Heart sample', 'Kidney sample',
                    'Liver sample','Lung or respiratory tract sample',
                    'Spleen sample','Urine sample', 'Mixed sample'),
  Testing_type = c('Morphological examination', 'Microbial culture',
                   'Immunological testing', 'Molecular testing'),
  Rodent_family = c('Caviidae', 'Chinchillidae', 'Cricetidae', 'Dipodidae', 'Hystricidae',
                    'Muridae', 'Myocastoridae', 'Nesomyidae', 'Sciuridae', 'Spalacidae'),
  Region = c('Central China', 'Eastern China', 'Northern China',
             'Northeastern China', 'Northwestern China',
             'Southern China', 'Southwestern China'),
  Time_interval = c('1950-2000', '2001-2010', '2011-2019', '2020-2023'),
  NDVI = c('Peak NDVI', 'Higher NDVI', 'Medium NDVI', 'Low NDVI')
)

# Define a function to filter data, ensuring each group has at least min_count observations
filter_data <- function(data, columns, min_count = 2) {
  for (column in columns) {
    if (column %in% names(data)) {
      data <- data %>%
        group_by(across(all_of(column))) %>%
        filter(n() >= min_count) %>%
        ungroup()
    } else {
      warning(paste("Column", column, "does not exist in the data."))
    }
  }
  return(data)
}

# Define datasets and columns to check
datasets <- list(
  data_a_original = data,  # All pathogens
  data_v_original = data %>% filter(Pathogen_type == "Viruses"),
  data_b_original = data %>% filter(Pathogen_type == "Bacteria"),
  data_p_original = data %>% filter(Pathogen_type == "Parasites")
)

columns_to_check <- c("Habitat_type", "Sample_source", "Testing_type",
                      "Rodent_family", "Region", "Time_interval", "NDVI")

# Apply factor level setting and filtering to each dataset, and define Study and Observation factors
datasets <- lapply(datasets, function(dataset) {
  for (col_name in names(levels_list)) {
    dataset <- define_factor_levels(dataset, col_name, levels_list[[col_name]])
  }
  # Define Study and Observation factors here to ensure uniqueness
  dataset$Observation <- factor(1:nrow(dataset))
  dataset$Study <- factor(dataset$Study)
  filter_data(dataset, columns_to_check)
})

## I² calculation function
i2 <- function(model) {
  # First check if the model is NULL or not a valid rma.mv object
  if (is.null(model) || !inherits(model, "rma.mv")) {
    warning("Model is NULL or not a valid 'rma.mv' object, returning NA for I2.")
    return(data.frame(I2_Overall = NA, I2_Level1 = NA, I2_Level2 = NA))
  }
  
  # Further check if the model contains the required components
  if (is.null(model$sigma2) || is.null(model$k) || is.null(model$p) || is.null(model$vi)) {
    warning("Model is missing essential components (sigma2, k, p, vi) for I2 calculation, returning NA.")
    return(data.frame(I2_Overall = NA, I2_Level1 = NA, I2_Level2 = NA))
  }
  
  # Safely attempt to get the model matrix. If it fails, it indicates a problem with the model structure.
  X <- tryCatch({
    model.matrix(model)
  }, error = function(e) {
    warning(paste("Could not extract model.matrix from model:", e$message, "Returning NULL for X and NA for I2."))
    return(NULL)
  })
  
  if (is.null(X)) {
    return(data.frame(I2_Overall = NA, I2_Level1 = NA, I2_Level2 = NA))
  }
  
  W <- diag(1 / model$vi)
  
  # Check if the matrix is invertible or has sufficient dimensions
  # Added check for ncol(X) == 0, as model.matrix can cause issues if only an intercept term
  if (ncol(X) == 0 || nrow(X) < ncol(X) || det(t(X) %*% W %*% X) == 0) {
    warning("Matrix for P is singular, ill-conditioned, or has insufficient dimensions. I2 cannot be calculated.")
    return(data.frame(I2_Overall = NA, I2_Level1 = NA, I2_Level2 = NA))
  }
  
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  
  # Ensure sigma2 is non-empty and non-zero to avoid division by zero or NaN
  if (length(model$sigma2) == 0 || all(model$sigma2 == 0, na.rm = TRUE)) { # Added na.rm for safety
    warning("sigma2 is empty or zero, I2 cannot be calculated.")
    return(data.frame(I2_Overall = NA, I2_Level1 = NA, I2_Level2 = NA))
  }
  
  denominator <- sum(model$sigma2) + (model$k - model$p) / sum(diag(P))
  
  if (denominator <= 0 || is.na(denominator) || is.infinite(denominator)) {
    I2_overall <- NA
    I2_levels <- rep(NA, length(model$sigma2))
    warning("Denominator for I2 calculation is non-positive or infinite, returning NA.")
  } else {
    I2_overall <- 100 * sum(model$sigma2) / denominator
    I2_levels <- 100 * model$sigma2 / denominator
  }
  
  # Prepare the returning data frame
  # Dynamically generate column names based on the number of random effects
  result_df <- data.frame(I2_Overall = round(I2_overall, 2))
  
  # Try to get the names of the random effects and use them for column names
  re_names <- if (!is.null(names(model$sigma2))) {
    names(model$sigma2)
  } else if (length(model$sigma2) == 2) { # Fallback to generic names if names are missing
    c("Study", "Observation")
  } else { # For any other number of levels
    paste0("Level", 1:length(model$sigma2))
  }
  
  for (j in seq_along(I2_levels)) {
    col_name <- re_names[j]
    # Handle cases where re_names might not be unique or valid, though less likely for metafor objects
    if (is.null(col_name) || is.na(col_name) || col_name == "") {
      col_name <- paste0("I2_Level", j)
    } else {
      col_name <- paste0("I2_", col_name)
    }
    result_df[[col_name]] <- round(I2_levels[j], 2)
  }
  
  # Ensure we have at least I2_Study and I2_Observation if they were expected but not present
  if (!("I2_Study" %in% names(result_df))) result_df$I2_Study <- NA
  if (!("I2_Observation" %in% names(result_df))) result_df$I2_Observation <- NA
  
  # Reorder columns for consistency if needed (optional, but good for output)
  result_df <- result_df %>%
    select(I2_Overall, starts_with("I2_Study"), starts_with("I2_Observation"), everything())
  
  return(result_df)
}

# Function to calculate effect size
calc_effect_size <- function(data) {
  # Ensure 'data' is a data frame
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  # Use escalc to calculate effect size and variance
  data <- data.frame(data, escalc(xi = data$Events, ni = data$Total, measure = "PLN"))
  data$original_ratio <- exp(data$yi)
  # Observation and Study factors are already handled in previous steps, no need to repeat here
  return(data)
}

# Calculate effect size for each dataset
datasets <- lapply(datasets, calc_effect_size)

# Reassign modified datasets
data_a_original <- datasets[[1]]
data_v_original <- datasets[[2]]
data_b_original <- datasets[[3]]
data_p_original <- datasets[[4]]

# --- Save datasets to different sub-sheets in Excel ---
output_excel_file <- file.path(output_path, paste0("filtered_datasets_", timestamp, ".xlsx"))

# Create a new Excel workbook
wb <- createWorkbook()

# Iterate through the datasets list and add each data frame to a different worksheet
for (name in names(datasets)) {
  addWorksheet(wb, sheetName = name) # Add worksheet
  writeData(wb, sheet = name, x = datasets[[name]]) # Write data
}

# Save the workbook
saveWorkbook(wb, output_excel_file, overwrite = TRUE)
message(paste0("All datasets successfully saved to: ", output_excel_file))


# --- Bootstrap resampling ---
boot.func <- function(dataset) {
  indices <- sample(1:nrow(dataset), replace = TRUE)
  dataset[indices, ] # Return the sampled data frame
}

# Perform Bootstrap to generate 1000 data frames
set.seed(123)
bootstrap_data_a <- lapply(1:1000, function(i) boot.func(data_a_original))
bootstrap_data_v <- lapply(1:1000, function(i) boot.func(data_v_original))
bootstrap_data_b <- lapply(1:1000, function(i) boot.func(data_b_original))
bootstrap_data_p <- lapply(1:1000, function(i) boot.func(data_p_original))

# --- VIF filtering and model formula generation ---
# Calculate GVIF and remove highly collinear variables
filter_high_vif <- function(data, threshold = 2.5) { # Threshold can be slightly higher, e.g., 2.5 or 5
  full_formula <- yi ~ Sample_source + Region + Habitat_type +
    Testing_type + Rodent_family + Time_interval + NDVI - 1 # -1 means no intercept, but meta-analysis usually needs an intercept
  
  # Check if enough data for lm fitting (at least one more observation than variables)
  if (nrow(data) < length(all.vars(full_formula)) - 1) {
    warning("Not enough data to fit lm model for VIF calculation.")
    return(NULL)
  }
  
  repeat {
    lm_model <- tryCatch(lm(full_formula, data = data, na.action = na.omit), # Ignore NA
                         error = function(e) {
                           message("Error in lm for VIF calculation: ", e$message)
                           return(NULL)
                         })
    if (is.null(lm_model)) return(NULL) # If lm fitting fails, return NULL directly
    
    # Check if VIF can be calculated (at least two variables)
    if (length(coef(lm_model)) < 2) break
    
    vif_values <- tryCatch(vif(lm_model), error = function(e) {
      message("Error in vif calculation: ", e$message)
      return(NULL)
    })
    if (is.null(vif_values) || (!is.numeric(vif_values) && !is.matrix(vif_values))) {
      warning("VIF calculation failed or returned invalid format.")
      break # Invalid VIF result, break loop
    }
    
    # Handle single-variable VIF (when vif returns a vector)
    if (is.vector(vif_values) && !is.null(names(vif_values))) {
      gvif <- vif_values
    } else if (is.matrix(vif_values) && ncol(vif_values) >= 2) {
      gvif <- vif_values[, 1]^(1 / (2 * vif_values[, 2]))
    } else {
      warning("Unexpected VIF output format, cannot extract GVIF.")
      break
    }
    
    high_vif_vars <- names(gvif[gvif > threshold])
    
    if (length(high_vif_vars) == 0) break  # If no high VIF variables, break loop
    
    # Remove only the variable with the highest VIF at each iteration (more robust approach)
    highest_vif_var <- names(gvif[which.max(gvif)])
    if (highest_vif_var %in% all.vars(full_formula)) { # Ensure variable is in the formula
      full_formula <- update(full_formula, paste(". ~ . -", highest_vif_var))
    } else { # If the highest VIF variable is not in the formula, re-evaluation may be needed
      break
    }
  }
  return(full_formula)
}

# Process bootstrap data, calculate GVIF and remove high collinearity variables
process_bootstrap_data <- function(bootstrap_data, threshold = 2.5) {
  results <- lapply(bootstrap_data, function(data) {
    formula <- tryCatch(
      filter_high_vif(data, threshold),
      error = function(e) {
        message("Error in filter_high_vif for a bootstrap sample: ", e$message)
        NULL
      }
    )
    list(data = data, formula = formula)
  })
  
  # Filter out data and formulas where formula is NULL
  filtered_results <- results[!sapply(results, function(x) is.null(x$formula))]
  
  list(
    data = lapply(filtered_results, function(x) x$data),
    formulas = lapply(filtered_results, function(x) x$formula)
  )
}

# Process bootstrap data with the updated function
message("Starting processing bootstrap data and VIF filtering...")
filtered_results_a <- process_bootstrap_data(bootstrap_data_a)
filtered_results_v <- process_bootstrap_data(bootstrap_data_v)
filtered_results_b <- process_bootstrap_data(bootstrap_data_b)
filtered_results_p <- process_bootstrap_data(bootstrap_data_p)
message("Bootstrap data and VIF filtering completed.")

# Get filtered data and formulas separately
filtered_data_a <- filtered_results_a$data
filtered_formulas_a <- filtered_results_a$formulas

filtered_data_v <- filtered_results_v$data
filtered_formulas_v <- filtered_results_v$formulas

filtered_data_b <- filtered_results_b$data
filtered_formulas_b <- filtered_results_b$formulas

filtered_data_p <- filtered_results_p$data
filtered_formulas_p <- filtered_results_p$formulas


# Function to fit the full model
fit_full_model <- function(data, formula) {
  tryCatch({
    # Perform model fitting
    fit <- rma.mv(
      yi, vi,
      mods = formula,
      random = ~1 | Study/Observation,
      method = "REML",
      data = data,
      control = control, # Use globally defined control object
      verbose = FALSE    # Can be TRUE for debugging
    )
    return(fit)  # Return the fitted model
  }, error = function(e) {
    # Gracefully handle errors
    message("Error in full model fitting: ", e$message)
    return(NULL)  # Return NULL on failure
  })
}

# Function to fit the null model
fit_null_model <- function(data) {
  tryCatch({
    fit <- rma.mv(
      yi, vi,
      mods = ~ -1, # Null model only fits intercept
      random = ~1 | Study/Observation,
      method = "REML",
      data = data,
      control = control, # Use globally defined control object
      verbose = FALSE
    )
    return(fit)
  }, error = function(e) {
    message("Error in null model fitting: ", e$message)
    return(NULL)
  })
}

# Progress printing function with timestamp
progress_message <- function(i, total, model_type_label, model_status = "") {
  ts <- format(Sys.time(), "%H:%M:%S")
  message(paste0("[", ts, "] Fitting ", model_type_label, " model ", i, " of ", total, " ", model_status))
}

# Iterate over data and formula lists using lapply, fit models separately, and store results separately
fit_list <- function(data_list, formula_list) {
  full_fits <- list()
  null_fits <- list()
  pseudo_r2s <- list()
  total_models <- length(data_list)
  
  for (i in seq_along(data_list)) {
    # Print full model fitting progress
    progress_message(i, total_models, "full")
    full_fit <- fit_full_model(data_list[[i]], formula_list[[i]])
    
    # Print null model fitting progress
    null_fit <- NULL # Set to NULL by default
    if (!is.null(full_fit)) { # Only try to fit null model if full model succeeds
      progress_message(i, total_models, "null")
      null_fit <- fit_null_model(data_list[[i]])
    }
    
    # Check if models fitted successfully
    if (!is.null(full_fit) && !is.null(null_fit)) {
      pseudo_r2 <- calculate_pseudoR2(full_fit, null_fit)
      full_fits[[i]] <- full_fit
      null_fits[[i]] <- null_fit
      pseudo_r2s[[i]] <- pseudo_r2
      progress_message(i, total_models, "both", "-> Success")
    } else {
      # Add more detailed error message
      message(sprintf("[WARN] Model fitting failed for element %d: Full fit: %s, Null fit: %s",
                      i,
                      ifelse(is.null(full_fit), "NULL", "Success"),
                      ifelse(is.null(null_fit), "NULL", "Success")))
      full_fits[[i]] <- NULL
      null_fits[[i]] <- NULL
      pseudo_r2s[[i]] <- NULL
    }
  }
  return(list(full_fits = full_fits, null_fits = null_fits, pseudo_r2s = pseudo_r2s))
}

# --- Execute model fitting ---
message("Starting fitting 'a' type models...")
fit_results_a <- fit_list(filtered_data_a, filtered_formulas_a)
message("Finished fitting 'a' type models.")

message("Starting fitting 'v' type models...")
fit_results_v <- fit_list(filtered_data_v, filtered_formulas_v)
message("Finished fitting 'v' type models.")

message("Starting fitting 'b' type models...")
fit_results_b <- fit_list(filtered_data_b, filtered_formulas_b)
message("Finished fitting 'b' type models.")

message("Starting fitting 'p' type models...")
fit_results_p <- fit_list(filtered_data_p, filtered_formulas_p)
message("Finished fitting 'p' type models.")

# Save RData image file of model results
save.image(file.path(output_path, paste0("fit_results_", timestamp, ".RData")))
message("Model fitting results saved to: ", file.path(output_path, paste0("fit_results_", timestamp, ".RData")))