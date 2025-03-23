library(metafor)
library(reshape2)
library(car)
library(openxlsx)
library(tidyverse)
library(boot)
library(aod)
library(gridExtra)

setwd("/data/home/wuhongyan/meta")
base_path <- "/data/home/wuhongyan/meta"
input_file <- file.path(base_path, "data/Supplementary Table 2.xlsx")
output_path <- file.path(base_path, "results")
set.seed(123)

# Load data-------
df <- read.xlsx("data/Supplementary Table 2.xlsx", 1)

# Remove unwanted columns using select()
dat <- df[, setdiff(names(df), c("Chinese_rodent_name", "Original_rodents",
                                 "Original_pathogens", "Standardization"))]
 
# Filter rows where pathogen_species is greater than or equal to 5
data <- dat %>%
  group_by(Pathogen_species) %>%
  filter(sum(Events) >= 5) %>%
  ungroup() 

# Keep only healthy rodents
data <- data %>%
  filter(Health_status == "Healthy") 

# Remove missing values ​​and unclear classification
data <- data[-which(is.na(data$NDVI)),,]
data <- data[-which(is.na(data$Time_interval)),,]
data <- data[-which(data$Rodents_type == 'Unknown'),,]
data <- data[-which(data$Testing_method =='Unidentified examination method'),,]
data <- data[-which(data$Sample_source%in%c('Unclassified sample')==TRUE),,]
nrow(data) #1214

# Define a function to set factor levels
define_factor_levels <- function(data, column_name, levels_order) {
  data[[column_name]] <- factor(data[[column_name]], levels = levels_order)
  return(data)
}

# Define the levels for each column
levels_list <- list(
  Rodents_type = c('Captive', 'Commensal', 'Wild'),
  Sample_source = c('Alimentary or intestinal sample', 'Blood or serum sample', 
                    'Body surface sample', 'Brain sample', 'Heart sample', 'Kidney sample', 
                    'Liver sample', 'Lung or respiratory tract sample', 'Spleen sample', 
                    'Mixed tissue sample'),
  Testing_method = c('Gross or microscopy examination', 'Molecular examination', 'Serological examination'),
  Rodent_family = c('Caviidae', 'Chinchillidae', 'Cricetidae', 'Dipodidae', 'Hystricidae', 
                    'Muridae', 'Myocastoridae', 'Nesomyidae', 'Sciuridae', 'Spalacidae'),
  Region = c('Central China', 'Eastern China', 'Northern China', 'Northeastern China', 
                            'Northwestern China', 'Southern China', 'Southwestern China'),
  Time_interval = c('1950-2000', '2001-2010', '2011-2019', '2020-2023'),
  NDVI = c('Peak NDVI', 'Higher NDVI', 'Medium NDVI', 'Low NDVI')
)

# Define a function to filter data by minimum group size
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
  data_a = data,
  data_v = data %>% filter(Pathogen_type == "Viruses"),
  data_b = data %>% filter(Pathogen_type == "Bacteria"),
  data_p = data %>% filter(Pathogen_type == "Parasites")
)
columns_to_check <- c("Rodents_type", "Sample_source", "Testing_method", 
                      "Rodent_family", "Region", "Time_interval", "NDVI")

# Apply transformations to each dataset
datasets <- lapply(datasets, function(dataset) {
  for (col_name in names(levels_list)) {
    dataset <- define_factor_levels(dataset, col_name, levels_list[[col_name]])
  }
  dataset$Observation <- factor(1:nrow(dataset))
  dataset$Study <- factor(dataset$Study)
  filter_data(dataset, columns_to_check)
})


## I² represents the percentage of variance in effect sizes that cannot be attributed to sampling error.
i2=function(model){
  
  ## metafor site code for I2
  W=diag(1/model$vi)  
  X=model.matrix(model)
  P=W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2=100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  I2=round(I2,2)
  
  ## summarize by each variance component
  allI2=100 * model$sigma2 / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P))) 
  return(list(I2=I2,allI2=allI2))
}

#  Calculating effect size
calc_effect_size <- function(data) {
  data <- data.frame(data, escalc(xi = data$Events, ni = data$Total, measure = "PLN"))
  data$original_ratio <- exp(data$yi)
  return(data)
}

datasets <- lapply(datasets, calc_effect_size)
# Assign back the modified datasets
data_a <- datasets[[1]]
data_v <- datasets[[2]]
data_b <- datasets[[3]]
data_p <- datasets[[4]]


###############  Full dataset; full_model
control <- list(optimizer="optim", optmethod="L-BFGS-B")

###############  Sub_dataset; Sub_models
boot.func <- function(dataset) {
  indices <- sample(1:nrow(dataset), replace = TRUE)  
  dataset[indices, ]  
}

# Bootstrap resampling
set.seed(123)
bootstrap_data_a <- lapply(1:999, function(i) boot.func(data_a))
bootstrap_data_v <- lapply(1:999, function(i) boot.func(data_v))
bootstrap_data_b <- lapply(1:999, function(i) boot.func(data_b))
bootstrap_data_p <- lapply(1:999, function(i) boot.func(data_p))

# Update data_list to contain multiple resampled data lists
resampled_data_list <- list(a = bootstrap_data_a, v = bootstrap_data_v, b = bootstrap_data_b, p = bootstrap_data_p)

# Define function: Calculate GVIF and remove highly collinear variables
filter_high_vif <- function(data, threshold = 2) {
  full_formula <- yi ~ Sample_source + Region + Rodents_type + 
    Testing_method + Rodent_family + Time_interval + NDVI - 1
  repeat {
    lm_model <- lm(full_formula, data = data)
    vif_values <- vif(lm_model)
    gvif <- vif_values[, 1]^(1 / (2 * vif_values[, 2]))  
    high_vif_vars <- names(gvif[gvif > threshold])  
    if (length(high_vif_vars) == 0) break  
    full_formula <- update(full_formula, paste(". ~ . -", high_vif_vars[1]))  
  }
  return(full_formula)  
}

# Remove highly collinear variables
process_bootstrap_data <- function(bootstrap_data, threshold = 2) {
  results <- lapply(bootstrap_data, function(data) {
    formula <- tryCatch(
      filter_high_vif(data, threshold),
      error = function(e) NULL  
    )
    list(data = data, formula = formula)  
  })
  
  # Filter out data and formulas whose formulas are not NULL
  filtered_results <- results[!sapply(results, function(x) is.null(x$formula))]
  
  # Returns a list of filtered data and formulas
  list(
    data = lapply(filtered_results, function(x) x$data),
    formulas = lapply(filtered_results, function(x) x$formula)
  )
}

# Using the updated function to process bootstrapped data
set.seed(123)
filtered_results_a <- process_bootstrap_data(bootstrap_data_a)
filtered_results_v <- process_bootstrap_data(bootstrap_data_v)
filtered_results_b <- process_bootstrap_data(bootstrap_data_b)
filtered_results_p <- process_bootstrap_data(bootstrap_data_p)

# Get the filtered data and formulas separately
filtered_data_a <- filtered_results_a$data
filtered_formulas_a <- filtered_results_a$formulas

filtered_data_v <- filtered_results_v$data
filtered_formulas_v <- filtered_results_v$formulas

filtered_data_b <- filtered_results_b$data
filtered_formulas_b <- filtered_results_b$formulas

filtered_data_p <- filtered_results_p$data
filtered_formulas_p <- filtered_results_p$formulas


# Define the function to run the analysis and capture results
run_analysis_results <- function(data, formula) {
  tryCatch({
    # Execute model fitting
    fit <- rma.mv(
      yi, vi,
      mods = as.formula(formula),
      random = ~1 | Study/Observation,
      method = "REML",
      data = data,
      control = control,
      verbose = TRUE
    )
    return(fit)  # Return the fitted model
  }, error = function(e) {
    # Handle errors gracefully
    message("Error in model fitting: ", e$message)
    return(NULL)  # Return NULL to indicate failure
  })
}

# Iterating model evaluation
fit_a <- mapply(run_analysis_results, filtered_data_a, filtered_formulas_a, SIMPLIFY = FALSE)
fit_v <- mapply(run_analysis_results, filtered_data_v, filtered_formulas_v, SIMPLIFY = FALSE)
fit_b <- mapply(run_analysis_results, filtered_data_b, filtered_formulas_b, SIMPLIFY = FALSE)
fit_p <- mapply(run_analysis_results, filtered_data_p, filtered_formulas_p, SIMPLIFY = FALSE)

# Filter out NULL objects in fit_results_
fit_results_a <- Filter(function(x) !is.null(x), fit_a)
fit_results_v <- Filter(function(x) !is.null(x), fit_v)
fit_results_b <- Filter(function(x) !is.null(x), fit_b)
fit_results_p <- Filter(function(x) !is.null(x), fit_p)

# Function to extract results from the model fit
extract_model_results <- function(fit) {
  if (!is.null(fit)) {
    model_results <- coef(summary(fit))  # Extract model coefficients
    return(data.frame(
      Term = rownames(model_results),
      Estimate = model_results[, "estimate"],
      StdError = model_results[, "se"],
      ZValue = model_results[, "zval"],
      PValue = model_results[, "pval"],
      CILower = model_results[, "ci.lb"],
      CIUpper = model_results[, "ci.ub"]
    ))
  } else {
    return(NULL)  # Return NULL if fit is empty
  }
}

# Extract and store model results for each dataset
results_a <- lapply(fit_results_a, extract_model_results)
results_v <- lapply(fit_results_v, extract_model_results)
results_b <- lapply(fit_results_b, extract_model_results)
results_p <- lapply(fit_results_p, extract_model_results)

# Get unique terms across all results
all_terms_a <- unique(unlist(lapply(results_a, function(x) if (!is.null(x)) x$Term)))
all_terms_v <- unique(unlist(lapply(results_v, function(x) if (!is.null(x)) x$Term)))
all_terms_b <- unique(unlist(lapply(results_b, function(x) if (!is.null(x)) x$Term)))
all_terms_p <- unique(unlist(lapply(results_p, function(x) if (!is.null(x)) x$Term)))

# Unified results for each type of results (a, v, b, p)
unified_results_a <- do.call(rbind, lapply(results_a, function(result) {
  if (is.null(result)) {
    return(data.frame(
      Term = all_terms_a,
      Estimate = NA, StdError = NA, ZValue = NA,
      PValue = NA, CILower = NA, CIUpper = NA
    ))
  }
  merged_result <- merge(data.frame(Term = all_terms_a), result, by = "Term", all.x = TRUE)
  return(merged_result)
}))

unified_results_v <- do.call(rbind, lapply(results_v, function(result) {
  if (is.null(result)) {
    return(data.frame(
      Term = all_terms_v,
      Estimate = NA, StdError = NA, ZValue = NA,
      PValue = NA, CILower = NA, CIUpper = NA
    ))
  }
  merged_result <- merge(data.frame(Term = all_terms_v), result, by = "Term", all.x = TRUE)
  return(merged_result)
}))

unified_results_b <- do.call(rbind, lapply(results_b, function(result) {
  if (is.null(result)) {
    return(data.frame(
      Term = all_terms_b,
      Estimate = NA, StdError = NA, ZValue = NA,
      PValue = NA, CILower = NA, CIUpper = NA
    ))
  }
  merged_result <- merge(data.frame(Term = all_terms_b), result, by = "Term", all.x = TRUE)
  return(merged_result)
}))

unified_results_p <- do.call(rbind, lapply(results_p, function(result) {
  if (is.null(result)) {
    return(data.frame(
      Term = all_terms_p,
      Estimate = NA, StdError = NA, ZValue = NA,
      PValue = NA, CILower = NA, CIUpper = NA
    ))
  }
  merged_result <- merge(data.frame(Term = all_terms_p), result, by = "Term", all.x = TRUE)
  return(merged_result)
}))

output_dir <- "model_results"  

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.xlsx(unified_results_a, file = file.path(output_dir, "unified_results_a.xlsx"), rowNames = FALSE)
write.xlsx(unified_results_v, file = file.path(output_dir, "unified_results_v.xlsx"), rowNames = FALSE)
write.xlsx(unified_results_b, file = file.path(output_dir, "unified_results_b.xlsx"), rowNames = FALSE)
write.xlsx(unified_results_p, file = file.path(output_dir, "unified_results_p.xlsx"), rowNames = FALSE)

# Compute mean results (delete NAs)
mean_results_a <- aggregate(
  . ~ Term,  
  data = na.omit(unified_results_a),  
  FUN = mean
)

mean_results_v <- aggregate(
  . ~ Term,  
  data = na.omit(unified_results_v),  
  FUN = mean
)

mean_results_b <- aggregate(
  . ~ Term,  
  data = na.omit(unified_results_b),  
  FUN = mean
)

mean_results_p <- aggregate(
  . ~ Term,  
  data = na.omit(unified_results_p),  
  FUN = mean
)

output_dir <- "model_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE)
}

# Save the mean results 
output_file_a <- file.path(output_dir, paste0("mean_results_a_", paste(sample(letters, 5, TRUE), collapse = ""), ".xlsx"))
tryCatch({
  write.xlsx(mean_results_a, output_file_a)
  message("Results successfully saved to ", output_file_a)
}, error = function(e) {
  message("Error in saving results for V: ", e$message)
})

# Save the mean results for unified_results_v
output_file_v <- file.path(output_dir, paste0("mean_results_v_", paste(sample(letters, 5, TRUE), collapse = ""), ".xlsx"))
tryCatch({
  write.xlsx(mean_results_v, output_file_v)
  message("Results successfully saved to ", output_file_v)
}, error = function(e) {
  message("Error in saving results for V: ", e$message)
})

# Save the mean results for unified_results_b
output_file_b <- file.path(output_dir, paste0("mean_results_b_", paste(sample(letters, 5, TRUE), collapse = ""), ".xlsx"))
tryCatch({
  write.xlsx(mean_results_b, output_file_b)
  message("Results successfully saved to ", output_file_b)
}, error = function(e) {
  message("Error in saving results for B: ", e$message)
})

# Save the mean results for unified_results_p
output_file_p <- file.path(output_dir, paste0("mean_results_p_", paste(sample(letters, 5, TRUE), collapse = ""), ".xlsx"))
tryCatch({
  write.xlsx(mean_results_p, output_file_p)
  message("Results successfully saved to ", output_file_p)
}, error = function(e) {
  message("Error in saving results for P: ", e$message)
})


# Save the model fitting and summary statistics (i2, QM)
i2_results_a <- do.call(rbind, lapply(fit_results_a, i2))
i2_results_v <- do.call(rbind, lapply(fit_results_v, i2))
i2_results_b <- do.call(rbind, lapply(fit_results_b, i2))
i2_results_p <- do.call(rbind, lapply(fit_results_p, i2))

# Create an Excel workbook to store the results
output_dir <- "model_results"
dir.create(output_dir, showWarnings = FALSE)  
output_file <- file.path(output_dir, paste0("i2_results_", paste(sample(letters, 5, TRUE), collapse = ""), ".xlsx"))
message("Output file path: ", output_file)


if (!file.access(output_dir, 2) == 0) stop("No write access to output directory.")
tryCatch({
  wb <- createWorkbook()
  add_non_empty_sheet <- function(wb, sheet_name, data) {
    if (!is.null(data) && nrow(data) > 0) {  
      addWorksheet(wb, sheet_name)          
      writeData(wb, sheet_name, data)      
    } else {
      message("Skipping ", sheet_name, " as it is empty or NULL.")
    }
  }
  
  add_non_empty_sheet(wb, "i2_a", i2_results_a)
  add_non_empty_sheet(wb, "i2_v", i2_results_v)
  add_non_empty_sheet(wb, "i2_b", i2_results_b)
  add_non_empty_sheet(wb, "i2_p", i2_results_p)
  
  saveWorkbook(wb, output_file, overwrite = TRUE)
  message("i2 results successfully saved to ", output_file)
}, error = function(e) {
  message("Error in saving i2 results to Excel: ", e$message)
})
