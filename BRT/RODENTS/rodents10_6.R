# Load packages
library(dplyr)
library(sf) 
library(blockCV) 
library(terra) 
library(dismo) 
library(gbm) 
library(precrec) 

# Read raster data
tif_file_name <- list.files("BRT/DATA/ResamplePost/Model/TrainRodents", 
                            pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)

predictors_file_name <- list.files("BRT/DATA/ResamplePost/Model/PredictorsRodents", 
                                   pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)

rasters <- terra::rast(tif_file_name)
predictors <- terra::rast(predictors_file_name)

points <- rast("BRT/DATA/ResamplePost/Model/PredictorsCases/RodentsPoints.tif")
points <- as.data.frame(points, xy=TRUE, na.rm=FALSE)
points[is.na(points)] <- 0
points <- rename(points, "Abspres" = "RodentsPoints")
rodentspoints <- sf::st_as_sf(points, coords = c("x", "y"), crs = CRS("+init=epsg:4326"))


# Extract point information and create test and training sets
model_data <- terra::extract(rasters, rodentspoints, df = TRUE, ID = FALSE, xy=TRUE)
model_data$Abspres <- as.numeric(rodentspoints$Abspres)

# Delete points that are all NA
miss_counts <- rowSums(is.na(model_data))
model_data$counts <- miss_counts
model_data <- subset(model_data, model_data$counts < 3)

rodents_table <- data.frame(model_data$x, model_data$y, model_data$Abspres)
rodents_table <- rename(rodents_table, "x" = "model_data.x")
rodents_table <- rename(rodents_table, "y" = "model_data.y")
rodents_table <- rename(rodents_table, "Abspres" = "model_data.Abspres")
rodents_table <- sf::st_as_sf(rodents_table, coords = c("x", "y"), crs = 4326)

# Create a table to store predictions
test_table <- rodents_table
test_table$preds <- NA
model_data <- model_data[, c(6:9, 11, 14, 15, 17, 19, 22, 24, 28)]
colnames(model_data)

# Create cross-validation blocks after removing points with missing values
set.seed(1242)
scv1 <- cv_spatial(
  x = rodents_table,
  column = "Abspres", # The response column (binary or multi-class)
  k = 10, # Number of folds
  size = 800000, # Size of the blocks in meters
  selection = "random", # Random blocks-to-fold assignment
  iteration = 50, # Find evenly dispersed folds
  progress = FALSE, # Turn off progress bar
  biomod2 = TRUE, # Also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # Options from cv_plot for better color contrast
) 

# Model building and AUC calculation
model_list <- list()
folds <- scv1$folds_list
for (k in seq_len(length(folds))) {
  trainSet <- unlist(folds[[k]][1]) # Training set indices; first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices; second element
  model_list[[k]] <- gbm.step(model_data[trainSet, ], gbm.x = 1:11, gbm.y = 12, family = "bernoulli", tree.complexity = 5,
                              learning.rate = 0.001, bag.fraction = 0.5, max.trees = 20000) 
  test_table$preds[testSet] <- predict.gbm(model_list[[k]], model_data[testSet, ], n.trees = model_list[[k]]$gbm.call$best.trees, type = "response") 
}
precrec_obj <- evalmod(scores = test_table$preds, labels = test_table$Abspres)
auc(precrec_obj)
write.csv(test_table, "BRT/OUTPUT/RODENTS/TEST/test6.csv")

# Write the average tif
rodents.risk <- rast(predictors)  # Use the same dimensions (number of rows, columns, and layers) as predictors
rodents.risk[] <- 0  # Initialize all values to 0
for (i in seq_along(model_list)) {
  risk <- predict(predictors, model_list[[i]], n.trees = model_list[[i]]$gbm.call$best.trees, type = "response")
  rodents.risk <- rodents.risk + risk # Accumulate the prediction results of the current model (risk) to the previous prediction results (rodents.risk)
}
rodents.risk <- rodents.risk / length(model_list)
writeRaster(rodents.risk, 
            filename = "BRT/OUTPUT/RODENTS/TIF/rodents_risk6.tif", overwrite = TRUE)

# Relative contribution
contribution.all <- data.frame()
for (i in seq_along(model_list)) {
  # Extract variable names and relative contribution information, and create a new dataframe a
  a <- data.frame(var = model_list[[i]]$contributions$var, contribution = model_list[[i]]$contributions$rel.inf)
  # Sort the dataframe a by the values in the var column
  a <- a[order(a$var),]
  # Rename the columns of the dataframe a to prevent duplication
  names(a)[2] <- paste("contribution", i, sep = "_")
  # If contribution.all is empty, assign a to it directly; otherwise, merge columns
  if (i == 1) {
    contribution.all <- a
  } else {
    contribution.all <- merge(contribution.all, a, by = "var", all = TRUE)
  }
}

write.csv(contribution.all, "BRT/OUTPUT/RODENTS/CONTRIBUTION/contribution6.csv")
save.image("BRT/OUTPUT/RODENTS/IMAGE/image6.RData")
