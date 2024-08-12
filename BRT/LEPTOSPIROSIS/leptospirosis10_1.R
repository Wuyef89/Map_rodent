library(dplyr)
library(sf) 
library(blockCV) 
library(terra) 
library(dismo) 
library(gbm) 
library(precrec)

# Function to process the shapefiles and PD data
process_province_data <- function(province_name, path_to_shp, path_to_pd_tif) {
  shp_path <- file.path(path_to_shp, paste0(province_name, ".shp"))
  pd_path <- file.path(path_to_pd_tif, "PD.tif")
  shp <- read_sf(shp_path)
  pd <- mask(rast(pd_path), shp)
  pd_df <- as.data.frame(pd, xy = TRUE) %>% 
    mutate(PD = PD / sum(PD, na.rm = TRUE)) %>% 
    na.omit() %>%
    mutate(id = rownames(.))
  pd_df$province <- province_name  # Add province name
  return(pd_df)
}

# Set the paths
path_to_shp <- "BRT/DATA/ProvincesSHP"
path_to_pd_tif <- "BRT/DATA/RasterRaw/PD/1995"

# Names of the provinces
provinces <- c("Beijing", "Tianjin", "Hebei", "Shanxi", "InnerMongolia", "Liaoning",
               "Jilin", "Heilongjiang", "Shanghai", "Jiangsu", "Zhejiang", "Anhui",
               "Fujian", "Jiangxi", "Shandong", "Henan", "Hubei", "Hunan", "Guangdong",
               "Guangxi", "Hainan", "Chongqing", "Sichuan", "Guizhou", "Yunnan",
               "Xizang", "Shaanxi", "Gansu", "Qinghai", "Ningxia", "Xinjiang")

# Process data for all provinces
province_data <- lapply(provinces, process_province_data, path_to_shp, path_to_pd_tif)

set.seed(0007)
data <- read.csv("BRT/DATA/case.csv")
# Create an empty list to store the results
final_province_data <- list()

# Loop through the provinces and perform the operation as described
for (i in seq_along(provinces)) {
  # Get the province name and corresponding data
  province_name <- provinces[i]
  province_df <- province_data[[i]]
  
  # Find the matching row in 'data'
  matching_row <- data[data$province == province_name, ]
  
  # Check the conditions and set Abspres accordingly
  if (matching_row$leptospirosis == 0) {
    province_df$Abspres <- 0
  } else if (nrow(province_df) < matching_row$leptospirosis) {
    province_df$Abspres <- 1
  } else {
    num_positive <- matching_row$leptospirosis
    positive_indices <- sample(province_df$id, num_positive, prob = province_df$PD, replace = FALSE)
    province_df$Abspres <- ifelse(province_df$id %in% positive_indices, 1, 0)
  }
  
  # Add the processed data frame to the results list
  final_province_data[[province_name]] <- province_df
}

# Combine all the data frames in the final_province_data list into one data frame
point <- do.call(rbind, final_province_data)

# Remove columns 3, 4, and 5
point <- point[,-c(3,4,5)]
point <- point[rowSums(is.na(point)) == 0, ] # Delete rows with missing values

tif_file_name <- list.files(r"(BRT/DATA/ResamplePost/Model/TrainCases)", 
                            pattern = ".tif$", full.names = TRUE, ignore.case = TRUE)
predictors_file_name <- list.files(r"(BRT/DATA/ResamplePost/Model/PredictorsCases)", 
                                   pattern = ".tif$", full.names = TRUE, ignore.case = TRUE)


rasters<-terra::rast(tif_file_name)
predictors<-terra::rast(predictors_file_name)

# Convert dataframe to sf object, specifying coordinates and coordinate reference system
point <- sf::st_as_sf(point, coords = c("x", "y"), crs = 4326)

# Create spatial validation blocks
set.seed(0008)

scv1 <- cv_spatial(
  x = point,  
  column = "Abspres", 
  k = 10, # number of folds 
  size = 800000, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  iteration = 50, # find evenly dispersed folds 
  progress = FALSE, 
  biomod2 = TRUE, 
  raster_colors = terrain.colors(10, rev = TRUE) 
) 

# Extract the values of the rasters (24 tif predictor variables) at each point to use as model data
model_data <- terra::extract(rasters, point, df = TRUE, ID = FALSE) 
model_data$Abspres <- as.numeric(point$Abspres)
model_data$RodentsPoints<-as.factor(model_data$RodentsPoints)
colnames(model_data)
model_data<-model_data[,c(8,17,14,15,21,11,10,1,25)] # Variables selected based on 5-fold cross-validation
colnames(model_data)

# Create a table to store the predictions
test_table <- point
test_table$preds <- NA  

# Build models and calculate AUC
model_list<-list()  
folds <- scv1$folds_list 
for(k in seq_len(length(folds))){  
  trainSet <- unlist(folds[[k]][1]) 
  testSet <- unlist(folds[[k]][2]) 
  model_list[[k]] <- gbm.step(model_data[trainSet, ], gbm.x=1:8,gbm.y=9,family="bernoulli",
                              tree.complexity=5,learning.rate=0.001,bag.fraction=0.5,max.trees=20000) 
  test_table$preds[testSet] <- predict.gbm(model_list[[k]], model_data[testSet, ], 
                                           n.trees=model_list[[k]]$gbm.call$best.trees, type="response")
}

precrec_obj <- evalmod(scores = test_table$preds, labels = test_table$Abspres)
auc(precrec_obj)

write.csv(test_table,"BRT/OUTPUT/LEPTOSPIROSIS/TEST/test1.csv")

#Write the average tif
leptospirosis.risk <- rast(predictors)  
leptospirosis.risk[] <- 0  
for (i in seq_along(model_list)) {
  risk <- predict(predictors, model_list[[i]], n.trees=model_list[[i]]$gbm.call$best.trees, type="response")
  leptospirosis.risk <- leptospirosis.risk + risk 
}
leptospirosis.risk <- leptospirosis.risk / length(model_list)
writeRaster(leptospirosis.risk, 
            filename = "BRT/OUTPUT/LEPTOSPIROSIS/TIF/leptospirosis_risk1.tif", overwrite=TRUE)

#Calculating Relative Importance(contribution) of predictors
contribution.all <- data.frame()
for (i in seq_along(model_list)) {
  a <- data.frame(var=model_list[[i]]$contributions$var, contribution=model_list[[i]]$contributions$rel.inf)
  a <- a[order(a$var),]
  names(a)[2] <- paste("contribution", i, sep="_")
  if (i == 1) {
    contribution.all <- a
  } else {
    contribution.all <- merge(contribution.all, a, by="var", all=TRUE)
  }
}

write.csv(contribution.all,"BRT/OUTPUT/LEPTOSPIROSIS/CONTRIBUTION/contribution1.csv")
save.image("BRT/OUTPUT/LEPTOSPIROSIS/IMAGE/image1.RData")
