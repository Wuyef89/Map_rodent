library(dplyr)
library(sf) 
library(blockCV) 
library(terra) 
library(dismo) 
library(gbm) 
library(precrec)

setwd("/data/home/wuhongyan/sourcedata/")
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

# Assuming province_data is a list of data frames and data is a data frame
# Set seed for reproducibility
set.seed(0021)
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
  if (matching_row$typhus == 0) {
    province_df$Abspres <- 0
  } else if (nrow(province_df) < matching_row$typhus) {
    province_df$Abspres <- 1
  } else {
    num_positive <- matching_row$typhus
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
set.seed(0022)

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
model_data<-model_data[,c(6,8,11,21,15,14,17,9,5,24,19,25)] # Variables selected based on 5-fold cross-validation
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
  model_list[[k]] <- gbm.step(model_data[trainSet, ], gbm.x=1:11,gbm.y=12,family="bernoulli",
                              tree.complexity=5,learning.rate=0.001,bag.fraction=0.5,max.trees=20000) 
  test_table$preds[testSet] <- predict.gbm(model_list[[k]], model_data[testSet, ], 
                                           n.trees=model_list[[k]]$gbm.call$best.trees, type="response")
}

precrec_obj <- evalmod(scores = test_table$preds, labels = test_table$Abspres)
auc(precrec_obj)

write.csv(test_table,"BRT/OUTPUT/TYPHUS/TEST/test8.csv")

#Write the average tif
typhus.risk <- rast(predictors)  
typhus.risk[] <- 0  
for (i in seq_along(model_list)) {
  risk <- predict(predictors, model_list[[i]], n.trees=model_list[[i]]$gbm.call$best.trees, type="response")
  typhus.risk <- typhus.risk + risk 
}
typhus.risk <- typhus.risk / length(model_list)
writeRaster(typhus.risk, 
            filename = "BRT/OUTPUT/TYPHUS/TIF/typhus_risk8.tif", overwrite=TRUE)

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

write.csv(contribution.all, "BRT/OUTPUT/TYPHUS/CONTRIBUTION/contribution8.csv")
save.image("BRT/OUTPUT/TYPHUS/IMAGE/image8.RData")

##############   Response curve of predictor (marginal effect)   #########################################

getVar = function(gbm.object, predictor_of_interest) {
  gbm.call <- gbm.object$gbm.call  
  gbm.x <- gbm.call$gbm.x  
  pred.names <- gbm.call$predictor.names  
  response.name <- gbm.call$response.name  
  data <- gbm.call$dataframe  
  k <- match(predictor_of_interest, pred.names)  
  var.name <- gbm.call$predictor.names[k]  
  pred.data <- data[, gbm.call$gbm.x[k]]  
  response.matrix <- gbm::plot.gbm(gbm.object, k, return.grid = TRUE)  
  data.frame(predictors = response.matrix[, 1],  
             responses = response.matrix[, 2] - mean(response.matrix[, 2]) 
  )
}

#1
Elevation8 <- lapply(1:length(model_list), function(i){
  data.frame(getVar(model_list[[i]],"Elevation"), model = i)
})
Elevation8 <- do.call(cbind, Elevation8)

#2
GDP8 <- lapply(1:length(model_list),function(i){
  data.frame(getVar(model_list[[i]],"GDP"), model = i)
})
GDP8 <- do.call(rbind, GDP8)

#3
MammalRichness8 <- lapply(1:length(model_list),function(i){
  data.frame(getVar(model_list[[i]],"MammalRichness"), model = i)
})
MammalRichness8 <- do.call(cbind, MammalRichness8)

#4
SSD8 <- lapply(1:length(model_list),function(i){
  data.frame(getVar(model_list[[i]],"SSD"), model = i)
})
SSD8 <- do.call(cbind, SSD8)

#5
PRS8 <- lapply(1:length(model_list),function(i){
  data.frame(getVar(model_list[[i]],"PRS"), model = i)
})
PRS8 <- do.call(cbind, PRS8)

#6
PRE8 <- lapply(1:length(model_list),function(i){
  data.frame(getVar(model_list[[i]],"PRE"), model = i)
})
PRE8 <- do.call(cbind, PRE8)

#7
RHU8 <- lapply(1:length(model_list),function(i){
  data.frame(getVar(model_list[[i]],"RHU"), model = i)
})
RHU8 <- do.call(cbind, RHU8)

#8
GST8 <- lapply(1:length(model_list),function(i){
  data.frame(getVar(model_list[[i]],"GST"), model = i)
})
GST8 <- do.call(cbind, GST8)

#9
EBT8 <- lapply(1:length(model_list),function(i){
  data.frame(getVar(model_list[[i]],"EBT"), model = i)
})
EBT8 <- do.call(cbind, EBT8)

#10
WIN8 <- lapply(1:length(model_list),function(i){
  data.frame(getVar(model_list[[i]],"WIN"), model = i)
})
WIN8 <- do.call(cbind, WIN8)

#11
Shrubs8 <- lapply(1:length(model_list),function(i){
  data.frame(getVar(model_list[[i]],"Shrubs"), model = i)
})
Shrubs8 <- do.call(cbind, Shrubs8)


write.csv(Elevation8,"BRT/OUTPUT/TYPHUS/VAR//Elevation8.csv")
write.csv(GDP8,"BRT/OUTPUT/TYPHUS/VAR//GDP8.csv")
write.csv(MammalRichness8,"BRT/OUTPUT/TYPHUS/VAR//MammalRichness8.csv")
write.csv(SSD8,"BRT/OUTPUT/TYPHUS/VAR//SSD8.csv")
write.csv(PRS8,"BRT/OUTPUT/TYPHUS/VAR//PRS8.csv")
write.csv(PRE8,"BRT/OUTPUT/TYPHUS/VAR//PRE8.csv")
write.csv(RHU8,"BRT/OUTPUT/TYPHUS/VAR//RHU8.csv")
write.csv(GST8,"BRT/OUTPUT/TYPHUS/VAR//GST8.csv")
write.csv(EBT8,"BRT/OUTPUT/TYPHUS/VAR//EBT8.csv")
write.csv(WIN8,"BRT/OUTPUT/TYPHUS/VAR//WIN8.csv")
write.csv(Shrubs8,"BRT/OUTPUT/TYPHUS/VAR//Shrubs8.csv")

