# Load necessary libraries
library(terra)
library(sf)
library(raster)

# File paths for raster data
Climate <- list.files(r"(BRT/DATA/ProjectPost/Climate)", 
                      pattern = ".tif$", full.names = TRUE, ignore.case = TRUE)
Elevation <- list.files(r"(BRT/DATA/RasterRaw/Elevation)", 
                        pattern = ".tif$", full.names = TRUE, ignore.case = TRUE)
GDP1995 <- list.files(r"(BRT/DATA/RasterRaw/GDP/1995)", 
                      pattern = ".tif$", full.names = TRUE, ignore.case = TRUE)
GDP2019 <- list.files(r"(BRT/DATA/RasterRaw/GDP/2019)", 
                      pattern = ".tif$", full.names = TRUE, ignore.case = TRUE)
PD1995 <- list.files(r"(BRT/DATA/RasterRaw/PD/1995)", 
                     pattern = ".tif$", full.names = TRUE, ignore.case = TRUE)
PD2019 <- list.files(r"(BRT/DATA/RasterRaw/PD/2019)", 
                     pattern = ".tif$", full.names = TRUE, ignore.case = TRUE)
LandCover_full <- list.files(r"(BRT/DATA/RasterRaw/Landcover/full_version)", 
                             pattern = ".tif$", full.names = TRUE, ignore.case = TRUE)
LandCover_reduced <- list.files(r"(BRT/DATA/RasterRaw/Landcover/reduced_version)", 
                                pattern = ".tif$", full.names = TRUE, ignore.case = TRUE)

# Read specific raster files
MammalRichness <- rast("BRT/DATA/RasterRaw/MammalRichness/MammalRichness.tif")
ReportEffort <- rast("BRT/DATA/RasterRaw/ReportEffort/ReportEffort.tif")
RodentsPoints <- rast("BRT/DATA/RasterRaw/RodentsPoints/RodentsPoints.tif")

# Convert file paths to raster stacks
Climate <- terra::rast(Climate)
Elevation <- terra::rast(Elevation)
GDP1995 <- terra::rast(GDP1995)
GDP2019 <- terra::rast(GDP2019)
PD1995 <- terra::rast(PD1995)
PD2019 <- terra::rast(PD2019)
LandCover_full <- terra::rast(LandCover_full)
LandCover_reduced <- terra::rast(LandCover_reduced)

# Load shapefile
china.shp <- st_read("BRT/DATA/ChinaSHP/China.shp", quiet = FALSE)
china.shp <- st_transform(china.shp, crs = 4326)

# Define a helper function to crop and average raster layers
crop_and_average <- function(rasters, extent_obj) {
  cropped <- lapply(rasters, crop, extent(extent_obj))
  mean_raster <- Reduce("+", cropped) / length(cropped)
  return(mean_raster)
}

# Crop and average climate data
EVP <- crop_and_average(list(Climate$EVP1960_ProjectRaster1, Climate$EVP1980_ProjectRaster11, 
                             Climate$EVP2000_ProjectRaster1, Climate$EVP2020_ProjectRaster1), china.shp)
GST <- crop_and_average(list(Climate$GST1960_ProjectRaster1, Climate$GST1980_ProjectRaster1, 
                             Climate$GST2000_ProjectRaster1, Climate$GST2020_ProjectRaster1), china.shp)
PRE <- crop_and_average(list(Climate$PRE1960_ProjectRaster1, Climate$PRE1980_ProjectRaster1, 
                             Climate$PRE2000_ProjectRaster1, Climate$PRE2020_ProjectRaster1), china.shp)
PRS <- crop_and_average(list(Climate$PRS1960_ProjectRaster1, Climate$PRS1980_ProjectRaster1, 
                             Climate$PRS2000_ProjectRaster1, Climate$PRS2020_ProjectRaster1), china.shp)
RHU <- crop_and_average(list(Climate$RHU1960_ProjectRaster1, Climate$RHU1980_ProjectRaster1, 
                             Climate$RHU2000_ProjectRaster1, Climate$RHU2020_ProjectRaster1), china.shp)
SSD <- crop_and_average(list(Climate$SSD1960_ProjectRaster1, Climate$SSD1980_ProjectRaster1, 
                             Climate$SSD2000_ProjectRaster1, Climate$SSD2020_ProjectRaster1), china.shp)
TEM <- crop_and_average(list(Climate$TEM1960_ProjectRaster1, Climate$TEM1980_ProjectRaster1, 
                             Climate$TEM2000_ProjectRaster1, Climate$TEM2020_ProjectRaster1), china.shp)
WIN <- crop_and_average(list(Climate$WIN1960_ProjectRaster1, Climate$WIN1980_ProjectRaster1, 
                             Climate$WIN2000_ProjectRaster1, Climate$WIN2020_ProjectRaster1), china.shp)

# Create mask
Mask <- EVP
values(Mask) <- NA

# Resample rasters to match Mask
Elevation <- resample(Elevation, Mask)
GDP1995 <- resample(GDP1995, Mask)
GDP2019 <- resample(GDP2019, Mask)
PD1995 <- resample(PD1995, Mask)
PD2019 <- resample(PD2019, Mask)

# Crop and resample land cover data
LandCover_types <- c("Barren", "CMV", "DBT", "DNT", "EBT", "HV", "MT", "OW", "RFV", "Shrubs", "Snow", "Urban")
LandCover_full_resampled <- list()
LandCover_reduced_resampled <- list()
for (type in LandCover_types) {
  full_cropped <- crop(LandCover_full[[type]], extent(china.shp))
  full_resampled <- resample(full_cropped, Mask)
  LandCover_full_resampled[[type]] <- full_resampled
  
  reduced_cropped <- crop(LandCover_reduced[[type]], extent(china.shp))
  reduced_resampled <- resample(reduced_cropped, Mask)
  LandCover_reduced_resampled[[type]] <- reduced_resampled
}

# Resample MammalRichness and ReportEffort
MammalRichness <- resample(crop(MammalRichness, extent(china.shp)), Mask)
ReportEffort <- resample(ReportEffort, Mask)
RodentsPoints <- resample(RodentsPoints, Mask)

# Save rasters to disk
output_dir <- "BRT/DATA/ResamplePost/Model/TrainRodents/"
output_dir_predictors <- "BRT/DATA/ResamplePost/Model/PredictorsRodents/"

# Write full land cover raster data
for (name in names(LandCover_full_resampled)) {
  writeRaster(LandCover_full_resampled[[name]], filename = paste0(output_dir, name, "_full.tif"))
}

# Write reduced land cover raster data
for (name in names(LandCover_reduced_resampled)) {
  writeRaster(LandCover_reduced_resampled[[name]], filename = paste0(output_dir_predictors, name, "_reduced.tif"))
}

# Write climate raster data
climate_rasters <- list(EVP = EVP, GST = GST, PRE = PRE, PRS = PRS, RHU = RHU, SSD = SSD, TEM = TEM, WIN = WIN)
for (name in names(climate_rasters)) {
  writeRaster(climate_rasters[[name]], filename = paste0(output_dir, name, ".tif"))
}

# Write additional raster data
additional_rasters <- list(Elevation = Elevation, GDP1995 = GDP1995, GDP2019 = GDP2019, PD1995 = PD1995, PD2019 = PD2019, MammalRichness = MammalRichness, ReportEffort = ReportEffort, RodentsPoints = RodentsPoints)
for (name in names(additional_rasters)) {
  writeRaster(additional_rasters[[name]], filename = paste0(output_dir, name, ".tif"))
}


