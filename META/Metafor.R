
library(metafor)
library(reshape2)
library(car)
library(openxlsx)
library(tidyverse)
library(aod)
library(gridExtra)


#################################### Estimating mean effect sizes of herbivores among types of restoration
# setup directory
setwd("D:/博士/文章1/1.7/2")
# Load data-------
dat <- read.xlsx("data/pathogen.xlsx", 1)
nrow(dat) #2266

#Delete all pathogens with pathogen positivity less than 5 rodents
data <- dat %>%
  group_by(Pathogen_species) %>%
  filter(sum(Events) >= 5) %>%
  ungroup() 
nrow(data) #2180
# write.xlsx(data, "filtered_excel_file.xlsx")

## Extract exact indexes of non-NA values (classification greater than or equal to 2)
data <- data[-which(is.na(data$Ndvi)),,]
nrow(data) #2173
# write.xlsx(data, "filtered_excel_file.xlsx")
data <- data[-which(is.na(data$Time_interval)),,]
nrow(data) #1723
# write.xlsx(data, "filtered_excel_file.xlsx")
data <- data[-which(data$Rodent_classification=='unknown source'),,]
nrow(data)#949
data <- data[-which(data$Detection_method=='unknown detection method'),,]
nrow(data) #947
data <- data[-which(data$Sample_source%in%c('unclassified sample')==TRUE),,]
nrow(data) #938
# write.xlsx(data, "filtered_excel_file.xlsx")


# Define the names of the columns to be filtered
columns_to_check <- c("Rodent_family", "Time_interval", "Rodent_classification", "Sample_source", 
                      "Detection_method", "Ndvi", "Geographic_division")

# Create a function that filters categories with at least 2 records per column.
filter_data <- function(data, columns) {
  for (column in columns) {
    if (column %in% names(data)) {
      data <- data %>%
        group_by(across(all_of(column))) %>%
        filter(n() >= 2) %>%
        ungroup()
    } else {
      warning(paste("Column", column, "does not exist in the data."))
    }
  }
  return(data)
}

data <- filter_data(data, columns_to_check)
nrow(data) #937
# write.xlsx(filtered_data, "filtered_excel_file1.xlsx")

colnames(data)
# [1] "NO."                      "Rodent_species"           "Rodent_genus"            
# [4] "Rodent_family"            "Events"                   "Total"                   
# [7] "Time_interval"            "Rodent_classification"   "Sample_source"           
# [10] "Detection_method"         "Pathogen_species"         "Ndvi"                    
# [13] "Pathogen_family"          "Pathogen_classiffication" "Province"                
# [16] "Geographic_division"    "Sampling_time"            "Year"                    
# [19] "Literature_sources"       "Author"                   "Title"                   
# [22] "Study"                   


# change order of factor levels
##Geographic_division
data$Geographic_division <- as.factor(data$Geographic_division)
data$Geographic_division <- factor(data$Geographic_division, levels=c('Central China',
                                                                      'East China',
                                                                      'North China',
                                                                      'Northeast China',
                                                                      'Northwest China',
                                                                      'South China',
                                                                      'Southwest China'))
#Ndvi
data$Ndvi <- as.factor(data$Ndvi)
data$Ndvi <- factor(data$Ndvi, levels=c('Peak ndvi','Higher ndvi',
                                        'Medium ndvi','Low ndvi'))

#Rodent_classification
data$Rodent_classification <- as.factor(data$Rodent_classification)
data$Rodent_classification <- factor(data$Rodent_classification, levels=c('Wild','Commensal','Captive'))
head(data)
##Sample_source
data$Sample_source <- as.factor(data$Sample_source)
data$Sample_source <- factor(data$Sample_source, levels=c('Alimentary or intestinal or faeces sample',
                                                          'Blood or serum sample',
                                                          'Body surface sample',
                                                          'Brain sample',
                                                          'Heart sample',
                                                          'Kidney sample',
                                                          'Liver sample',
                                                          'Lung or respiratory tracts',
                                                          'Spleen sample',
                                                          'Mixed tissues'))

##Detection_method
data$Detection_method <- as.factor(data$Detection_method)
data$Detection_method <- factor(data$Detection_method, levels=c('Gross observation or microbes',
                                                                'Molecular diagnostics',
                                                                'Serological examination'))
##Rodent_family
data$Rodent_family <- as.factor(data$Rodent_family)
data$Rodent_family <- factor(data$Rodent_family, levels=c('Caviidae','Chinchillidae','Cricetidae',
                                                          'Dipodidae','Hystricidae','Muridae',
                                                          'Myocastoridae','Nesomyidae','Sciuridae',
                                                          'Spalacidae'))

##Time_interval
data$Time_interval <- as.factor(data$Time_interval)
data$Time_interval <- factor(data$Time_interval, levels=c('1950-2000',
                                                          '2001-2010',
                                                          '2011-2019',
                                                          '2020-2023'))
#calculateI2
i2=function(model){
  
  ## metafor site code for I2
  W=diag(1/model$vi)  ##
  X=model.matrix(model)
  P=W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2=100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  I2=round(I2,2)
  
  ## summarize by each variance component
  allI2=100 * model$sigma2 / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P))) ##"sigma"（σ）通常指的是标准差
  return(list(I2=I2,allI2=allI2))
}

data=data.frame(data,escalc(xi=data$Events,ni=data$Total,measure="PLN"))
data$original_ratio <- exp(data$yi)
# write.xlsx(data,"data_z.xlsx")

# make observation and study-level random effect
data$observation=factor(1:nrow(data))
data$Study=factor(data$Study)

# Fit a mixed-effects model
fit_null <- rma.mv(yi, vi, random=~1|Study/observation, 
                   method="REML", data=data,
                   control=list(optimizer="optim", optmethod="BFGS")) 
fit_null

# Compare models (add first covariate)
fit_var1 <- update(fit_null, mods = ~ Detection_method-1)
fit_var2 <- update(fit_null, mods = ~ Geographic_division-1)
fit_var3 <- update(fit_null, mods = ~ Ndvi-1)
fit_var4 <- update(fit_null, mods = ~ Rodent_classification-1)
fit_var5 <- update(fit_null, mods = ~ Rodent_family-1)
fit_var6 <- update(fit_null, mods = ~ Time_interval-1)
fit_var7 <- update(fit_null, mods = ~ Sample_source-1)

# Check for variable covariance, when GVIF^(1/(2Df)) is greater than 2, it indicates a more severe covariance problem
vif_values <- vif(lm(yi ~ Detection_method + Geographic_division + Ndvi + Rodent_classification + 
                       Rodent_family + Time_interval + Sample_source,
                     data = data))
vif_values

# Setting up the optimizer and optimization methods
control <- list(optimizer="optim", optmethod="L-BFGS-B")

# Setting up the model
fit <- rma.mv(yi, vi, mods=~Detection_method + Geographic_division + Ndvi + Rodent_classification + 
                            Rodent_family + Time_interval + Sample_source-1,
              random=~1|Study/observation, 
              method="REML", data=data,
              control=control,
              verbose=TRUE)

fit 

# Testing the significance of fixed effects in the model
anova(fit)

###############   viruses    ###################################################

data_v <- data %>% filter(Pathogen_classiffication=="Viruses")
nrow(data_v) # 427
# write.xlsx(data_v,"data_v.xlsx")

#Study records for each categorical level of screening predictor variables were no less than 2
data_v <- filter_data(data_v, columns_to_check)
nrow(data_v) # 426         
# write.xlsx(data_v,"data_v1.xlsx")

# change order of factor levels 
##Geographic_division
data_v$Geographic_division <- as.factor(data_v$Geographic_division)
data_v$Geographic_division <- factor(data_v$Geographic_division, levels=c('Central China',
                                                                          'East China',
                                                                          'North China',
                                                                          'Northeast China',
                                                                          'Northwest China',
                                                                          'South China',
                                                                          'Southwest China'))
#Ndvi
data_v$Ndvi <- as.factor(data_v$Ndvi)
data_v$Ndvi <- factor(data_v$Ndvi, levels=c('Peak ndvi','Higher ndvi',
                                            'Medium ndvi','Low ndvi'))
##Rodent_classification
data_v$Rodent_classification <- as.factor(data_v$Rodent_classification)
data_v$Rodent_classification <- factor(data_v$Rodent_classification, levels=c('Wild','Commensal','Captive'))

##Sample_source
data_v$Sample_source <- as.factor(data_v$Sample_source)
data_v$Sample_source <- factor(data_v$Sample_source, levels=c('Alimentary or intestinal or faeces sample',
                                                              'Blood or serum sample',
                                                              'Body surface sample',
                                                              'Brain sample',
                                                              'Heart sample',
                                                              'Kidney sample',
                                                              'Liver sample',
                                                              'Lung or respiratory tracts',
                                                              'Spleen sample',
                                                              'Mixed tissues'))

##Detection_method
data_v$Detection_method <- as.factor(data_v$Detection_method)
data_v$Detection_method <- factor(data_v$Detection_method, levels=c('Gross observation or microbes',
                                                                    'Molecular diagnostics',
                                                                    'Serological examination'))
##Rodent_family
data_v$Rodent_family <- as.factor(data_v$Rodent_family)
data_v$Rodent_family <- factor(data_v$Rodent_family, levels=c('Caviidae','Chinchillidae','Cricetidae',
                                                              'Dipodidae','Hystricidae','Muridae',
                                                              'Myocastoridae','Nesomyidae','Sciuridae',
                                                              'Spalacidae'))

##Time_interval
data_v$Time_interval <- as.factor(data_v$Time_interval)
data_v$Time_interval <- factor(data_v$Time_interval, levels=c('1950-2000',
                                                              '2001-2010',
                                                              '2011-2019',
                                                              '2020-2023'))

# make observation and study-level random effect
data_v$observation=factor(1:nrow(data_v))
data_v$Study=factor(data_v$Study)

## Fit a mixed-effects model
fit_null_v <- rma.mv(yi, vi, random=~1|Study/observation, 
                     method="REML", data=data_v,
                     control=list(optimizer="optim", optmethod="BFGS"))  ## Using experiment nested in publication as the random effect
fit_null_v

# Compare models (add first covariate)
fit_var_v1 <- update(fit_null_v, mods = ~ Detection_method-1)
fit_var_v2 <- update(fit_null_v, mods = ~ Geographic_division-1)
fit_var_v3 <- update(fit_null_v, mods = ~ Ndvi-1)
fit_var_v4 <- update(fit_null_v, mods = ~ Rodent_classification-1)
fit_var_v5 <- update(fit_null_v, mods = ~ Rodent_family-1)
fit_var_v6 <- update(fit_null_v, mods = ~ Time_interval-1)
fit_var_v7 <- update(fit_null_v, mods = ~ Sample_source-1)

# Check for variable covariance, when GVIF^(1/(2Df)) is greater than 2, it indicates a more severe covariance problem
vif_values_v <- vif(lm(yi ~ Detection_method + Geographic_division + Ndvi + Rodent_classification + 
                       Rodent_family + Time_interval + Sample_source,
                     data = data_v))
vif_values_v

fit_v <- rma.mv(yi, vi, mods=~Detection_method + Geographic_division + Rodent_classification + 
                              Rodent_family + Time_interval + Sample_source-1,
                random=~1|Study/observation, 
                method="REML", data=data_v,
                control=control,
                verbose=TRUE) 
fit_v 

# Testing the significance of fixed effects in the model
anova(fit_v) 

## #################      bacteria    ##########################################
data_b <- data %>% filter(Pathogen_classiffication=="Bacteria")
# write.xlsx(data_b,"data_b.xlsx")
nrow(data_b)  

#Study records for each categorical level of screening predictor variables were no less than 2
data_b <- filter_data(data_b, columns_to_check)
nrow(data_b) #  
# write.xlsx(data_b,"data_b1.xlsx")

# change order of factor levels 
##Geographic_division
data_b$Geographic_division <- as.factor(data_b$Geographic_division)
data_b$Geographic_division <- factor(data_b$Geographic_division, levels=c('Central China',
                                                                          'East China',
                                                                          'North China',
                                                                          'Northeast China',
                                                                          'Northwest China',
                                                                          'South China',
                                                                          'Southwest China'))
#Ndvi
data_b$Ndvi <- as.factor(data_b$Ndvi)
data_b$Ndvi <- factor(data_b$Ndvi, levels=c('Peak ndvi','Higher ndvi',
                                            'Medium ndvi','Low ndvi'))
##Rodent_classification
data_b$Rodent_classification <- as.factor(data_b$Rodent_classification)
data_b$Rodent_classification <- factor(data_b$Rodent_classification, levels=c('Wild','Commensal','Captive'))

##Sample_source
data_b$Sample_source <- as.factor(data_b$Sample_source)
data_b$Sample_source <- factor(data_b$Sample_source, levels=c('Alimentary or intestinal or faeces sample',
                                                              'Blood or serum sample',
                                                              'Body surface sample',
                                                              'Brain sample',
                                                              'Heart sample',
                                                              'Kidney sample',
                                                              'Liver sample',
                                                              'Lung or respiratory tracts',
                                                              'Spleen sample',
                                                              'Mixed tissues'))

##Detection_method
data_b$Detection_method <- as.factor(data_b$Detection_method)
data_b$Detection_method <- factor(data_b$Detection_method, levels=c('Gross observation or microbes',
                                                                    'Molecular diagnostics',
                                                                    'Serological examination'))
##Rodent_family
data_b$Rodent_family <- as.factor(data_b$Rodent_family)
data_b$Rodent_family <- factor(data_b$Rodent_family, levels=c('Caviidae','Chinchillidae','Cricetidae',
                                                              'Dipodidae','Hystricidae','Muridae',
                                                              'Myocastoridae','Nesomyidae','Sciuridae',
                                                              'Spalacidae'))

##Time_interval
data_b$Time_interval <- as.factor(data_b$Time_interval)
data_b$Time_interval <- factor(data_b$Time_interval, levels=c('1950-2000',
                                                              '2001-2010',
                                                              '2011-2019',
                                                              '2020-2023'))
# make observation and study-level random effect
data_b$observation=factor(1:nrow(data_b))
data_b$Study=factor(data_b$Study)

## Fit a mixed-effects model
fit_null_b <- rma.mv(yi, vi, random=~1|Study/observation, 
                     method="REML", data=data_b,
                     control=list(optimizer="optim", optmethod="BFGS"))  ## Using experiment nested in publication as the random effect
fit_null_b

# Compare models (add first covariate)
fit_var_b1 <- update(fit_null_b, mods = ~ Detection_method-1)
fit_var_b2 <- update(fit_null_b, mods = ~ Geographic_division-1)
fit_var_b3 <- update(fit_null_b, mods = ~ Ndvi-1)
fit_var_b4 <- update(fit_null_b, mods = ~ Rodent_classification-1)
fit_var_b5 <- update(fit_null_b, mods = ~ Rodent_family-1)
fit_var_b6 <- update(fit_null_b, mods = ~ Time_interval-1)
fit_var_b7 <- update(fit_null_b, mods = ~ Sample_source-1)

# Check for variable covariance, when GVIF^(1/(2Df)) is greater than 2, it indicates a more severe covariance problem
vif_values_b <- vif(lm(yi ~ Detection_method + Geographic_division + Ndvi + Rodent_classification + 
                       Rodent_family + Time_interval + Sample_source,
                     data = data_b))
vif_values_b

fit_b<- rma.mv(yi, vi, mods=~Detection_method + Geographic_division + Rodent_classification + 
                             Rodent_family + Time_interval + Sample_source-1,
               random=~1|Study/observation, 
               method="REML", data=data_b,
               control=control,
               verbose=TRUE) 
fit_b 

# Testing the significance of fixed effects in the model
anova(fit_b) 

#######################    parasites    ########################################
data_p <- data %>% filter(Pathogen_classiffication=="Parasites")
# write.xlsx(data_p,"data_p.xlsx")
nrow(data_p)  ## 226

#Study records for each categorical level of screening predictor variables were no less than 2
data_p <- filter_data(data_p, columns_to_check)
# write.xlsx(data_p,"data_p1.xlsx")
nrow(data_p) # 225 

# change order of factor levels
##Geographic_division
data_p$Geographic_division <- as.factor(data_p$Geographic_division)
data_p$Geographic_division <- factor(data_p$Geographic_division, levels=c('Central China',
                                                                          'East China',
                                                                          'North China',
                                                                          'Northeast China',
                                                                          'Northwest China',
                                                                          'South China',
                                                                          'Southwest China'))
#Ndvi
data_p$Ndvi <- as.factor(data_p$Ndvi)
data_p$Ndvi <- factor(data_p$Ndvi, levels=c('Peak ndvi','Higher ndvi',
                                            'Medium ndvi','Low ndvi'))
##Rodent_classification
data_p$Rodent_classification <- as.factor(data_p$Rodent_classification)
data_p$Rodent_classification <- factor(data_p$Rodent_classification, levels=c('Wild','Commensal','Captive'))

##Sample_source
data_p$Sample_source <- as.factor(data_p$Sample_source)
data_p$Sample_source <- factor(data_p$Sample_source, levels=c('Alimentary or intestinal or faeces sample',
                                                              'Blood or serum sample',
                                                              'Body surface sample',
                                                              'Brain sample',
                                                              'Heart sample',
                                                              'Kidney sample',
                                                              'Liver sample',
                                                              'Lung or respiratory tracts',
                                                              'Spleen sample',
                                                              'Mixed tissues'))

##Detection_method
data_p$Detection_method <- as.factor(data_p$Detection_method)
data_p$Detection_method <- factor(data_p$Detection_method, levels=c('Gross observation or microbes',
                                                                    'Molecular diagnostics',
                                                                    'Serological examination'))
##Rodent_family
data_p$Rodent_family <- as.factor(data_p$Rodent_family)
data_p$Rodent_family <- factor(data_p$Rodent_family, levels=c('Caviidae','Chinchillidae','Cricetidae',
                                                              'Dipodidae','Hystricidae','Muridae',
                                                              'Myocastoridae','Nesomyidae','Sciuridae',
                                                              'Spalacidae'))

##Time_interval
data_p$Time_interval <- as.factor(data_p$Time_interval)
data_p$Time_interval <- factor(data_p$Time_interval, levels=c('1950-2000',
                                                              '2001-2010',
                                                              '2011-2019',
                                                              '2020-2023'))
# make observation and study-level random effect
data_p$observation=factor(1:nrow(data_p))
data_p$Study=factor(data_p$Study)


## Fit a mixed-effects model
fit_null_p <- rma.mv(yi, vi, random=~1|Study/observation, 
                     method="REML", data=data_p,
                     control=list(optimizer="optim", optmethod="BFGS"))  ## Using experiment nested in publication as the random effect
fit_null_p

# Compare models (add first covariate)
fit_var_p1 <- update(fit_null_p, mods = ~ Detection_method-1)
fit_var_p2 <- update(fit_null_p, mods = ~ Geographic_division-1)
fit_var_p3 <- update(fit_null_p, mods = ~ Ndvi-1)
fit_var_p4 <- update(fit_null_p, mods = ~ Rodent_classification-1)
fit_var_p5 <- update(fit_null_p, mods = ~ Rodent_family-1)
fit_var_p6 <- update(fit_null_p, mods = ~ Time_interval-1)
fit_var_p7 <- update(fit_null_p, mods = ~ Sample_source-1)

vif_values_p <- vif(lm(yi ~ Detection_method + Geographic_division + Ndvi + Rodent_classification + 
                         Rodent_family + Time_interval + Sample_source,
                       data = data_p))
vif_values_p

fit_p<- rma.mv(yi, vi, mods=~ Geographic_division  + Ndvi + Rodent_family + 
                              Time_interval + Sample_source-1,
               random=~1|Study/observation, 
               method="REML", data=data_p,
               control=control,
               verbose=TRUE) 
fit_p 

# Testing the significance of fixed effects in the model
anova(fit_p) 
