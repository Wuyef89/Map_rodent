##  Load R package
library(metafor)
library(reshape2)
library(car)
library(openxlsx)
library(tidyverse)
library(aod)
library(gridExtra)


# Load data-------
########################     Total pathogens     ###############################

dat <- read.xlsx("META/Pathogens_China.xlsx", 1)
nrow(dat) #2266

# Remove all pathogens with a pathogen positive count of less than 5
data <- dat %>%
  group_by(Pathogen_species) %>%
  filter(sum(Events) >= 5) %>%
  ungroup() 
nrow(data) #2180
# write.xlsx(data, "filtered_excel_file.xlsx")

## 删除数据中的缺失值
data <- data[-which(is.na(data$Ndvi)),,]
nrow(data) #2173
# write.xlsx(data, "filtered_excel_file.xlsx")
data <- data[-which(is.na(data$Time_interval)),,]
nrow(data) #1723
# write.xlsx(data, "filtered_excel_file.xlsx")

data <- data[-which(data$Rodents_classification=='unknown source'),,]
nrow(data)#949
data <- data[-which(data$Detection_method=='unknown detection method'),,]
nrow(data) #947
data <- data[-which(data$Tissue_source%in%c('unclassified sample')==TRUE),,]
nrow(data) #938
# write.xlsx(data, "filtered_excel_file.xlsx")


# Define the names of the columns to be filtered
columns_to_check <- c("Rodent_family", "Time_interval", "Rodents_classification", "Tissue_source", 
                      "Detection_method", "Ndvi", "Geographical_division")

# No fewer than 2 research records at the categorical level for each predictor of the screening predictor variables
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
# [7] "Time_interval"            "Rodents_classification"   "Tissue_source"           
# [10] "Detection_method"         "Pathogen_species"         "Ndvi"                    
# [13] "Pathogen_family"          "Pathogen_classiffication" "Province"                
# [16] "Geographical_division"    "Sampling_time"            "Year"                    
# [19] "Literature_sources"       "Author"                   "Title"                   
# [22] "Study"                   


# change order of factor levels 
##Geographical_division
data$Geographical_division <- as.factor(data$Geographical_division)
data$Geographical_division <- factor(data$Geographical_division, levels=c('Central China',
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
##Rodents_classification
data$Rodents_classification <- as.factor(data$Rodents_classification)
data$Rodents_classification <- factor(data$Rodents_classification, levels=c('wild','commensal','captive'))

##Tissue_source
data$Tissue_source <- as.factor(data$Tissue_source)
data$Tissue_source <- factor(data$Tissue_source, levels=c('Alimentary or intestinal or faeces sample',
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
data$Rodent_family <- factor(data$Rodent_family, levels=c('Caviidae','Chinchillidae','Cricetidae','
                                                           Dipodidae','Hystricidae','Muridae',
                                                          'Myocastoridae','Nesomyidae','Sciuridae',
                                                          'Spalacidae'))

##Time_interval
data$Time_interval <- as.factor(data$Time_interval)
data$Time_interval <- factor(data$Time_interval, levels=c('1950-2000',
                                                          '2001-2010',
                                                          '2011-2019',
                                                          '2020-2023'))

## function for I2 for rma.mv
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

## Determining the link function
events <- data$Events
n <- data$Total
rate <- transform(data,p = events/n,log = log(events/n),
                  logit = log((events/n)/(1-events/n)),
                  arcsin.size = asin(sqrt(events/(n + 1))),
                  darcsin = 0.5 * (asin(sqrt(events/(n + 1))) +
                                     asin((sqrt(events + 1)/(n + 1)))))
shapiro.test(rate$p)
shapiro.test(rate$log)
shapiro.test(rate$logit)
shapiro.test(rate$arcsin.size)
shapiro.test(rate$darcsin)

#  The codes corresponding to each of the above formulas
# "PR" for the raw proportion,
# "PLN" for the log transformed proportion,
# "PLO" for the logit transformed proportion (i.e., log odds),
# "PAS" for the arcsine square root transformed proportion (i.e., the angular transformation),
# "PFT" for the Freeman-Tukey double arcsine transformed proportion (Freeman & Tukey, 1950).

#   本例数据最适合"PLN" 
# > shapiro.test(rate$p)
# Shapiro-Wilk normality test
# data:  rate$p
# W = 0.70658, p-value < 2.2e-16

# > shapiro.test(rate$log)
# Shapiro-Wilk normality test
# data:  rate$log
# W = 0.98869, p-value = 5.104e-12

# > shapiro.test(rate$logit)
# Shapiro-Wilk normality test
# data:  rate$logit
# W = NaN, p-value = NA

# > shapiro.test(rate$arcsin.size)
# Shapiro-Wilk normality test
# data:  rate$arcsin.size
# W = 0.89088, p-value < 2.2e-16

# > shapiro.test(rate$darcsin)
# Shapiro-Wilk normality test
# data:  rate$darcsin
# W = 0.84098, p-value < 2.2e-16

## make observation and study-level random effect
data$observation=factor(1:nrow(data))
data$Study=factor(data$Study)

## Fit a mixed-effects model（Empty model without intercepts）
fit_null <- rma.mv(yi, vi, random=~1|Study/observation, 
              method="REML", data=data,
              control=list(optimizer="optim", optmethod="BFGS"))  ## Using experiment nested in publication as the random effect
fit_null

## Test the significance of fixed effects for each covariate in the model
fit_var1 <- update(fit_null, mods = ~ Tissue_source-1)
fit_var2 <- update(fit_null, mods = ~ Geographical_division-1)
fit_var3 <- update(fit_null, mods = ~ Rodents_classification-1)
fit_var4 <- update(fit_null, mods = ~ Detection_method-1)
fit_var5 <- update(fit_null, mods = ~ Rodent_family-1)
fit_var6 <- update(fit_null, mods = ~ Time_interval-1)
fit_var7 <- update(fit_null, mods = ~ Ndvi-1)

# Set up optimizers and optimization methods
control <- list(optimizer="optim", optmethod="L-BFGS-B")

# Setting up the model
fit <- rma.mv(yi, vi, mods=~Tissue_source + Geographical_division + Rodents_classification + 
                Detection_method + Rodent_family + Time_interval + Ndvi-1,
              random=~1|Study/observation, 
              method="REML", data=data,
              control=control,
              verbose=TRUE)

fit 

# Check for variable covariance, when GVIF^(1/(2Df)) is greater than 2, it indicates a more severe covariance problem
vif_values <- vif(lm(yi ~ Tissue_source + Geographical_division + Rodents_classification + 
                       Detection_method + Rodent_family + Time_interval + Ndvi,
                     data = data))
print(vif_values)


## Calculate the value of i2
i2_results = i2(fit_null)
print(i2_results$I2)
print(i2_results$allI2)

## Testing the significance of fixed effects in the model
anova(fit)

###############      viruses     ###############################################

data_v <- data %>% filter(Pathogen_classiffication=="Viruses")
nrow(data_v) # 427
# write.xlsx(data_v,"data_v.xlsx")

# No fewer than 2 research records at the categorical level for each predictor of the screening predictor variables
data_v <- filter_data(data_v, columns_to_check)
nrow(data_v) # 426         
# write.xlsx(data_v,"data_v1.xlsx")

# change order of factor levels 
##Geographical_division
data_v$Geographical_division <- as.factor(data_v$Geographical_division)
data_v$Geographical_division <- factor(data_v$Geographical_division, levels=c('Central China',
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
##Rodents_classification
data_v$Rodents_classification <- as.factor(data_v$Rodents_classification)
data_v$Rodents_classification <- factor(data_v$Rodents_classification, levels=c('wild','commensal','captive'))

##Tissue_source
data_v$Tissue_source <- as.factor(data_v$Tissue_source)
data_v$Tissue_source <- factor(data_v$Tissue_source, levels=c('Alimentary or intestinal or faeces sample',
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
data_v$Rodent_family <- factor(data_v$Rodent_family, levels=c('Caviidae','Chinchillidae','Cricetidae','
                                                           Dipodidae','Hystricidae','Muridae',
                                                              'Myocastoridae','Nesomyidae','Sciuridae',
                                                              'Spalacidae'))

##Time_interval
data_v$Time_interval <- as.factor(data_v$Time_interval)
data_v$Time_interval <- factor(data_v$Time_interval, levels=c('1950-2000',
                                                              '2001-2010',
                                                              '2011-2019',
                                                              '2020-2023'))

## 设置随机效应
data_v$observation=factor(1:nrow(data_v))
data_v$Study=factor(data_v$Study)

## Fit a mixed-effects model（Empty model without intercepts）
fit_null_v <- rma.mv(yi, vi, random=~1|Study/observation, 
                   method="REML", data=data_v,
                   control=list(optimizer="optim", optmethod="BFGS"))  ## Using experiment nested in publication as the random effect
fit_null_v

## Test the significance of fixed effects for each covariate in the model
fit_var_v1 <- update(fit_null_v, mods = ~ Tissue_source-1)
fit_var_v2 <- update(fit_null_v, mods = ~ Geographical_division-1)
fit_var_v3 <- update(fit_null_v, mods = ~ Rodents_classification-1)
fit_var_v4 <- update(fit_null_v, mods = ~ Detection_method-1)
fit_var_v5 <- update(fit_null_v, mods = ~ Rodent_family-1)
fit_var_v6 <- update(fit_null_v, mods = ~ Time_interval-1)
fit_var_v7 <- update(fit_null_v, mods = ~ Ndvi-1)


fit_v <- rma.mv(yi, vi, mods=~Tissue_source + Geographical_division + Rodents_classification + 
                  Detection_method + Rodent_family + Time_interval-1,
                random=~1|Study/observation, 
                method="REML", data=data_v,
                control=control,
                verbose=TRUE) 
fit_v 

# Check for variable covariance
vif_values <- vif(lm(yi ~ Tissue_source + Geographical_division + Rodents_classification + 
                       Detection_method + Rodent_family + Time_interval,
                     data = data_v))
print(vif_values)

i2_results_v = i2(fit_v)
print(i2_results_v$I2)
print(i2_results_v$allI2)
anova(fit_v) 

###################     bacteria       #########################################

data_b <- data %>% filter(Pathogen_classiffication=="Bacteria")
# write.xlsx(data_b,"data_b.xlsx")
nrow(data_b)  

# No fewer than 2 research records at the categorical level for each predictor of the screening predictor variables
data_b <- filter_data(data_b, columns_to_check)
nrow(data_b) #  
# write.xlsx(data_b,"data_b1.xlsx")

# change order of factor levels 
##Geographical_division
data_b$Geographical_division <- as.factor(data_b$Geographical_division)
data_b$Geographical_division <- factor(data_b$Geographical_division, levels=c('Central China',
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
##Rodents_classification
data_b$Rodents_classification <- as.factor(data_b$Rodents_classification)
data_b$Rodents_classification <- factor(data_b$Rodents_classification, levels=c('wild','commensal','captive'))

##Tissue_source
data_b$Tissue_source <- as.factor(data_b$Tissue_source)
data_b$Tissue_source <- factor(data_b$Tissue_source, levels=c('Alimentary or intestinal or faeces sample',
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
data_b$Rodent_family <- factor(data_b$Rodent_family, levels=c('Caviidae','Chinchillidae','Cricetidae','
                                                           Dipodidae','Hystricidae','Muridae',
                                                              'Myocastoridae','Nesomyidae','Sciuridae',
                                                              'Spalacidae'))

##Time_interval
data_b$Time_interval <- as.factor(data_b$Time_interval)
data_b$Time_interval <- factor(data_b$Time_interval, levels=c('1950-2000',
                                                              '2001-2010',
                                                              '2011-2019',
                                                              '2020-2023'))
##设置随机效应
data_b$observation=factor(1:nrow(data_b))
data_b$Study=factor(data_b$Study)

## Fit a mixed-effects model（Empty model without intercepts）
fit_null_b <- rma.mv(yi, vi, random=~1|Study/observation, 
                     method="REML", data=data_b,
                     control=list(optimizer="optim", optmethod="BFGS"))  ## Using experiment nested in publication as the random effect
fit_null_b

## Test the significance of fixed effects for each covariate in the model
fit_var_b1 <- update(fit_null_b, mods = ~ Tissue_source-1)
fit_var_b2 <- update(fit_null_b, mods = ~ Geographical_division-1)
fit_var_b3 <- update(fit_null_b, mods = ~ Rodents_classification-1)
fit_var_b4 <- update(fit_null_b, mods = ~ Detection_method-1)
fit_var_b5 <- update(fit_null_b, mods = ~ Rodent_family-1)
fit_var_b6 <- update(fit_null_b, mods = ~ Time_interval-1)
fit_var_b7 <- update(fit_null_b, mods = ~ Ndvi-1)

fit_b<- rma.mv(yi, vi, mods=~Tissue_source + Geographical_division + Rodents_classification + 
                 Detection_method + Rodent_family + Time_interval-1,
               random=~1|Study/observation, 
               method="REML", data=data_b,
               control=control,
               verbose=TRUE) 
fit_b 

# Check for variable covariance
vif_values <- vif(lm(yi ~ Tissue_source + Geographical_division + Rodents_classification + 
                       Detection_method + Rodent_family + Time_interval,
                     data = data_b))
vif_values

i2_results_b = i2(fit_b)
print(i2_results_b$I2)
print(i2_results_b$allI2)
anova(fit_b) 

####################     parasites     ###########################################

data_p <- data %>% filter(Pathogen_classiffication=="Parasites")
# write.xlsx(data_p,"data_p.xlsx")
nrow(data_p)  ## 226

#No fewer than 2 research records at the categorical level for each predictor of the screening predictor variables
data_p <- filter_data(data_p, columns_to_check)
# write.xlsx(data_p,"data_p1.xlsx")
nrow(data_p) # 225 

# change order of factor levels 
##Geographical_division
data_p$Geographical_division <- as.factor(data_p$Geographical_division)
data_p$Geographical_division <- factor(data_p$Geographical_division, levels=c('Central China',
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
##Rodents_classification
data_p$Rodents_classification <- as.factor(data_p$Rodents_classification)
data_p$Rodents_classification <- factor(data_p$Rodents_classification, levels=c('wild','commensal','captive'))

##Tissue_source
data_p$Tissue_source <- as.factor(data_p$Tissue_source)
data_p$Tissue_source <- factor(data_p$Tissue_source, levels=c('Alimentary or intestinal or faeces sample',
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
data_p$Rodent_family <- factor(data_p$Rodent_family, levels=c('Caviidae','Chinchillidae','Cricetidae','
                                                           Dipodidae','Hystricidae','Muridae',
                                                              'Myocastoridae','Nesomyidae','Sciuridae',
                                                              'Spalacidae'))

##Time_interval
data_p$Time_interval <- as.factor(data_p$Time_interval)
data_p$Time_interval <- factor(data_p$Time_interval, levels=c('1950-2000',
                                                              '2001-2010',
                                                              '2011-2019',
                                                              '2020-2023'))
## Setting up random effects
data_p$observation=factor(1:nrow(data_p))
data_p$Study=factor(data_p$Study)


## Fit a mixed-effects model（null model）
fit_null_p <- rma.mv(yi, vi, random=~1|Study/observation, 
                     method="REML", data=data_p,
                     control=list(optimizer="optim", optmethod="BFGS"))  ## Using experiment nested in publication as the random effect
fit_null_p

## Test the significance of fixed effects for each covariate in the model
fit_var_p1 <- update(fit_null_p, mods = ~ Tissue_source-1)
fit_var_p2 <- update(fit_null_p, mods = ~ Geographical_division-1)
fit_var_p3 <- update(fit_null_p, mods = ~ Rodents_classification-1) #
fit_var_p4 <- update(fit_null_p, mods = ~ Detection_method-1) #
fit_var_p5 <- update(fit_null_p, mods = ~ Rodent_family-1) 
fit_var_p6 <- update(fit_null_p, mods = ~ Time_interval-1)
fit_var_p7 <- update(fit_null_p, mods = ~ Ndvi-1)

fit_p<- rma.mv(yi, vi, mods=~Tissue_source + Geographical_division + Rodent_family + 
                             Time_interval + Ndvi-1,
               random=~1|Study/observation, 
               method="REML", data=data_p,
               control=control,
               verbose=TRUE) 
fit_p 

# Check for variable covariance, when GVIF^(1/(2Df)) is greater than 2, it indicates a more severe covariance problem
vif_values <- vif(lm(yi ~ Tissue_source + Geographical_division + 
                       Rodent_family + Time_interval + Ndvi,
                     data = data_p))
print(vif_values)

i2_results_p = i2(fit_p)
print(i2_results_p$I2)
print(i2_results_p$allI2)

## Test the significance of fixed effects in the model
anova(fit_p) 

