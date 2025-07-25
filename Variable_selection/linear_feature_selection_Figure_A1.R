source('Functions/manuscript_functions.R')

packages = c("rlist","dplyr","tidyverse",
             "raster", "sf", "stringr", 
             "ggplot2", "ggpubr","RColorBrewer", "colorspace","MASS","plot3D","RobustGaSP","elevatr",
             "lattice", "viridisLite", "latticeExtra", "grid", "reshape", "reshape2")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


##1. Loading data rasters 

# necessary values
desired_layers <- length(2003:2020) * 12 #number of months between 2003 and 2020
first_layer <- length(1895:2020)*12 - desired_layers+1 #first layer is 2003 January

admin <- st_read('Data/fourcorners_states.shp')  # state boundaries

# NDVI raster
modis <- stack('Data/modis_cube_processed.grd') # loads NDVI raster. Available only for years 2003-2020

# Precipitation raster
precip_long <- stack('Data/precip_reproj_trim_res.grd') # loads precipitation raster
precip <- precip_long[[first_layer:nlayers(precip_long)]] # keep only 2003-2020


#VPD raster
vpdmax_long <- stack('Data/vpdmax_reproj_trim_res.grd') # loads VPD raster
vpdmax <- vpdmax_long[[first_layer:nlayers(vpdmax_long)]] # keep only 2003-2020

##2. Transform rasters to arrays so easier to work with. 100x140x216 (216 = 18 years x 12 months)
y <- as.array(modis) # NDVI
x1 <- as.array(precip) # precipitation
x2 <- as.array(vpdmax) # VPD

##3. Fill NA values

loc_ne=get_ne(y) ##y is ndvi
output_NDVI_all_mat_missing_filled=fill_NA(output=y)
precipitatian_all_mat_missing_filled=fill_NA(output=x1)
VPD_all_mat_missing_filled=fill_NA(output=x2)

#0 since we got rid of NA
sum(is.na(output_NDVI_all_mat_missing_filled))
sum(is.na(precipitatian_all_mat_missing_filled))
sum(is.na(VPD_all_mat_missing_filled))

##4. Compute monthly averages

list_eda_ndvi=eda_corr(output_ndvi=output_NDVI_all_mat_missing_filled,
                       input_vpd=VPD_all_mat_missing_filled,
                       input_precipitation=precipitatian_all_mat_missing_filled,num_yrs_all=18)

num_loc=dim(output_NDVI_all_mat_missing_filled)[1] # number of locations/pixels/grids
num_obs_all=dim(output_NDVI_all_mat_missing_filled)[2] # number of months total 
num_yrs_all=num_obs_all/12 # number of years total 
Aug_index=8+12*(0:(num_yrs_all-1))


# R2 Heatmaps

#Linear model for NDVI with vpd and precipitation input, attribution model
lm_one_year_ahead <-function(testing_index=11:18,output_all,
                             input_1=NA,input_2=NA,input_3=NA
){
  sse_training = rep(NA, num_loc)
  l2_pred_error_lm_model_record = matrix(NA,num_loc,length(testing_index))
  pred_lm_model_record=matrix(NA,num_loc,length(testing_index))
  interval_length_lm_model_record=matrix(NA,num_loc,length(testing_index))
  prop95_lm_model_record=matrix(NA,num_loc,length(testing_index))

  i_test_start=testing_index[1]
  input=cbind(input_1[1,1:(i_test_start-1)])
  if(!is.na(input_2)[1]){
    input=cbind(input,input_2[1,1:(i_test_start-1)])
  }
  if(!is.na(input_3)[1]){
    input=cbind(input,input_3[1,1:(i_test_start-1)])
  }
  p=dim(input)[2]

  for(i_loc in 1:num_loc){

    for(i_test in testing_index){
      output=as.vector(output_all[i_loc,1:(i_test-1)])

      input=cbind(output, input_1[i_loc,1:(i_test-1)])
      testing_input=cbind(input_1[i_loc,i_test])
      if(p==2){
        input=cbind(input,input_2[i_loc,1:(i_test-1)])
        testing_input=cbind(testing_input,input_2[i_loc,i_test])
      }else if(p==3){ ##p=3
        input=cbind(input,input_2[i_loc,1:(i_test-1)])
        input=cbind(input,input_3[i_loc,1:(i_test-1)])
        testing_input=cbind(testing_input,input_2[i_loc,i_test])
        testing_input=cbind(testing_input,input_3[i_loc,i_test])
      }
      # make a data frame for linear models
      input = as.data.frame(input)
      colnames(input) = c('output', paste('input',1:p, sep='_'))
      testing_input = as.data.frame(testing_input)
      colnames(testing_input) = paste('input',1:p, sep='_')
      m_i_loc=lm(output~., data=input)

      if(i_test == i_test_start){
        sse_training[i_loc] = sum(residuals(m_i_loc)^2)
      }

      m_i_loc_pred=predict(m_i_loc, newdata = testing_input,
                           interval='predict')

      pred_lm_model_record[i_loc, i_test-i_test_start+1]=m_i_loc_pred[,1]
      l2_pred_error_lm_model_record[i_loc, i_test-i_test_start+1]= (m_i_loc_pred[,1]-output_all[i_loc,i_test])^2
      interval_length_lm_model_record[i_loc, i_test-i_test_start+1] = m_i_loc_pred[,3]-m_i_loc_pred[,2]
      prop95_lm_model_record[i_loc,i_test-i_test_start+1] = (m_i_loc_pred[,3]>output_all[i_loc,i_test])*(m_i_loc_pred[,2]<output_all[i_loc,i_test])
    }
  }

  return.list=list()
  return.list$pred_lm_model= pred_lm_model_record
  return.list$l2_pred_error_lm_model=l2_pred_error_lm_model_record
  return.list$interval_length_lm_model=interval_length_lm_model_record
  return.list$prop95_lm_model=prop95_lm_model_record
  return.list$sse_training_vector = sse_training
  return(return.list)
}

#every month index
num_yrs_all=18
Aug_index=8+12*(0:(num_yrs_all-1))
July_index=7+12*(0:(num_yrs_all-1))
June_index=6+12*(0:(num_yrs_all-1))
May_index=5+12*(0:(num_yrs_all-1))
April_index=4+12*(0:(num_yrs_all-1))
March_index=3+12*(0:(num_yrs_all-1))
Feb_index=2+12*(0:(num_yrs_all-1))
Jan_index=1+12*(0:(num_yrs_all-1))
num_loc=dim(output_NDVI_all_mat_missing_filled)[1]

output_ndvi= output_NDVI_all_mat_missing_filled

input_vpd= VPD_all_mat_missing_filled
input_precip = precipitatian_all_mat_missing_filled


VPD_loc_JanAug_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_FebAug_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_MarchAug_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_AprilAug_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_MayAug_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_JuneAug_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_JulyAug_mean_record=matrix(NA,num_loc,num_yrs_all)

Precip_loc_JanAug_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_FebAug_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_MarchAug_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_AprilAug_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_MayAug_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_JuneAug_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_JulyAug_mean_record=matrix(NA,num_loc,num_yrs_all)

##vpd averages, all including august
for(i_loc in 1:num_loc){
  ##VPD
  #Jan to August
  VPD_i_loc_Jan_to_Aug_mean=rowMeans(cbind(input_vpd[i_loc,Jan_index],input_vpd[i_loc,Feb_index],
                                           input_vpd[i_loc,March_index],input_vpd[i_loc,April_index],
                                           input_vpd[i_loc,May_index],input_vpd[i_loc,June_index],
                                           input_vpd[i_loc,July_index],input_vpd[i_loc,Aug_index]))
  VPD_loc_JanAug_mean_record[i_loc,]=VPD_i_loc_Jan_to_Aug_mean

  #Feb to August
  VPD_i_loc_Feb_to_Aug_mean=rowMeans(cbind(input_vpd[i_loc,Feb_index],
                                           input_vpd[i_loc,March_index],input_vpd[i_loc,April_index],
                                           input_vpd[i_loc,May_index],input_vpd[i_loc,June_index],
                                           input_vpd[i_loc,July_index],input_vpd[i_loc,Aug_index]))
  VPD_loc_FebAug_mean_record[i_loc,]=VPD_i_loc_Feb_to_Aug_mean

  #March to August
  VPD_i_loc_March_to_Aug_mean=rowMeans(cbind(input_vpd[i_loc,March_index],input_vpd[i_loc,April_index],
                                             input_vpd[i_loc,May_index],input_vpd[i_loc,June_index],
                                             input_vpd[i_loc,July_index],input_vpd[i_loc,Aug_index]))
  VPD_loc_MarchAug_mean_record[i_loc,]=VPD_i_loc_March_to_Aug_mean

  #April to August
  VPD_i_loc_April_to_Aug_mean=rowMeans(cbind(input_vpd[i_loc,April_index],
                                             input_vpd[i_loc,May_index],input_vpd[i_loc,June_index],
                                             input_vpd[i_loc,July_index],input_vpd[i_loc,Aug_index]))
  VPD_loc_AprilAug_mean_record[i_loc,]=VPD_i_loc_April_to_Aug_mean

  #May to August
  VPD_i_loc_May_to_Aug_mean=rowMeans(cbind(input_vpd[i_loc,May_index],input_vpd[i_loc,June_index],
                                           input_vpd[i_loc,July_index],input_vpd[i_loc,Aug_index]))
  VPD_loc_MayAug_mean_record[i_loc,]=VPD_i_loc_May_to_Aug_mean

  #June to August
  VPD_i_loc_June_to_Aug_mean=rowMeans(cbind(input_vpd[i_loc,June_index],
                                            input_vpd[i_loc,July_index],input_vpd[i_loc,Aug_index]))
  VPD_loc_JuneAug_mean_record[i_loc,]=VPD_i_loc_June_to_Aug_mean

  #July to August
  VPD_i_loc_July_to_Aug_mean=rowMeans(cbind(input_vpd[i_loc,July_index],input_vpd[i_loc,Aug_index]))
  VPD_loc_JulyAug_mean_record[i_loc,]=VPD_i_loc_July_to_Aug_mean

  #results
  vpd.mean.list1=list()
  vpd.mean.list1$VPD_loc_JanAug_mean_record = VPD_loc_JanAug_mean_record
  vpd.mean.list1$VPD_loc_FebAug_mean_record = VPD_loc_FebAug_mean_record
  vpd.mean.list1$VPD_loc_MarchAug_mean_record = VPD_loc_MarchAug_mean_record
  vpd.mean.list1$VPD_loc_AprilAug_mean_record = VPD_loc_AprilAug_mean_record
  vpd.mean.list1$VPD_loc_MayAug_mean_record = VPD_loc_MayAug_mean_record
  vpd.mean.list1$VPD_loc_JuneAug_mean_record = VPD_loc_JuneAug_mean_record
  vpd.mean.list1$VPD_loc_JulyAug_mean_record = VPD_loc_JulyAug_mean_record

}


## precipitation averages all including august
for(i_loc in 1:num_loc){
  ##Preci
  #Jan to August
  Precip_i_loc_Jan_to_Aug_mean=rowMeans(cbind(input_precip[i_loc,Jan_index],input_precip[i_loc,Feb_index],
                                              input_precip[i_loc,March_index],input_precip[i_loc,April_index],
                                              input_precip[i_loc,May_index],input_precip[i_loc,June_index],
                                              input_precip[i_loc,July_index],input_precip[i_loc,Aug_index]))
  Precip_loc_JanAug_mean_record[i_loc,]=Precip_i_loc_Jan_to_Aug_mean

  #Feb to August
  Precip_i_loc_Feb_to_Aug_mean=rowMeans(cbind(input_precip[i_loc,Feb_index],
                                              input_precip[i_loc,March_index],input_precip[i_loc,April_index],
                                              input_precip[i_loc,May_index],input_precip[i_loc,June_index],
                                              input_precip[i_loc,July_index],input_precip[i_loc,Aug_index]))
  Precip_loc_FebAug_mean_record[i_loc,]=Precip_i_loc_Feb_to_Aug_mean


  #March to August
  Precip_i_loc_March_to_Aug_mean=rowMeans(cbind(input_precip[i_loc,March_index],input_precip[i_loc,April_index],
                                                input_precip[i_loc,May_index],input_precip[i_loc,June_index],
                                                input_precip[i_loc,July_index],input_precip[i_loc,Aug_index]))
  Precip_loc_MarchAug_mean_record[i_loc,]=Precip_i_loc_March_to_Aug_mean

  #April to August
  Precip_i_loc_April_to_Aug_mean=rowMeans(cbind(input_precip[i_loc,April_index],
                                                input_precip[i_loc,May_index],input_precip[i_loc,June_index],
                                                input_precip[i_loc,July_index],input_precip[i_loc,Aug_index]))
  Precip_loc_AprilAug_mean_record[i_loc,]=Precip_i_loc_April_to_Aug_mean

  #May to August
  Precip_i_loc_May_to_Aug_mean=rowMeans(cbind(input_precip[i_loc,May_index],input_precip[i_loc,June_index],
                                              input_precip[i_loc,July_index],input_precip[i_loc,Aug_index]))
  Precip_loc_MayAug_mean_record[i_loc,]=Precip_i_loc_May_to_Aug_mean

  #June to August
  Precip_i_loc_June_to_Aug_mean=rowMeans(cbind(input_precip[i_loc,June_index],
                                               input_precip[i_loc,July_index],input_precip[i_loc,Aug_index]))
  Precip_loc_JuneAug_mean_record[i_loc,]=Precip_i_loc_June_to_Aug_mean

  #July to August
  Precip_i_loc_July_to_Aug_mean=rowMeans(cbind(input_precip[i_loc,July_index],input_precip[i_loc,Aug_index]))
  Precip_loc_JulyAug_mean_record[i_loc,]=Precip_i_loc_July_to_Aug_mean

  #results
  precip.mean.list1=list()
  precip.mean.list1$Precip_loc_JanAug_mean_record = Precip_loc_JanAug_mean_record
  precip.mean.list1$Precip_loc_FebAug_mean_record = Precip_loc_FebAug_mean_record
  precip.mean.list1$Precip_loc_MarchAug_mean_record = Precip_loc_MarchAug_mean_record
  precip.mean.list1$Precip_loc_AprilAug_mean_record = Precip_loc_AprilAug_mean_record
  precip.mean.list1$Precip_loc_MayAug_mean_record = Precip_loc_MayAug_mean_record
  precip.mean.list1$Precip_loc_JuneAug_mean_record = Precip_loc_JuneAug_mean_record
  precip.mean.list1$Precip_loc_JulyAug_mean_record = Precip_loc_JulyAug_mean_record

}


results = list()
r2_training = list()
r2_test = list()
###  long stretch of months with linear model
combinations <- expand.grid(vpd = names(vpd.mean.list1), precip = names(precip.mean.list1))

#counter
count <- 0

  
# apply for each combination
run_model <- function(vpd, precip) {
  count <<- count + 1

  lm_result <- lm_one_year_ahead(testing_index = 11:18,
                                 output_all = output_NDVI_all_mat_missing_filled[, Aug_index],
                                 input_1 = vpd.mean.list1[[vpd]][, 1:18],
                                 input_2 = precip.mean.list1[[precip]][, 1:18],
                                 input_3 = NA)
  
  r2_train_function <- function(i){
    1-(lm_result$sse_training_vector[i]/sum((output_NDVI_all_mat_missing_filled[,Aug_index][i,1:10]-mean(output_NDVI_all_mat_missing_filled[,Aug_index][i,1:10]))^2))
  }
  
  grand_r2 = mean(sapply(1:14000, r2_train_function))
 
  print(paste("Model run", count, "of", nrow(combinations), ": R-squared Train =", grand_r2))

  return(c(r2_train = grand_r2))
}

#apply the function to each combination
results <- mapply(run_model, combinations$vpd, combinations$precip, SIMPLIFY = FALSE)

# results to a data frame
results_df <- do.call(rbind, results)

results_df <- as.data.frame(results_df)

#column names for the new data frame
colnames(results_df) <- c("r2_train")

# results to matrix format for heatmap plotting
r2_train_matrix <- matrix(results_df$r2_train, nrow = length(vpd.mean.list1), byrow = FALSE)

#dimension names for the matrices
dimnames(r2_train_matrix) <- list(names(vpd.mean.list1), names(precip.mean.list1))

# matrices to data frames for plotting
df_train_heatmap <- melt(r2_train_matrix, varnames = c("VPD", "Precip"))

lm_aug_R2_train <- df_train_heatmap %>% arrange(desc(value))

#Figure A1a
p_train <- ggplot(df_train_heatmap, aes(x = VPD, y = Precip, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(c(0,mean(df_train_heatmap$value))), limit = c(0, max(df_train_heatmap$value))) +
  scale_x_discrete(labels = c("Jan-Aug", "Feb-Aug", "March-Aug", "April-Aug", "May-Aug", "June-Aug", "July-Aug")) +
  scale_y_discrete(labels = c("Jan-Aug", "Feb-Aug", "March-Aug", "April-Aug", "May-Aug", "June-Aug", "July-Aug")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Monthly Ranges of VPD")+
  ylab("Monthly Ranges of Precipitation")+
  theme(axis.text.x = element_text( size =14, ),
        axis.text.y = element_text( size = 14),  
        axis.title.x = element_text( size = 18),
        axis.title.y = element_text( size = 18))+
  theme(legend.text = element_text(size = 18))+
  theme(legend.box.just = "center")+
  theme(legend.key.height = unit(2, "cm"))+
  geom_text(aes(label = round(value, 2))) 


#alternative heatmaps

#VPD_loc_Jan_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_Feb_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_March_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_April_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_May_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_June_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_July_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_Aug_mean_record=matrix(NA,num_loc,num_yrs_all)


#Precip_loc_Jan_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_Feb_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_March_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_April_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_May_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_June_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_July_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_Aug_mean_record=matrix(NA,num_loc,num_yrs_all)

#just a single month Range
for(i_loc in 1:num_loc){
  ##VPD
  #Feb 
  VPD_i_loc_Feb_mean=rowMeans(cbind(input_vpd[i_loc,Feb_index]))
  VPD_loc_Feb_mean_record[i_loc,]=VPD_i_loc_Feb_mean
  
  #March 
  VPD_i_loc_March_mean=rowMeans(cbind(input_vpd[i_loc,March_index]))
  VPD_loc_March_mean_record[i_loc,]=VPD_i_loc_March_mean
  
  #April
  VPD_i_loc_April_mean=rowMeans(cbind(input_vpd[i_loc,April_index]))
  VPD_loc_April_mean_record[i_loc,]=VPD_i_loc_April_mean
  
  #May
  VPD_i_loc_May_mean=rowMeans(cbind(input_vpd[i_loc,May_index]))
  VPD_loc_May_mean_record[i_loc,]=VPD_i_loc_May_mean
  
  #June 
  VPD_i_loc_June_mean=rowMeans(cbind(input_vpd[i_loc,June_index]))
  VPD_loc_June_mean_record[i_loc,]=VPD_i_loc_June_mean
  
  #July 
  VPD_i_loc_July_mean=rowMeans(cbind(input_vpd[i_loc,July_index]))
  VPD_loc_July_mean_record[i_loc,]=VPD_i_loc_July_mean
  
  #Aug 
  VPD_i_loc_Aug_mean=rowMeans(cbind(input_vpd[i_loc,Aug_index]))
  VPD_loc_Aug_mean_record[i_loc,]=VPD_i_loc_Aug_mean
  
  #results
  vpd.mean.list2=list()
  vpd.mean.list2$VPD_loc_Feb_mean_record = VPD_loc_Feb_mean_record
  vpd.mean.list2$VPD_loc_March_mean_record = VPD_loc_March_mean_record
  vpd.mean.list2$VPD_loc_April_mean_record = VPD_loc_April_mean_record
  vpd.mean.list2$VPD_loc_May_mean_record = VPD_loc_May_mean_record
  vpd.mean.list2$VPD_loc_June_mean_record = VPD_loc_June_mean_record
  vpd.mean.list2$VPD_loc_July_mean_record = VPD_loc_July_mean_record
  vpd.mean.list2$VPD_loc_Aug_mean_record = VPD_loc_Aug_mean_record
}


for(i_loc in 1:num_loc){
  ##Preci

  #Feb
  Precip_i_loc_Feb_mean=rowMeans(cbind(input_precip[i_loc,Feb_index]))
  Precip_loc_Feb_mean_record[i_loc,]=Precip_i_loc_Feb_mean
  
  
  #March 
  Precip_i_loc_March_mean=rowMeans(cbind(input_precip[i_loc,March_index]))
  Precip_loc_March_mean_record[i_loc,]=Precip_i_loc_March_mean
  
  #April
  Precip_i_loc_April_mean=rowMeans(cbind(input_precip[i_loc,April_index]))
  Precip_loc_April_mean_record[i_loc,]=Precip_i_loc_April_mean
  
  #May
  Precip_i_loc_May_mean=rowMeans(cbind(input_precip[i_loc,May_index]))
  Precip_loc_May_mean_record[i_loc,]=Precip_i_loc_May_mean
  
  #June 
  Precip_i_loc_June_mean=rowMeans(cbind(input_precip[i_loc,June_index]))
  Precip_loc_June_mean_record[i_loc,]=Precip_i_loc_June_mean
  
  #July 
  Precip_i_loc_July_mean=rowMeans(cbind(input_precip[i_loc,July_index]))
  Precip_loc_July_mean_record[i_loc,]=Precip_i_loc_July_mean
  
  #Aug 
  Precip_i_loc_Aug_mean=rowMeans(cbind(input_precip[i_loc,Aug_index]))
  Precip_loc_Aug_mean_record[i_loc,]=Precip_i_loc_Aug_mean
  
  #results
  precip.mean.list2=list()
  #precip.mean.list2$Precip_loc_Jan_mean_record = Precip_loc_Jan_mean_record
  precip.mean.list2$Precip_loc_Feb_mean_record = Precip_loc_Feb_mean_record
  precip.mean.list2$Precip_loc_March_mean_record = Precip_loc_March_mean_record
  precip.mean.list2$Precip_loc_April_mean_record = Precip_loc_April_mean_record
  precip.mean.list2$Precip_loc_May_mean_record = Precip_loc_May_mean_record
  precip.mean.list2$Precip_loc_June_mean_record = Precip_loc_June_mean_record
  precip.mean.list2$Precip_loc_July_mean_record = Precip_loc_July_mean_record
  precip.mean.list2$Precip_loc_Aug_mean_record = Precip_loc_Aug_mean_record
  
}


results = list()
r2_training = list()
r2_test = list()
#  long stretch of months with linear model
combinations2 <- expand.grid(vpd = names(vpd.mean.list2), precip = names(precip.mean.list2))

#counter
count <- 0

# apply for each combination
run_model_2 <- function(vpd, precip) {
  count <<- count + 1
  
  lm_result <- lm_one_year_ahead(testing_index = 11:18,
                                 output_all = output_NDVI_all_mat_missing_filled[, Aug_index],
                                 input_1 = vpd.mean.list2[[vpd]][, 1:18],
                                 input_2 = precip.mean.list2[[precip]][, 1:18],
                                 input_3 = NA)

  
  r2_train_function <- function(i){
    1-(lm_result$sse_training_vector[i]/sum((output_NDVI_all_mat_missing_filled[,Aug_index][i,1:10]-mean(output_NDVI_all_mat_missing_filled[,Aug_index][i,1:10]))^2))
  }
  
  grand_r2 = mean(sapply(1:14000, r2_train_function))
  
  print(paste("Model run", count, "of", nrow(combinations2), ": R-squared Train =", grand_r2))
  
  return(c(r2_train = grand_r2))
}

#apply the function to each combination
results_2 <- mapply(run_model_2, combinations2$vpd, combinations2$precip, SIMPLIFY = FALSE)

results_df2 <- do.call(rbind, results_2)

results_df2 <- as.data.frame(results_df2)

colnames(results_df2) <- c("r2_train")

# results to matrix format for heatmap plotting
r2_train_matrix2 <- matrix(results_df2$r2_train, nrow = length(vpd.mean.list2), byrow = FALSE)

#dimension names for the matrices
dimnames(r2_train_matrix2) <- list(names(vpd.mean.list2), names(precip.mean.list2))

df_train_heatmap2 <- melt(r2_train_matrix2, varnames = c("VPD", "Precip"))

lm_R2_train <- df_train_heatmap2 %>% arrange(desc(value))

#Figure A1b
p_train_2 <- ggplot(df_train_heatmap2, aes(x = VPD, y = Precip, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(c(0,mean(df_train_heatmap$value))), limit = c(0, max(df_train_heatmap$value))) +
  #scale_x_discrete(labels = c("Jan", "Feb", "March", "April", "May", "June", "July")) +
  #scale_y_discrete(labels = c("Jan", "Feb", "March", "April", "May", "June", "July")) +
  scale_x_discrete(labels = c( "Feb", "March", "April", "May", "June", "July", "Aug")) +
  scale_y_discrete(labels = c( "Feb", "March", "April", "May", "June", "July", "Aug")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #ggtitle("Heatmap of R-squared for Linear Model Training Data")+
  xlab("Monthly Ranges of VPD")+
  ylab("Monthly Ranges of Precipitation")+
  #labs(fill="R-squared")+
  theme(axis.text.x = element_text( size =14, ),
        axis.text.y = element_text( size = 14),  
        axis.title.x = element_text( size = 18),
        axis.title.y = element_text( size = 18))+
  theme(legend.text = element_text(size = 18))+
  #theme(legend.title = element_text(size = 18))+
  theme(legend.box.just = "center")+
  theme(legend.key.height = unit(2, "cm"))+
  geom_text(aes(label = round(value, 2))) 

#more averages for vpd and precip
VPD_loc_JanFeb_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_FebMarch_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_MarchApril_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_AprilMay_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_MayJune_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_JuneJuly_mean_record=matrix(NA,num_loc,num_yrs_all)
VPD_loc_JulyAug_mean_record=matrix(NA,num_loc,num_yrs_all)

Precip_loc_JanFeb_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_FebMarch_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_MarchApril_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_AprilMay_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_MayJune_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_JuneJuly_mean_record=matrix(NA,num_loc,num_yrs_all)
Precip_loc_JulyAug_mean_record=matrix(NA,num_loc,num_yrs_all)


for(i_loc in 1:num_loc){
  ##VPD
  #Jan to Feb
  VPD_i_loc_Jan_to_Feb_mean=rowMeans(cbind(input_vpd[i_loc,Jan_index],input_vpd[i_loc,Feb_index]))
  VPD_loc_JanFeb_mean_record[i_loc,]=VPD_i_loc_Jan_to_Feb_mean
  
  #Feb to March
  VPD_i_loc_Feb_to_March_mean=rowMeans(cbind(input_vpd[i_loc,Feb_index],
                                           input_vpd[i_loc,March_index]))
  VPD_loc_FebMarch_mean_record[i_loc,]=VPD_i_loc_Feb_to_March_mean
  
  #March to April
  VPD_i_loc_March_to_April_mean=rowMeans(cbind(input_vpd[i_loc,March_index],input_vpd[i_loc,April_index]))
  VPD_loc_MarchApril_mean_record[i_loc,]=VPD_i_loc_March_to_April_mean
  
  #April to May
  VPD_i_loc_April_to_May_mean=rowMeans(cbind(input_vpd[i_loc,April_index],
                                             input_vpd[i_loc,May_index]))
  VPD_loc_AprilMay_mean_record[i_loc,]=VPD_i_loc_April_to_May_mean
  
  #May to June
  VPD_i_loc_May_to_June_mean=rowMeans(cbind(input_vpd[i_loc,May_index],input_vpd[i_loc,June_index]))
  VPD_loc_MayJune_mean_record[i_loc,]=VPD_i_loc_May_to_June_mean
  
  #June to July
  VPD_i_loc_June_to_July_mean=rowMeans(cbind(input_vpd[i_loc,June_index],
                                            input_vpd[i_loc,July_index]))
  VPD_loc_JuneJuly_mean_record[i_loc,]=VPD_i_loc_June_to_July_mean
  
  #July to August
  VPD_i_loc_July_to_Aug_mean=rowMeans(cbind(input_vpd[i_loc,July_index],input_vpd[i_loc,Aug_index]))
  VPD_loc_JulyAug_mean_record[i_loc,]=VPD_i_loc_July_to_Aug_mean
  
  #results
  vpd.mean.list3=list()
  vpd.mean.list3$VPD_loc_JanFeb_mean_record = VPD_loc_JanFeb_mean_record
  vpd.mean.list3$VPD_loc_FebMarch_mean_record = VPD_loc_FebMarch_mean_record
  vpd.mean.list3$VPD_loc_MarchApril_mean_record = VPD_loc_MarchApril_mean_record
  vpd.mean.list3$VPD_loc_AprilMay_mean_record = VPD_loc_AprilMay_mean_record
  vpd.mean.list3$VPD_loc_MayJune_mean_record = VPD_loc_MayJune_mean_record
  vpd.mean.list3$VPD_loc_JuneJuly_mean_record = VPD_loc_JuneJuly_mean_record
  vpd.mean.list3$VPD_loc_JulyAug_mean_record = VPD_loc_JulyAug_mean_record
  
}


for(i_loc in 1:num_loc){
  ##Preci
  #Jan to Feb
  Precip_i_loc_Jan_to_Feb_mean=rowMeans(cbind(input_precip[i_loc,Jan_index],input_precip[i_loc,Feb_index]))
  Precip_loc_JanFeb_mean_record[i_loc,]=Precip_i_loc_Jan_to_Feb_mean
  
  #Feb to March
  Precip_i_loc_Feb_to_March_mean=rowMeans(cbind(input_precip[i_loc,Feb_index],
                                              input_precip[i_loc,March_index]))
  Precip_loc_FebMarch_mean_record[i_loc,]=Precip_i_loc_Feb_to_March_mean
  
  
  #March to April
  Precip_i_loc_March_to_April_mean=rowMeans(cbind(input_precip[i_loc,March_index],input_precip[i_loc,April_index]))
  Precip_loc_MarchApril_mean_record[i_loc,]=Precip_i_loc_March_to_April_mean
  
  #April to May
  Precip_i_loc_April_to_May_mean=rowMeans(cbind(input_precip[i_loc,April_index],
                                                input_precip[i_loc,May_index]))
  Precip_loc_AprilMay_mean_record[i_loc,]=Precip_i_loc_April_to_May_mean
  
  #May to June
  Precip_i_loc_May_to_June_mean=rowMeans(cbind(input_precip[i_loc,May_index],input_precip[i_loc,June_index]))
  Precip_loc_MayJune_mean_record[i_loc,]=Precip_i_loc_May_to_June_mean
  
  #June to July
  Precip_i_loc_June_to_July_mean=rowMeans(cbind(input_precip[i_loc,June_index],
                                               input_precip[i_loc,July_index]))
  Precip_loc_JuneJuly_mean_record[i_loc,]=Precip_i_loc_June_to_July_mean
  
  #July to August
  Precip_i_loc_July_to_Aug_mean=rowMeans(cbind(input_precip[i_loc,July_index],input_precip[i_loc,Aug_index]))
  Precip_loc_JulyAug_mean_record[i_loc,]=Precip_i_loc_July_to_Aug_mean
  
  #results
  precip.mean.list3=list()
  precip.mean.list3$Precip_loc_JanFeb_mean_record = Precip_loc_JanFeb_mean_record
  precip.mean.list3$Precip_loc_FebMarch_mean_record = Precip_loc_FebMarch_mean_record
  precip.mean.list3$Precip_loc_MarchApril_mean_record = Precip_loc_MarchApril_mean_record
  precip.mean.list3$Precip_loc_AprilMay_mean_record = Precip_loc_AprilMay_mean_record
  precip.mean.list3$Precip_loc_MayJune_mean_record = Precip_loc_MayJune_mean_record
  precip.mean.list3$Precip_loc_JuneJuly_mean_record = Precip_loc_JuneJuly_mean_record
  precip.mean.list3$Precip_loc_JulyAug_mean_record = Precip_loc_JulyAug_mean_record
  
}


results = list()
r2_training = list()
r2_test = list()


#  long stretch of months with linear model
combinations3 <- expand.grid(vpd = names(vpd.mean.list3), precip = names(precip.mean.list3))

#counter
count <- 0


# apply for each combination
run_model_3 <- function(vpd, precip) {
  count <<- count + 1
  
  lm_result <- lm_one_year_ahead(testing_index = 11:18,
                                 output_all = output_NDVI_all_mat_missing_filled[, Aug_index],
                                 input_1 = vpd.mean.list3[[vpd]][, 1:18],
                                 input_2 = precip.mean.list3[[precip]][, 1:18],
                                 input_3 = NA)
  
  
  r2_train_function <- function(i){
    1-(lm_result$sse_training_vector[i]/sum((output_NDVI_all_mat_missing_filled[,Aug_index][i,1:10]-mean(output_NDVI_all_mat_missing_filled[,Aug_index][i,1:10]))^2))
  }
  
  grand_r2 = mean(sapply(1:14000, r2_train_function))
  
  print(paste("Model run", count, "of", nrow(combinations3), ": R-squared Train =", grand_r2))
  
  return(c(r2_train = grand_r2))
}

#apply the function to each combination
results_3 <- mapply(run_model_3, combinations3$vpd, combinations3$precip, SIMPLIFY = FALSE)

results_df3 <- do.call(rbind, results_3)

results_df3 <- as.data.frame(results_df3)

colnames(results_df3) <- c("r2_train")

# results to matrix format for heatmap plotting
r2_train_matrix3 <- matrix(results_df3$r2_train, nrow = length(vpd.mean.list3), byrow = FALSE)

#dimension names for the matrices
dimnames(r2_train_matrix3) <- list(names(vpd.mean.list3), names(precip.mean.list3))

df_train_heatmap3 <- melt(r2_train_matrix3, varnames = c("VPD", "Precip"))

lm_R2_two_month <- df_train_heatmap3 %>% arrange(desc(value))

#Figure A1c
p_train_3 <- ggplot(df_train_heatmap3, aes(x = VPD, y = Precip, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(c(0,mean(df_train_heatmap$value))), limit = c(0, max(df_train_heatmap$value))) +
  scale_x_discrete(labels = c("Jan-Feb", "Feb-March", "March-April", "April-May", "May-June", "June-July", "July-Aug")) +
  scale_y_discrete(labels = c("Jan-Feb", "Feb-March", "March-April", "April-May", "May-June", "June-July", "July-Aug")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Monthly Ranges of VPD")+
  ylab("Monthly Ranges of Precipitation")+
  theme(axis.text.x = element_text( size =14, ),
        axis.text.y = element_text( size = 14),  
        axis.title.x = element_text( size = 18),
        axis.title.y = element_text( size = 18))+
  theme(legend.text = element_text(size = 18))+
  #theme(legend.title = element_text(size = 18))+
  theme(legend.box.just = "center")+
  theme(legend.key.height = unit(2, "cm"))+
  geom_text(aes(label = round(value, 2))) 


###  long stretch of months with linear model
combinations4 <- expand.grid(vpd = names(vpd.mean.list3), precip = names(precip.mean.list1))

#counter
count <- 0


# apply for each combination
run_model_4 <- function(vpd, precip) {
  count <<- count + 1
  
  lm_result <- lm_one_year_ahead(testing_index = 11:18,
                                 output_all = output_NDVI_all_mat_missing_filled[, Aug_index],
                                 input_1 = vpd.mean.list3[[vpd]][, 1:18],
                                 input_2 = precip.mean.list1[[precip]][, 1:18],
                                 input_3 = NA)
  
  
  r2_train_function <- function(i){
    1-(lm_result$sse_training_vector[i]/sum((output_NDVI_all_mat_missing_filled[,Aug_index][i,1:10]-mean(output_NDVI_all_mat_missing_filled[,Aug_index][i,1:10]))^2))
  }
  
  grand_r2 = mean(sapply(1:14000, r2_train_function))
  
  print(paste("Model run", count, "of", nrow(combinations4), ": R-squared Train =", grand_r2))
  
  return(c(r2_train = grand_r2))
}

#apply the function to each combination
results_4 <- mapply(run_model_4, combinations4$vpd, combinations4$precip, SIMPLIFY = FALSE)

results_df4 <- do.call(rbind, results_4)

results_df4 <- as.data.frame(results_df4)

colnames(results_df4) <- c("r2_train")

# results to matrix format for heatmap plotting
r2_train_matrix4 <- matrix(results_df4$r2_train, nrow = length(vpd.mean.list3), byrow = FALSE)

#dimension names for the matrices
dimnames(r2_train_matrix4) <- list(names(vpd.mean.list3), names(precip.mean.list1))

df_train_heatmap4 <- melt(r2_train_matrix4, varnames = c("VPD", "Precip"))

lm_R2_two_month_3_1 <- df_train_heatmap4 %>% arrange(desc(value))

#Figure A1d
p_train_4 <- ggplot(df_train_heatmap4, aes(x = VPD, y = Precip, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(c(0,mean(df_train_heatmap$value))), limit = c(0, max(df_train_heatmap$value))) +
  scale_x_discrete(labels = c("Jan-Feb", "Feb-March", "March-April", "April-May", "May-June", "June-July", "July-Aug")) +
  scale_y_discrete(labels = c("Jan-Aug", "Feb-Aug", "March-Aug", "April-Aug", "May-Aug", "June-Aug", "July-Aug")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Monthly Ranges of VPD")+
  ylab("Monthly Ranges of Precipitation")+
  theme(axis.text.x = element_text( size =14, ),
        axis.text.y = element_text( size = 14),  
        axis.title.x = element_text( size = 18),
        axis.title.y = element_text( size = 18))+
  theme(legend.text = element_text(size = 18))+
  theme(legend.box.just = "center")+
  theme(legend.key.height = unit(2, "cm"))+
  geom_text(aes(label = round(value, 2))) 
