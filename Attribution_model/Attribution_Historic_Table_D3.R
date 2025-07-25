## Table S10
source('Functions/manuscript_functions.R')

packages = c("rlist","dplyr","tidyverse",
             "raster", "sf", "stringr", 
             "ggplot2", "ggpubr","RColorBrewer", "colorspace","MASS","plot3D","RobustGaSP","elevatr", 
             "lattice", "viridisLite", "latticeExtra", "grid","tseries") 
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

## Data loading
##1. Loading data rasters 

# necessary values
desired_layers <- length(2003:2020) * 12 #number of months between 2003 and 2020
first_layer <- length(1895:2020)*12 - desired_layers+1 #first layer is 2003 January
admin <- st_read('Data/fourcorners_states.shp') 

# NDVI raster
modis <- stack('Data/modis_cube_processed.grd') 

# Precipitation raster
precip_long <- stack('Data/precip_reproj_trim_res.grd') 
precip <- precip_long[[first_layer:nlayers(precip_long)]] # keep only 2003-2020

#VPD raster
vpdmax_long <- stack('Data/vpdmax_reproj_trim_res.grd')
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

#validation function for gppgp
rmse_gppgp_validation<-function(param,VPD_samples,precip_samples,NDVI_output, d=2){ ###one can add more argument
  if(is.null(dim(VPD_samples))){
    VPD_samples = matrix(VPD_samples, nrow=1)
    precip_samples = matrix(precip_samples, nrow=1)
    NDVI_output = matrix(NDVI_output, nrow=1)
  }
  beta_samples=exp(param[-length(param)])
  nugget_samples=exp(param[length(param)])
  n_years = ncol(VPD_samples)
  sse_gross = 0
  for(i_loc in 1:dim(VPD_samples)[1]){
    input=cbind(as.vector(VPD_samples[i_loc,1:(n_years-d)]),as.vector(precip_samples[i_loc,1:(n_years-d)]))
    output=as.vector(NDVI_output[i_loc,1:(n_years-d)])
    
    m_rgasp_samples=rgasp(design=input,response=output,nugget.est =F,
                          range.par = 1/beta_samples,nugget=nugget_samples )
    
    testing_input=cbind(as.vector(VPD_samples[i_loc,(n_years-d+1):n_years]),as.vector(precip_samples[i_loc,(n_years-d+1):n_years]))
    pred_rgasp_i_grid=predict(m_rgasp_samples,testing_input)
    sse_loc = sum((NDVI_output[i_loc,(n_years-d+1):n_years] - pred_rgasp_i_grid$mean)^2)
    
    sse_gross=sse_gross+sse_loc
  }
  ##add prior
  return(sqrt(sse_gross/(dim(VPD_samples)[1]*d)))
}

## for the historical value
set.seed(1234)
###subsampling some grids to predict
num_samples=500 ##number of subsamples
sample_grid_indices=sample(num_loc,size=num_samples) 
record_prediction_NDVI_prev=matrix(NA,nrow=num_loc,ncol = num_yrs_all)


lower95_prev = matrix(NA,nrow=num_loc,ncol = num_yrs_all)
upper95_prev = matrix(NA,nrow=num_loc,ncol = num_yrs_all)



## historic covariates

for(time_t in 11:num_yrs_all){
  print(time_t)
  
  param_ini=c(-7, -7, -9)
  
  m_est=optim(param_ini,rmse_gppgp_validation,VPD_samples=list_eda_ndvi$VPD_loc_JulyAug_mean_record[sample_grid_indices,1:(time_t-2)],
              precip_samples=list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[sample_grid_indices,1:(time_t-2)],
              NDVI_output=output_NDVI_all_mat_missing_filled[,Aug_index][sample_grid_indices,2:(time_t-1)]
  ) ###if num of grids are large, do not use BFGS as the gradient is too large
  print(m_est$par)
  print(m_est$value)
  
  
  
  
  for(i_grid in 1:num_loc){
    if(i_grid%%1000==0){
      print(i_grid)
    }
    input=cbind(list_eda_ndvi$VPD_loc_JulyAug_mean_record[i_grid,1:(time_t-2)],
                list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[i_grid,1:(time_t-2)])
    
    output=output_NDVI_all_mat_missing_filled[,Aug_index][i_grid,2:(time_t-1)]
    m_rgasp_i_grid=rgasp(design=input,response=output,
                         range.par = 1/exp(m_est$par[-length(m_est$par)]), nugget=exp(m_est$par[length(m_est$par)]),nugget.est =F )
    
    ###previous year obs
    testing_input_prev=cbind(list_eda_ndvi$VPD_loc_JulyAug_mean_record[i_grid,(time_t-1)],
                             list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[i_grid,(time_t-1)]
    )
    
    
    pred_rgasp_i_grid_prev=predict(m_rgasp_i_grid,testing_input_prev)
    
    record_prediction_NDVI_prev[i_grid,time_t]=pred_rgasp_i_grid_prev$mean
    
    lower95_prev[i_grid,time_t] = pred_rgasp_i_grid_prev$lower95
    upper95_prev[i_grid,time_t] = pred_rgasp_i_grid_prev$upper95
    
  }
}


###by  G-PPGP historic data
sqrt(mean((record_prediction_NDVI_prev[,11:num_yrs_all]-output_NDVI_all_mat_missing_filled[,Aug_index][,11:18])^2))
# 0.0369637 with c(-7, -7, -9), size=500, seed=1234
mean((output_NDVI_all_mat_missing_filled[,Aug_index][,11:18] > lower95_prev[,11:18]) * (output_NDVI_all_mat_missing_filled[,Aug_index][,11:18] < upper95_prev[,11:18]) ) 
# 95% coverage 0.9319464
mean(upper95_prev[,11:18]-lower95_prev[,11:18])
# L95% 0.1365863


#LM method with only historic data
lm_update_each_year_ndvi_historic<-function(testing_index=11:18,output_all,
                                            input_1=NA,input_2=NA,
                                            pred_vpd_JulyAug_testing_est_par = NULL, 
                                            pred_precip = NULL){
  
  pred_lm=matrix(NA,num_loc,length(testing_index))
  interval_length_lm_model_record=matrix(NA,num_loc,length(testing_index))
  prop95_lm_model_record=matrix(NA,num_loc,length(testing_index))
  
  for(i_test in testing_index){
    
    for(i_loc in 1:num_loc){
      
      #ndvi output at training years at specific location
      output=as.vector(output_all[i_loc,2:(i_test-1)])
      
      
      #get training input of both vpd and precip
      input_vpd=cbind(input_1[i_loc,1:(i_test-2)])
      input_precip=cbind(input_2[i_loc,1:(i_test-2)])
      #input_ndvi=cbind(output_all[i_loc,1:(i_test-2)])
      
      #predicted VPD at test location (pred_vpd_JulyAug_testing_est_par is only length of testing_index)
      #testing_input_vpd = pred_vpd_JulyAug_testing_est_par[i_loc, i_test - min(testing_index) + 1]
      testing_input_vpd = input_1[i_loc, i_test-1]
      
      #predicted precip at test location (pred_precip is only length of testing_index)
      #testing_input_precip = pred_precip[i_loc, i_test - min(testing_index) + 1]
      testing_input_precip = input_2[i_loc, i_test-1]
      #testing_input_ndvi = output_all[i_loc, i_test-1]
      
      df_pixel <- data.frame(August_NDVI = output,
                             avg_vpd_Jul_Aug = input_vpd,
                             avg_precip_Jan_Aug = input_precip
      )
      
      # using training years fit model
      model_pixel <- lm(August_NDVI ~ avg_vpd_Jul_Aug + avg_precip_Jan_Aug, data = df_pixel)
      
      df_pixel_test <- data.frame(avg_vpd_Jul_Aug = testing_input_vpd,
                                  avg_precip_Jan_Aug = testing_input_precip
      )
      
      # predict using the linear model on the test years
      predictions <- predict(model_pixel, newdata = df_pixel_test, interval='predict')
      
      #save results
      pred_lm[i_loc,i_test-testing_index[1]+1]=predictions[,1]
      interval_length_lm_model_record[i_loc,i_test-testing_index[1]+1] = predictions[,3]-predictions[,2]
      prop95_lm_model_record[i_loc,i_test-testing_index[1]+1] = (predictions[,3]>output_all[i_loc,i_test])*(predictions[,2]<output_all[i_loc,i_test])
      
    }
  }
  
  return.list=list()
  return.list$pred_lm=pred_lm
  return.list$interval_length_lm_model_record =interval_length_lm_model_record
  return.list$prop95_lm_model_record =prop95_lm_model_record
  
  return(return.list)
}

lm_historic = lm_update_each_year_ndvi_historic(testing_index=11:18,
                                               output_all=output_NDVI_all_mat_missing_filled[,Aug_index],
                                               input_1=list_eda_ndvi$VPD_loc_JulyAug_mean_record,
                                               input_2=list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record,
                                               pred_vpd_JulyAug_testing_est_par = NULL, 
                                               pred_precip = NULL)
sqrt(mean((lm_historic$pred_lm-output_NDVI_all_mat_missing_filled[,Aug_index][,11:18])^2))
# 0.04254147
mean(lm_historic$prop95_lm_model_record) 
# 95% coverage 0.9416875
mean(lm_historic$interval_length_lm_model_record)
# L95% 0.1581126
