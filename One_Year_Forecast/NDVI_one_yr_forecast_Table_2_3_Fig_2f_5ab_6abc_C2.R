source('One_Year_Forecast/vpd_forecast_TableS14_Figs2e_7def_S6.R')
source('One_Year_Forecast/precip_forecast_Table_S13_Figs2b_7abc.R')
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


##5. Table 2, Location Mean
baseline_NDVI <- matrix(NA, nrow=14000, ncol=18)


for (t in 1:18){
  for (i in 1:14000){
    baseline_NDVI[i,t] = mean(output_NDVI_all_mat_missing_filled[,Aug_index][i,1:t], na.rm=TRUE)
  }
}
#basline rmse for ndvi
sqrt(mean((output_NDVI_all_mat_missing_filled[,Aug_index][,11:18]-baseline_NDVI[,10:17])^2))

##6. Table 2, Last Year
sqrt(mean((output_NDVI_all_mat_missing_filled[,Aug_index][,11:18]-output_NDVI_all_mat_missing_filled[,Aug_index][,10:17])^2))

##7. Table 2, LM

#lm vpd predict
vpd_lin_predict <- lm_update_each_year_covariate(testing_index=11:18,output_all= list_eda_ndvi$VPD_loc_JulyAug_mean_record)
#lm precip predict
precip_lin_predict <- lm_update_each_year_covariate(testing_index=11:18,output_all= list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record)
#lm ndvi predict
lm_results <- lm_update_each_year_ndvi(testing_index=11:18,output_all=output_NDVI_all_mat_missing_filled[,Aug_index],
                                       input_1=list_eda_ndvi$VPD_loc_JulyAug_mean_record,
                                       input_2=list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record,
                                       pred_vpd_JulyAug_testing_est_par =vpd_lin_predict$pred_lm_covariate,
                                       pred_precip = precip_lin_predict$pred_lm_covariate)
#rmse
sqrt(mean((output_NDVI_all_mat_missing_filled[,Aug_index][,11:18]-lm_results$pred_lm)^2))
#L95
mean(lm_results$interval_length_lm_model_record)
#P95
sum(lm_results$prop95_lm_model_record)/(14000*8)


##8. Table 2 AR(1), by hand, using kalman filter

aug_ndvi = output_NDVI_all_mat_missing_filled[,Aug_index]
delta_aug_ndvi = aug_ndvi - rowMeans(aug_ndvi[, 1:10]) %*% matrix(1, ncol=num_yrs_all)

# First we need to get the profile_likelihood
neg_prof_log_lik = function(param, y){
  # construct the KF to get the likelihood (this is not for the data just for the quadratic/matrix inversion)
  n_obs = dim(y)[1]
  n_yr = dim(y)[2]
  m_KF_grid = matrix(NA, nrow=n_obs, ncol=n_yr) # the given mean
  a_KF_grid = matrix(NA, nrow=n_obs, ncol=n_yr) # the updated mean 
  f_KF_grid = matrix(NA, nrow=n_obs, ncol=n_yr) # the projected mean for observation, but would be the same as 'a'
  C_KF_grid = matrix(NA, nrow=n_obs, ncol=n_yr) # the given covariance
  R_KF_grid = matrix(NA, nrow=n_obs, ncol=n_yr) # the updated covariance
  Q_KF_grid = matrix(NA, nrow=n_obs, ncol=n_yr) # the projected covariance -- R_t + V_t
  
  a_KF_grid[,1]=0 # set to be zero, to allow the innovation term
  R_KF_grid[,1]=1 # let us set this as just one
  f_KF_grid[,1]=a_KF_grid[,1] # cuz the coeff is just 1 for OU process
  #Q_KF_grid[,1]=R_KF_grid[,1]+param[2] # set V_t as just 1 for now 
  Q_KF_grid[,1]=R_KF_grid[,1]+0
  Q_KF_inv=1/Q_KF_grid[1,1]    ##1d ##solve(Q_KF[1])
  m_KF_grid[,1]=a_KF_grid[,1]+R_KF_grid[,1]*Q_KF_inv*(y[,1]-f_KF_grid[,1])
  C_KF_grid[,1]=R_KF_grid[,1]-R_KF_grid[,1]*Q_KF_inv*R_KF_grid[,1]
  
  rho_KF= param[1]
  # Now time to filter
  for(i in 2:n_yr){
    a_KF_grid[,i]=rho_KF*m_KF_grid[,i-1]
    R_KF_grid[,i]=rho_KF^2*C_KF_grid[,i-1]+(1-rho_KF^2)*R_KF_grid[,1]
    
    f_KF_grid[,i]=a_KF_grid[,i]
    Q_KF_grid[,i]=R_KF_grid[,i]+0
    
    Q_KF_inv=1/(Q_KF_grid[,i])  ##solve
    m_KF_grid[,i]=a_KF_grid[,i]+R_KF_grid[,i]*Q_KF_inv*(y[, i]-f_KF_grid[,i])
    C_KF_grid[,i]=R_KF_grid[,i]-R_KF_grid[,i]*Q_KF_inv*R_KF_grid[,i]
  }
  
  quad_term = sum(log(rowMeans(((y-f_KF_grid)^2)/Q_KF_grid)))
  det_term = sum(log(Q_KF_grid))
  neg_prof_log_lik = n_yr*quad_term/2 + det_term/2
  return(neg_prof_log_lik)
}

neg_prof_log_lik(c(0.5, 1), delta_aug_ndvi[,1:10])

rho_vec <- rep(NA, 8)

for(i_yr in 1:8){
  optimized = optim(par = c(-.5), fn = neg_prof_log_lik, y = delta_aug_ndvi[,1:(i_yr+9)],
                    method = 'L-BFGS-B', lower=c(-0.99), upper = c(0.99))
  rho_vec[i_yr] = optimized$par[1]
}

# get the sigma_hat_mle as a quadratic form
rho_vec
eta_vec<- rep(0, 8) #for AR(1) equivalence

num_grid = nrow(delta_aug_ndvi)
sigma_2_hat_mle = matrix(NA, nrow=num_grid, ncol=8)

for(i_grid in 1:num_grid){
  for(i_yr in 1:8){
    y = delta_aug_ndvi[i_grid,1:(i_yr+9)]
    rho = rho_vec[i_yr]
    n_y = length(y)
    m_KF <- rep(NA, n_y) 
    a_KF <- rep(NA, n_y)
    f_KF <- rep(NA, n_y) 
    C_KF <- rep(NA, n_y)
    R_KF <- rep(NA, n_y)
    Q_KF <- rep(NA, n_y)
    
    #rho_KF <- rho 
    a_KF[1] <- 0 # set to be zero, to allow the innovation term
    R_KF[1] <- 1 # assuming C0 is equal to the sigma_2 which is 1 for prof_lik
    f_KF[1] <- a_KF[1] # cuz the coeff is just 1 for OU process
    Q_KF[1] <- R_KF[1]  + eta_vec[i_yr]# set V_t as just 1 for now and param is the eta (nugget)
    Q_KF_inv <- 1/Q_KF[1]    ##1d ##solve(Q_KF[1])
    m_KF[1] <- a_KF[1] +R_KF[1]*Q_KF_inv*(y[1]- f_KF[1])
    C_KF[1] <- R_KF[1] -R_KF[1]*Q_KF_inv*R_KF[1]
    
    # Now time to filter
    for(i in 2:n_y){
      
      rho_KF_iterate <- rho
      a_KF[i] <- rho_KF_iterate*m_KF[i - 1]
      R_KF[i] <- rho_KF_iterate^2*C_KF[i - 1] + (1 -rho_KF_iterate^2)*R_KF[1]
      
      f_KF[i] <- a_KF[i]
      Q_KF[i] <- R_KF[i] + eta_vec[i_yr]
      
      Q_KF_inv <- 1/(Q_KF[i])
      m_KF[i] <- a_KF[i] + R_KF[i]*Q_KF_inv*(y[i] - f_KF[i])
      C_KF[i] <- R_KF[i] - R_KF[i]*Q_KF_inv*R_KF[i]
    }
    
    #sigma_2_hat_mle_loc = sum(((y-f_KF)^2)/Q_KF)/n_y
    sigma_2_hat_mle_loc = mean(((y-f_KF)^2)/Q_KF)
    sigma_2_hat_mle[i_grid, i_yr] = sigma_2_hat_mle_loc
    
  }
}
sum(is.na(sigma_2_hat_mle))

# now with updated eta and sigma we predict the KF
all_upper95_updated  <- matrix(NA, 14000, 8)
all_lower95_updated <- matrix(NA, 14000, 8)
dlm_pred_ndvi_mat = matrix(NA, nrow=num_grid, ncol=8)

for(i_yr in 1:8){
  m_KF_grid = matrix(NA, nrow=num_loc, ncol=i_yr+10) # the given mean
  a_KF_grid = matrix(NA, nrow=num_loc, ncol=i_yr+10) # the updated mean 
  f_KF_grid = matrix(NA, nrow=num_loc, ncol=i_yr+10) # the projected mean for observation, but would be the same as 'a'
  C_KF_grid = matrix(NA, nrow=num_loc, ncol=i_yr+10) # the given covariance
  R_KF_grid = matrix(NA, nrow=num_loc, ncol=i_yr+10) # the updated covariance
  Q_KF_grid = matrix(NA, nrow=num_loc, ncol=i_yr+10) # the projected covariance -- R_t + V_t
  
  rho_KF= rho_vec[i_yr]  
  a_KF_grid[,1]=0 # set to be zero, to allow the innovation term
  R_KF_grid[,1]=sigma_2_hat_mle[,i_yr] # let us set this as just one
  f_KF_grid[,1]=a_KF_grid[,1] # cuz the coeff is just 1 for OU process
  Q_KF_grid[,1]=R_KF_grid[,1]+eta_vec[i_yr]*sigma_2_hat_mle[,i_yr] # sigma0^2 = eta * sigma^2
  Q_KF_inv=1/Q_KF_grid[,1]    ##1d ##solve(Q_KF[1])
  m_KF_grid[,1]=a_KF_grid[,1]+R_KF_grid[,1]*Q_KF_inv*(delta_aug_ndvi[,1]-f_KF_grid[,1])
  C_KF_grid[,1]=R_KF_grid[,1]-R_KF_grid[,1]*Q_KF_inv*R_KF_grid[,1]
  
  # now time to filter
  for(i in 2:(i_yr+10)){
    a_KF_grid[,i]=rho_KF*m_KF_grid[,i-1]
    R_KF_grid[,i]=rho_KF^2*C_KF_grid[,i-1]+(1-rho_KF^2)*R_KF_grid[,1]
    
    f_KF_grid[,i]=a_KF_grid[,i]
    Q_KF_grid[,i]=R_KF_grid[,i]+eta_vec[i_yr]*sigma_2_hat_mle[,i_yr]
    
    Q_KF_inv=1/(Q_KF_grid[,i])  ##solve
    m_KF_grid[,i]=a_KF_grid[,i]+R_KF_grid[,i]*Q_KF_inv*(delta_aug_ndvi[,i]-f_KF_grid[,i])
    C_KF_grid[,i]=R_KF_grid[,i]-R_KF_grid[,i]*Q_KF_inv*R_KF_grid[,i]
  }
  dlm_pred_ndvi_mat[,i_yr] = a_KF_grid[,i_yr+10]
  all_upper95_updated[, i_yr] = f_KF_grid[,i_yr+10] + 1.96*sqrt(Q_KF_grid[,i_yr+10])
  all_lower95_updated[, i_yr] = f_KF_grid[,i_yr+10] - 1.96*sqrt(Q_KF_grid[,i_yr+10])
}

#rmse
sqrt(mean((delta_aug_ndvi[, 11:18] - dlm_pred_ndvi_mat)^2)) 

lower95 = all_lower95_updated + rowMeans(aug_ndvi[, 1:10]) %*% matrix(1, ncol=8)
upper95 = all_upper95_updated + rowMeans(aug_ndvi[, 1:10]) %*% matrix(1, ncol=8)
l2_pred_error_dlm_updated = (dlm_pred_ndvi_mat-delta_aug_ndvi[,11:18])^2
interval_length_dlm_updated = all_upper95_updated - all_lower95_updated
prop95_dlm_updated = (all_upper95_updated > delta_aug_ndvi[, 11:18])*(all_lower95_updated<delta_aug_ndvi[, 11:18])

#L95
mean(interval_length_dlm_updated)
#P95
sum(prop95_dlm_updated)/(14000*8)

#ar prediction scaled back
ar_ndvi <- dlm_pred_ndvi_mat+ rowMeans(aug_ndvi[, 1:10]) %*% matrix(1, ncol=8)


##9. G-PPGP, attribution and one-year-ahead

##this is for matern 2.5, a function to add up the log marginal post
##using TWO covariates in the past year, VPD_samples,precip_samples without past NDVI

#rmse validation function for hyperparameter optimization
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
  return(sqrt(sse_gross/(dim(VPD_samples)[1]*d)))
}

# example
rmse_gppgp_validation(param=c(.1,.1,.1),
                      VPD_samples=list_eda_ndvi$VPD_loc_JulyAug_mean_record[1:10,1:10],
                      precip_samples=list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[1:10,1:10],
                      NDVI_output=output_NDVI_all_mat_missing_filled[,Aug_index][1:10,1:10]
)

set.seed(1234) #tried various seeds
###subsampling some grids to predict
num_samples=500 ##number of subsamples
sample_grid_indices=sample(num_loc,size=num_samples)
record_prediction_NDVI_attribution=matrix(NA,nrow=num_loc,ncol = num_yrs_all)
record_prediction_NDVI_two_stage=matrix(NA,nrow=num_loc,ncol = num_yrs_all)
lower95_attribution = matrix(NA,nrow=num_loc,ncol = num_yrs_all)
upper95_attribution = matrix(NA,nrow=num_loc,ncol = num_yrs_all)
lower95_two_stage = matrix(NA,nrow=num_loc,ncol = num_yrs_all)
upper95_two_stage = matrix(NA,nrow=num_loc,ncol = num_yrs_all)


#G-PPGP for both attribution and two stage model
for(time_t in 11:num_yrs_all){
  print(time_t)
  
  param_ini=c(-7, -7, -9)
  
  m_est=optim(param_ini,rmse_gppgp_validation,VPD_samples=list_eda_ndvi$VPD_loc_JulyAug_mean_record[sample_grid_indices,1:(time_t-1)],
              precip_samples=list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[sample_grid_indices,1:(time_t-1)],
              NDVI_output=output_NDVI_all_mat_missing_filled[,Aug_index][sample_grid_indices,1:(time_t-1)]
  ) ###if num of grids are large, do not use BFGS as the gradient is too large
  print(m_est$par)
  print(m_est$value)
  
  
  
  
  for(i_grid in 1:num_loc){
    if(i_grid%%1000==0){
      print(i_grid)
    }
    input=cbind(list_eda_ndvi$VPD_loc_JulyAug_mean_record[i_grid,1:(time_t-1)],
                list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[i_grid,1:(time_t-1)])
    
    output=output_NDVI_all_mat_missing_filled[,Aug_index][i_grid,1:(time_t-1)]
    m_rgasp_i_grid=rgasp(design=input,response=output,
                         range.par = 1/exp(m_est$par[-length(m_est$par)]), nugget=exp(m_est$par[length(m_est$par)]),nugget.est =F )
    
    ###previous year obs
    testing_input_attribution=cbind(list_eda_ndvi$VPD_loc_JulyAug_mean_record[i_grid,(time_t)],
                        list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[i_grid,(time_t)]
    )
    
    testing_input_two_stage=cbind(
      pred_vpd_JulyAug_testing_est_par[i_grid,(time_t-10)],
      dlm_pred_preci_mat[i_grid,(time_t-10)]
    )
    
    pred_rgasp_i_grid_attribution=predict(m_rgasp_i_grid,testing_input_attribution)
    pred_rgasp_i_grid_two_stage=predict(m_rgasp_i_grid,testing_input_two_stage)
    
    record_prediction_NDVI_attribution[i_grid,time_t]=pred_rgasp_i_grid_attribution$mean
    record_prediction_NDVI_two_stage[i_grid,time_t]=pred_rgasp_i_grid_two_stage$mean
    lower95_attribution[i_grid,time_t] = pred_rgasp_i_grid_attribution$lower95
    upper95_attribution[i_grid,time_t] = pred_rgasp_i_grid_attribution$upper95
    lower95_two_stage[i_grid,time_t] = pred_rgasp_i_grid_two_stage$lower95
    upper95_two_stage[i_grid,time_t] = pred_rgasp_i_grid_two_stage$upper95
  }
}


#by G-PPGP, attribution  
sqrt(mean((record_prediction_NDVI_attribution[,11:num_yrs_all]-output_NDVI_all_mat_missing_filled[,Aug_index][,11:18])^2))
# 0.02911757 with c(-7, -7, -9), size=500, seed=1234
#0.02912207 with seed = 0 
#0.02915573 with seed =1998

mean((output_NDVI_all_mat_missing_filled[,Aug_index][,11:18] > lower95_attribution[,11:18]) * (output_NDVI_all_mat_missing_filled[,Aug_index][,11:18] < upper95_attribution[,11:18]) ) 
# 95% coverage 0.9492768

mean(upper95_attribution[,11:18]-lower95_attribution[,11:18])
# L95% 0.1098873

#by G-PPGP, one-year-ahead
sqrt(mean((record_prediction_NDVI_two_stage[,11:num_yrs_all]-output_NDVI_all_mat_missing_filled[,Aug_index][,11:18])^2))
# 0.03517
mean((output_NDVI_all_mat_missing_filled[,Aug_index][,11:18] > lower95_two_stage[,11:18]) * (output_NDVI_all_mat_missing_filled[,Aug_index][,11:18] < upper95_two_stage[,11:18]) ) 
# 95% coverage 0.883
mean(upper95_two_stage[,11:18]-lower95_two_stage[,11:18])
# L95% 0.102843


##10. Figure 2f
loc =5346
test_output <- record_prediction_NDVI_two_stage[,11:num_yrs_all][loc,]
plot(2003:2020, output_NDVI_all_mat_missing_filled[,Aug_index][loc,],xlab='Year',ylab='August NDVI', cex.main = 2.5, ylim= c(0.15,0.3), type = 'p', pch = 16,cex.lab=2.5, cex.axis=2.5, cex =1.5)
polygon( c(2013:2020,rev(2013:2020)),c(lower95_two_stage[,11:num_yrs_all][loc,],rev(upper95_two_stage[,11:num_yrs_all][loc,])),col = "darkseagreen2", border = F)
lines(2013:2020, test_output,type='l', col = 'darkgreen' , lwd = 4)
lines(2013:2020, output_NDVI_all_mat_missing_filled[,Aug_index][loc,11:18],type='p', col="black", pch = 16, cex =1.5)
abline(v = 2013 , lty = 2, lwd =2)
legend("topleft", legend= c("G-PPGP"), col ='darkgreen' , lty = 1, lwd =4 , cex=2.5)




## 11, Figure 6b-c (6a is in data figure file)

#Figure 6b, predicted 2020 NDVI
admin.sf <- st_crop(admin, extent(modis))

plot_column_as_image_ndvi(record_prediction_NDVI_two_stage[,18],
                          admin.sf,
                          modis,
                          min_limit =0 ,
                          max_limit =1, color= 'Greens', legend_name = "NDVI")

#Figure 6c
#read in fno output
fno <- read.matrix('Data/pred_reshaped.txt')

loc = 5346
FNO_reconfigure <- cbind(as.vector(t(fno[1:140,])),as.vector(t(fno[141:280,])),
                         as.vector(t(fno[281:420,])), as.vector(t(fno[421:560,])),
                         as.vector(t(fno[561:700,])), as.vector(t(fno[701:840,])),
                         as.vector(t(fno[841:980,])),
                         as.vector(t(fno[981:1120,])))

test_output <- record_prediction_NDVI_two_stage[,11:18][loc,]
par(mar=c(5,8,6,2))
plot(2003:2020, output_NDVI_all_mat_missing_filled[,Aug_index][loc,],xlab='Year',ylab='August NDVI', cex.main =3, ylim= c(0.15,0.3), type = 'p', pch = 16,cex.lab=2, cex.axis=2, cex = 1.5)
polygon( c(2013:2020,rev(2013:2020)),c(lower95_two_stage[,11:18][loc,],rev(upper95_two_stage[,11:18][loc,])),col = "darkseagreen2", border = F)
lines(2013:2020, test_output,type='l', col = 'darkgreen' , lwd = 4)
lines(2013:2020, test_output,type='l', col = 'darkgreen' , lwd = 4)
lines(2013:2020, output_NDVI_all_mat_missing_filled[,Aug_index][loc,11:18],type='p', col="black", pch = 16, cex = 1.5)
abline(v = 2013 , lty = 2)
legend("topleft", legend= c("G-PPGP", "FNO", "AR(1)"), col =c('darkgreen','steelblue', "deeppink") , lty = c(1,2,3), lwd =c(4,4,4) , cex=2)
lines(2013:2020, FNO_reconfigure[loc,],lty=2, col = 'steelblue' , lwd = 4)
lines(2013:2020, ar_ndvi[loc,],lty=3, col = 'deeppink' , lwd = 4)


##12. Figure 5. Gross NDVI. 
RNN <- as.matrix(read.table('Data/all_predictions_RNN.txt'))

Years <- 2013:2020
real_test_data <- colSums(output_NDVI_all_mat_missing_filled[,Aug_index][,11:18])
FNO_data <- c(sum(as.vector(t(fno[1:140,]))),sum(as.vector(t(fno[141:280,]))),
              sum(as.vector(t(fno[281:420,]))), sum(as.vector(t(fno[421:560,]))),
              sum(as.vector(t(fno[561:700,]))), sum(as.vector(t(fno[701:840,]))),
              sum(as.vector(t(fno[841:980,]))),
              sum(as.vector(t(fno[981:1120,]))))
Last_year_data <- colSums(output_NDVI_all_mat_missing_filled[,Aug_index][,10:17])
PPGP_data <- colSums(record_prediction_NDVI_two_stage[,11:18])
LM_data <- colSums(lm_results$pred_lm)
baseline_data<- colSums(baseline_NDVI)[11:18]
AR_Data <- colSums(ar_ndvi)
RNN_Data <- colSums(RNN)

all_data <- data.frame(
  years = Years,
  real = real_test_data,
  last_year = Last_year_data,
  ppgp = PPGP_data,
  fno_test = FNO_data,
  lin_data =LM_data,
  average = baseline_data,
  ar = AR_Data,
  rnn_dat = RNN_Data
)


#Figure 5a
ggplot(all_data, aes(x=years)) + 
  geom_line(aes(y = real, color = "Truth", linetype ="Truth"), size = 1.75) + 
  geom_line(aes(y = fno_test, color="FNO",linetype = "FNO"), size = 1.75) +
  geom_line(aes(y = ppgp, color="G-PPGP", linetype = "G-PPGP"), size = 1.75) +
  #geom_line(aes(y = last_year, color="Previous Yr", linetype = "Previous Yr"), size = 1.75)+
  #geom_line(aes(y = lin_data, color="Linear", linetype = "Linear"), size = 1.75)+
  geom_line(aes(y = average, color="Location Mean", linetype = "Location Mean"), size = 1.75)+
  #geom_ribbon(aes(ymin = as.vector(lower_95_gppgp[,11:18]), ymax = as.vector(upper_95_gppgp[,11:18]) ), alpha = 0.5) +
  geom_line(aes(y=ar, color="AR(1)", linetype ="AR(1)"), size = 1.75)+
  #geom_line(aes(y=rnn_dat, color="RNN", linetype ="RNN"), size = 1.75)+
  
  geom_point(aes(y = real, color = "Truth",  shape = "Truth"), size = 4)+
  geom_point(aes(y = fno_test, color = "FNO", shape = "FNO"), size = 4)+
  geom_point(aes(y = ppgp, color = "G-PPGP",  shape = "G-PPGP"), size = 4)+
  #geom_point(aes(y = last_year, color = "Previous Yr",  shape = "Previous Yr"), size = 4)+
  #geom_point(aes(y = lin_data, color = "Linear", shape = "Linear"), size = 4)+
  geom_point(aes(y = average, color = "Location Mean", shape = "Location Mean"), size = 4)+
  geom_point(aes(y = ar, color = "AR(1)", shape = "AR(1)"), size = 4)+
  #geom_point(aes(y=rnn_dat, color="RNN", shape ="RNN"), size = 4)+
  
  labs(title = "One Year Prediction of Gross NDVI",  y = "Gross NDVI", x = "Years", color = "Methods", linetype = "Methods")+
  scale_color_manual(name = "Methods", values = c("Truth" = "black", "FNO" = "steelblue", "G-PPGP" = "darkgreen", "Previous Yr" = "darkred", "Linear" = "coral1", "Location Mean" = "purple", "AR(1)"="deeppink", "RNN" = "saddlebrown"))+
  scale_linetype_manual(name = "Methods",values = c("Truth" = "solid", "FNO" = "dashed", "G-PPGP" = "twodash", "Previous Yr" = "longdash", "Linear" = "4C88C488", "Location Mean" = "dotdash", "AR(1)"="dotted", "RNN" ="324B644B"))+
  scale_shape_manual(name = "Methods",values = c("Truth" = 16, "FNO" = 17, "G-PPGP" = 15, "Previous Yr" = 18, "Linear" = 0, "Location Mean" = 1, "AR(1)" = 6, "RNN" = 5))+
  theme( text = element_text(size = 22),axis.text = element_text(size = 22),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.key=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))


#residuals
FNO_data_residuals <- real_test_data - FNO_data
Last_year_data_residuals <- real_test_data - Last_year_data
PPGP_data_residuals <- real_test_data - PPGP_data
LM_data_residuals <- real_test_data - LM_data 
baseline_data_residuals <- real_test_data - baseline_data
ar_residuals <- real_test_data - AR_Data
rnn_residuals <- real_test_data - RNN_Data

all_data_residuals <- data.frame(
  years_res = Years,
  last_year_res = Last_year_data_residuals,
  ppgp_res = PPGP_data_residuals,
  fno_test_res = FNO_data_residuals,
  lin_data_res =LM_data_residuals,
  average_res = baseline_data_residuals,
  AR_res = ar_residuals,
  RNN_res = rnn_residuals
)


colnames(all_data_residuals) <- c("Years", "Previous Yr", "G-PPGP", 
                                  "FNO","Linear", "Location Mean", "AR(1)", "RNN") 

all_data_residuals <-all_data_residuals %>%
  pivot_longer(!Years, names_to = "Method", values_to = "Residual")

#Table 2
as.data.frame(all_data_residuals)

sqrt(mean((all_data$ppgp-all_data$real)^2)) #rmse gppgp
sqrt(mean((all_data$last_year-all_data$real)^2)) #previous year
sqrt(mean((all_data$fno_test-all_data$real)^2)) #fno
sqrt(mean((all_data$lin_data-all_data$real)^2)) #LM
sqrt(mean((all_data$average-all_data$real)^2)) #location mean
sqrt(mean((all_data$ar-all_data$real)^2)) #ar
sqrt(mean((all_data$rnn_dat-all_data$real)^2)) #RNN

all_data_residuals_top4 <- filter(all_data_residuals,all_data_residuals$Method %in% c("G-PPGP", "Location Mean", "FNO", "AR(1)"))
#Figure 5b
vertical_lines = (Years + 0.5)
ggplot(all_data_residuals_top4, aes(fill=Method, y= Residual, x= Years)) + 
  geom_bar(position="dodge", stat="identity")+
  geom_vline(xintercept = vertical_lines, linetype ="dashed", lwd= 0.75 ) +
  geom_hline(yintercept = c(0))+
  scale_fill_manual(name = "Methods", values = c( "FNO" = "steelblue", "G-PPGP" = "darkgreen", "Previous Yr" = "darkred", "Linear" = "coral1", "Location Mean" = "purple", "AR(1)" = "deeppink", "RNN" = "saddlebrown"))+
  ylim(c(-625, 550))+
  theme( text = element_text(size = 22),axis.text = element_text(size = 22),
         panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"),legend.key=element_blank(), panel.grid = element_line(color = "lightgrey"))+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "Gross NDVI Residuals")


## Supplement Gross NDVI Figure S4a
ggplot(all_data, aes(x=years)) + 
  geom_line(aes(y = real, color = "Truth", linetype ="Truth"), size = 1.75) + 
  geom_line(aes(y = fno_test, color="FNO",linetype = "FNO"), size = 1.75) +
  geom_line(aes(y = ppgp, color="G-PPGP", linetype = "G-PPGP"), size = 1.75) +
  geom_line(aes(y = last_year, color="Previous Yr", linetype = "Previous Yr"), size = 1.75)+
  geom_line(aes(y = lin_data, color="Linear", linetype = "Linear"), size = 1.75)+
  geom_line(aes(y = average, color="Location Mean", linetype = "Location Mean"), size = 1.75)+
  #geom_ribbon(aes(ymin = as.vector(lower_95_gppgp[,11:18]), ymax = as.vector(upper_95_gppgp[,11:18]) ), alpha = 0.5) +
  geom_line(aes(y=ar, color="AR(1)", linetype ="AR(1)"), size = 1.75)+
  geom_line(aes(y=rnn_dat, color="RNN", linetype ="RNN"), size = 1.75)+
  
  geom_point(aes(y = real, color = "Truth",  shape = "Truth"), size = 4)+
  geom_point(aes(y = fno_test, color = "FNO", shape = "FNO"), size = 4)+
  geom_point(aes(y = ppgp, color = "G-PPGP",  shape = "G-PPGP"), size = 4)+
  geom_point(aes(y = last_year, color = "Previous Yr",  shape = "Previous Yr"), size = 4)+
  geom_point(aes(y = lin_data, color = "Linear", shape = "Linear"), size = 4)+
  geom_point(aes(y = average, color = "Location Mean", shape = "Location Mean"), size = 4)+
  geom_point(aes(y = ar, color = "AR(1)", shape = "AR(1)"), size = 4)+
  geom_point(aes(y=rnn_dat, color="RNN", shape ="RNN"), size = 4)+
  
  labs(title = "One Year Prediction of Gross NDVI",  y = "Gross NDVI", x = "Years", color = "Methods", linetype = "Methods")+
  scale_color_manual(name = "Methods", values = c("Truth" = "black", "FNO" = "steelblue", "G-PPGP" = "darkgreen", "Previous Yr" = "darkred", "Linear" = "coral1", "Location Mean" = "purple", "AR(1)"="deeppink", "RNN" = "saddlebrown"))+
  scale_linetype_manual(name = "Methods",values = c("Truth" = "solid", "FNO" = "dashed", "G-PPGP" = "twodash", "Previous Yr" = "longdash", "Linear" = "4C88C488", "Location Mean" = "dotdash", "AR(1)"="dotted", "RNN" ="324B644B"))+
  scale_shape_manual(name = "Methods",values = c("Truth" = 16, "FNO" = 17, "G-PPGP" = 15, "Previous Yr" = 18, "Linear" = 0, "Location Mean" = 1, "AR(1)" = 6, "RNN" = 5))+
  theme( text = element_text(size = 22),axis.text = element_text(size = 22),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.key=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))

## SUpplement Residual Figure S4b
vertical_lines = (Years + 0.5)
ggplot(all_data_residuals, aes(fill=Method, y= Residual, x= Years)) + 
  geom_bar(position="dodge", stat="identity")+
  geom_vline(xintercept = vertical_lines, linetype ="dashed", lwd= 0.75 ) +
  geom_hline(yintercept = c(0))+
  scale_fill_manual(name = "Methods", values = c( "FNO" = "steelblue", "G-PPGP" = "darkgreen", "Previous Yr" = "darkred", "Linear" = "coral1", "Location Mean" = "purple", "AR(1)" = "deeppink", "RNN" = "saddlebrown"))+
  ylim(c(-625, 550))+
  theme( text = element_text(size = 22),axis.text = element_text(size = 22),
         panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"),legend.key=element_blank(), panel.grid = element_line(color = "lightgrey"))+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "Gross NDVI Residuals")