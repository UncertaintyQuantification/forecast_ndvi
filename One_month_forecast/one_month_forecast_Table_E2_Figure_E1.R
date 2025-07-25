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

n_training = 12*10 ## how many month, total is 216
n_testing = dim(output_NDVI_all_mat_missing_filled)[2]-n_training


#### delta data generation
###use the training data to form the trend 
output_training = output_NDVI_all_mat_missing_filled[,1:n_training]
output_testing = output_NDVI_all_mat_missing_filled[,(n_training+1):(n_training+n_testing)]
period=12
num_period_training=floor(ncol(output_training)/period)
num_period_testing=floor(ncol(output_testing)/period)

output_training_baseline=matrix(NA, nrow=num_loc, ncol=period)

for(i_loc in 1:num_loc){
  for(i in 1:period){
    index_here=i+(0:(num_period_training-1))*period
    
    output_training_baseline[i_loc, i] = mean(output_training[i_loc,index_here])
    
  }  
}

plot(output_training[14,])
plot(output_training_baseline[14,])

output_training_delta = matrix(NA, nrow=num_loc, ncol=n_training)

for(i_loc in 1:num_loc){
  for(i in 1:period){
    index_here=i+(0:(num_period_training-1))*period
    output_training_delta[i_loc, index_here]=output_training[i_loc, index_here]-output_training_baseline[i_loc,i]
  }  
}

output_testing_delta = matrix(NA, nrow=num_loc, ncol=n_testing)


for(i_loc in 1:num_loc){
  for(i in 1:period){
    index_here=i+(0:(num_period_testing-1))*period
    output_testing_delta[i_loc, index_here]=output_testing[i_loc, index_here]-output_training_baseline[i_loc,i]
  }
  if( (num_period_testing*period)<ncol(output_testing)){
    diff=ncol(output_testing)-num_period_testing*period
    for(j in 1:diff){
      output_testing_delta[i_loc, num_period_testing*period+j]=output_testing[i_loc, num_period_testing*period+j]-output_training_baseline[i_loc,j]
    }
  }
  
}
sum(is.na(output_testing_delta))

output_delta = cbind(output_training_delta, output_testing_delta)

#baseline for just using historic aug ndvi
baseline_NDVI <- matrix(NA, nrow=14000, ncol=18)


for (t in 1:18){
  for (i in 1:num_loc){
    baseline_NDVI[i,t] = mean(output_delta[,Aug_index][i,1:t], na.rm=TRUE)
  }
}
#basline rmse for ndvi
sqrt(mean((output_delta[,Aug_index][,11:18]-baseline_NDVI[,10:17])^2))
# 0.03634604

#pred ndvi next step, just using time as input

set.seed(1234)
sample_loc = sample(1:14000, size=10, replace = F)
one_step_pred_nrmse_vpd<-function(param,i_train_last,d,kernel_type='matern_5_2',alpha=c(1)){
  
  #assign the input parameters to new name
  est_all_par_vpd=param
  
  ##predict each step
  pred_vpd_JulyAug_testing=matrix(NA,num_loc,1)
  lower95_vpd_JulyAug_testing=matrix(NA,num_loc,1)
  upper95_vpd_JulyAug_testing=matrix(NA,num_loc,1)
  
  #training indices up to the last training observation minus d (a lag)
  training_index_aug=as.numeric(1:(i_train_last-d))
  
  # testing input indices, which are the next d observations after the last training observation
  testing_input=as.matrix(as.numeric( (i_train_last-d+1):(i_train_last) ))
  
  # make model using ppgasp
  ppgasp.model.vpd=ppgasp(design=as.matrix(as.numeric(training_index_aug)),
                          response=t(output_delta[sample_loc,training_index_aug]),
                          nugget.est=F,nugget=exp(est_all_par_vpd[2]),range.par=1/exp(est_all_par_vpd[1]),
                          kernel_type=kernel_type,alpha=alpha )
  
  
  # predict vpd values
  pred.model.vpd=predict(ppgasp.model.vpd,testing_input=testing_input )
  
  #metrics
  sd_baseline=sd(output_delta[,testing_input])
  NRMSE_vpd=sqrt(mean((output_delta[sample_loc,testing_input]-t(pred.model.vpd$mean) )^2))/sd_baseline
  return(NRMSE_vpd)
  
}

#ppgp
training_index=as.numeric(1:120)
testing_index=as.numeric(121:216)
kernel_type='matern_5_2'
est_all_par_vpd=c(-4,-10)


param_ini_vpd=c(-2,-6)

pred_ppgp_one_step=matrix(NA,num_loc,n_testing)
lower95_ppgp_one_step=matrix(NA,num_loc,n_testing)
upper95_ppgp_one_step=matrix(NA,num_loc,n_testing)

for(i_test in testing_index){ 
  # optimize parameters for VPD using the optim function
  # The training data goes up to the index right before the current test index
  #validation set size is d
  m_vpd=optim(param_ini_vpd,one_step_pred_nrmse_vpd,method='BFGS',i_train_last=i_test-1,d=12)

  # set the parameters obtained by optim as nugget and range param
  est_all_par_vpd=m_vpd$par

  # training indices up to the point right before the current test index
  training_index_aug=as.numeric(1:(i_test-1))
  
  # testing input is just the current test index
  testing_input=as.matrix(as.numeric(i_test))
  
  # train ppgasp model for vpd using parameters obtained in optim
  ppgasp.model.vpd=ppgasp(design=as.matrix(as.numeric(training_index_aug)),
                          response=t(output_delta[,training_index_aug]),
                          nugget.est=F,nugget=exp(est_all_par_vpd[2]),range.par=1/exp(est_all_par_vpd[1]) )
  
  print(i_test)
  print(c(1/exp(est_all_par_vpd[1]), exp(est_all_par_vpd[2])))
  #predict for the test set using the model
  pred.model.vpd=predict(ppgasp.model.vpd,testing_input=testing_input )
  
  # save metrics
  pred_ppgp_one_step[,i_test-testing_index[1]+1]= pred.model.vpd$mean
  lower95_ppgp_one_step[,i_test-testing_index[1]+1]= pred.model.vpd$lower95
  upper95_ppgp_one_step[,i_test-testing_index[1]+1]= pred.model.vpd$upper95
}

#one step ahead ppgp 0.04177811 with all data. 
sqrt(mean( (pred_ppgp_one_step-output_testing_delta)^2)) ##PPGP just with time information, range parameter didn't update
# 0.04547828
sqrt(mean( (pred_ppgp_one_step[,Aug_index[1:8]]-output_testing_delta[,Aug_index[1:8]])^2)) ##PPGP just with time information, range parameter didn't update
# 0.03369699


#linear model, just august month
testing_index= 11:18
pred_lm=matrix(NA,num_loc,8)
interval_length_lm_model_record=matrix(NA,num_loc,8)
prop95_lm_model_record=matrix(NA,num_loc,8)

for(i_test in 11:18){
  
  for(i_loc in 1:num_loc){
    
    #want to predict month 128, which is august of first test year, index 11. 
    output_training_lin=output_delta[i_loc,1:(Aug_index[i_test]-1)]
    #vector of 14000 x 1. NDVI values of month 128
    output_testing_lin=output_delta[i_loc,(Aug_index[i_test]):(Aug_index[i_test])]
    num_loc=dim(output_NDVI_all_mat_missing_filled)[1]
    
    num_obs=length(output_training_lin)
    
    #1:127 for time input for linear model. 
    input_training_lm=as.numeric(1:num_obs)
    
    df_pixel <- data.frame(August_NDVI = output_training_lin,
                           time = input_training_lm)
    
    # using training years fit model
    model_pixel <- lm(August_NDVI ~ time, data = df_pixel)
    
    df_pixel_test <- data.frame(time = as.numeric(num_obs+1))
    
    # predict using the linear model on the test years
    predictions <- predict(model_pixel, newdata = df_pixel_test, interval='predict')
    
    #save results
    pred_lm[i_loc,i_test-testing_index[1]+1]=predictions[,1]
    interval_length_lm_model_record[i_loc,i_test-testing_index[1]+1] = predictions[,3]-predictions[,2]
    prop95_lm_model_record[i_loc,i_test-testing_index[1]+1] = (predictions[,3]>output_NDVI_all_mat_missing_filled[,Aug_index][i_loc,i_test])*(predictions[,2]<output_NDVI_all_mat_missing_filled[,Aug_index][i_loc,i_test])
    
  }
  print(i_test)
}

#aug month pred
sqrt(mean((output_delta[,Aug_index][,11:18]-pred_lm)^2))

#LM with all months
testing_index = 121:216

pred_lm_all_months = matrix(NA, num_loc, length(testing_index))
interval_length_lm_all_months = matrix(NA, num_loc, length(testing_index))
prop95_lm_all_months = matrix(NA, num_loc, length(testing_index))

for(i_test in testing_index){
  
  for(i_loc in 1:num_loc){
    
    output_training_lin = output_delta[i_loc, 1:(i_test - 1)]
    output_testing_true = output_delta[i_loc, i_test]
    
    num_obs = length(output_training_lin)
    input_training_lm = as.numeric(1:num_obs)
    
    df_pixel <- data.frame(NDVI = output_training_lin,
                           time = input_training_lm)
    
    model_pixel <- lm(NDVI ~ time, data = df_pixel)
    
    df_pixel_test <- data.frame(time = num_obs + 1)
    predictions <- predict(model_pixel, newdata = df_pixel_test, interval = 'predict')
    
    col_i = i_test - testing_index[1] + 1
    pred_lm_all_months[i_loc, col_i] = predictions[1]
    interval_length_lm_all_months[i_loc, col_i] = predictions[3] - predictions[2]
    prop95_lm_all_months[i_loc, col_i] = (predictions[2] < output_testing_true) & (predictions[3] > output_testing_true)
    print(paste(i_loc, i_test))
  }
}


#linear model monthly prediction
sqrt(mean( (pred_lm_all_months-output_testing_delta)^2)) 
# 0.04155294
sqrt(mean( (pred_lm_all_months[,Aug_index[1:8]]-output_testing_delta[,Aug_index[1:8]])^2)) 
# 0.03769784




## how to recover the prediction with seasonal effect
seasonal_effect = kronecker(matrix(1, nrow=1, ncol=8), output_training_baseline)
pred_ppgp_one_step_season = pred_ppgp_one_step + seasonal_effect
lower95_ppgp_one_step_season = lower95_ppgp_one_step + seasonal_effect
upper95_ppgp_one_step_season = upper95_ppgp_one_step + seasonal_effect


pred_lm_season = pred_lm_all_months + seasonal_effect


#Figure S5
set.seed(100)
#plot 6 random locations. 
testing_month=121:(121+95)
for(i in 1:6){
  plot_index=sample(1:num_loc,size=1)
  plot(testing_month,output_testing[plot_index,],type='l',xlab='Month',ylab='NDVI',ylim=c(min(lower95_ppgp_one_step_season[plot_index,],output_testing[plot_index,])
                                                                                          ,max(upper95_ppgp_one_step_season[plot_index,],output_testing[plot_index,])),
       cex.axis = 2, cex.lab = 1.5)
  polygon( c(testing_month,rev(testing_month)),c(lower95_ppgp_one_step_season[plot_index,],
                                                 rev(upper95_ppgp_one_step_season[plot_index,])),col = "grey80", border = F)
  lines(testing_month,output_testing[plot_index,],type='l',xlab='month',ylab='NDVI', lwd= 2)
  
  lines(testing_month,pred_ppgp_one_step_season[plot_index,],type='l',col='blue',lty=6, lwd = 2)
  lines(testing_month,pred_lm_season[plot_index,],type='l',col='red',lty=3, lwd = 2)
  
  # legend("topright", legend = c("True", "PPGP", "Linear Model", "95% CI"),
  #        col = c("black", "blue", "red", "grey80"),
  #        lty = c(1, 6, 3, NA),
  #        lwd = c(2, 2, 2, NA),
  #        pch = c(NA, NA, NA, 15),
  #        pt.cex = 2.5,
  #        cex = 2,
  #        x.intersp =1,
  #        y.intersp =1,
  #        bty = "n")

}
