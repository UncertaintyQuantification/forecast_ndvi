source('Functions/manuscript_functions.R')

packages = c("rlist","dplyr","tidyverse",
             "raster", "sf", "stringr", 
             "ggplot2", "ggpubr","RColorBrewer", "colorspace","MASS","plot3D","RobustGaSP","elevatr", "tseries") 
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

##5. Table S13, Precip forecasting

baseline_precipitation <-matrix(NA, nrow=14000, ncol=18)

#baseline/location mean method for precipitation
for (t in 1:18){
  for (i in 1:14000){
    baseline_precipitation[i,t] = mean(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[i,1:t], na.rm=TRUE)
  }
}

#basline/location mean rmse for precipitation
sqrt(mean((list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[,11:18]-baseline_precipitation[,10:17])^2))

#Previous Year method for precipitation
sqrt(mean((list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[,11:18]-list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[,10:17])^2))

#LM
#lm precip predict
precip_lin_predict <- lm_update_each_year_covariate(testing_index=11:18,output_all= list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record)
#rmse
sqrt(mean((list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[,11:18]-precip_lin_predict$pred_lm)^2))
#L95
mean(precip_lin_predict$interval_length_lm_model_record)
#P95
sum(precip_lin_predict$prop95_lm_model_record)/(14000*8)


#PPGP Method (non-empirical)

preci_janaug_grid = list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record


delta_preci_janaug_grid = preci_janaug_grid - rowMeans(preci_janaug_grid[, 1:10]) %*% matrix(1, ncol=num_yrs_all)


# parameter tuning for the ou process
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
set.seed(1234) # this is the default seed
sample_loc = sample(1:14000, size=10, replace = F)

rho_vec <- rep(NA, 8)
eta_vec <- rep(NA, 8)
dlm_pred_delta_preci_mat = matrix(NA, nrow=num_loc, ncol=8)
all_upper95_updated  <- matrix(NA, 14000, 8)
all_lower95_updated <- matrix(NA, 14000, 8)

for(train_idx in 1:8){
  train_yrs = train_idx + 7
  optimize = optim(par = c(-.5,.2), fn = preci_parameter_tuning, t = train_yrs,
                   method = 'L-BFGS-B', lower=c(-0.99,0), upper = c(0.99,Inf))  
  rho_vec[train_idx] = optimize$par[1]
  eta_vec[train_idx] = optimize$par[2]
  y_pred = preci_predict(param = optimize$par, t = train_yrs+2)
  dlm_pred_delta_preci_mat[,train_idx] = y_pred$pred
  all_upper95_updated[,train_idx] = y_pred$pred + 1.96 * sqrt(y_pred$var)
  all_lower95_updated[,train_idx] = y_pred$pred - 1.96 * sqrt(y_pred$var)
  print(optimize$par)
}

# get the sigma_hat_mle as a quadratic form
rho_vec
eta_vec
plot(rho_vec)

sqrt(mean((delta_preci_janaug_grid[, 11:18] - dlm_pred_delta_preci_mat)^2))
# 8.240735 with the 10 subsampling validation

#properly scaled precip
dlm_pred_preci_mat = dlm_pred_delta_preci_mat + rowMeans(preci_janaug_grid[, 1:10]) %*% matrix(1, ncol=8)
sqrt(mean((preci_janaug_grid[, 11:18] - dlm_pred_preci_mat)^2))

#lower and upper bounds of CI
lower95 = all_lower95_updated + rowMeans(preci_janaug_grid[, 1:10]) %*% matrix(1, ncol=8)
upper95 = all_upper95_updated + rowMeans(preci_janaug_grid[, 1:10]) %*% matrix(1, ncol=8)
interval_length_dlm_updated = upper95 - lower95
prop95_dlm_updated = (upper95 > preci_janaug_grid[, 11:18])*(lower95<dlm_pred_preci_mat)

#L95
mean(interval_length_dlm_updated)
#P95
sum(prop95_dlm_updated)/(14000*8)



##6. Figures

#Figure 2b
loc_1=5346
plot(2003:2020, list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[loc_1,],xlab='Year',ylab='Jan-Aug Precipitation (cm)', cex.main = 2.5, ylim= c(0,45), type = 'p', pch = 16,cex.lab=2.5, cex.axis=2.5, cex =1.5)
polygon( c(2013:2020,rev(2013:2020)),c(lower95[loc_1,],rev(upper95[loc_1,])),col = "lightblue", border = F)
lines(2013:2020,dlm_pred_preci_mat[loc_1,],type='l', col = 'darkblue' , lwd = 4)
lines(2013:2020, list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[loc_1,11:18],type='p', col="black", pch = 16, cex =1.5)
abline(v = 2013 , lty = 2, lwd =2)
legend("topleft", legend= c("PPGP"), col ='darkblue' , lty = 1, lwd =4 , cex=2.5)

#Figure 7a

#2020 true precip
admin.sf <- st_crop(admin, extent(precip))
plot_column_as_image_ndvi(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[,18],
                          admin.sf,
                          precip,
                          min_limit =0 ,
                          max_limit =153, color= 'Blues', legend_name = "Precip\n(cm)")

#Figure 7b
#2020 ppgp prediction of precip
admin.sf <- st_crop(admin, extent(precip))
plot_column_as_image_ndvi(dlm_pred_preci_mat[,8],
                          admin.sf,
                          precip,
                          min_limit =0 ,
                          max_limit =153, color= 'Blues', legend_name = "Precip\n(cm)")

#Figure 7c
#abs difference between truth and prediction
plot_column_as_image_ndvi(abs(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[,18]-dlm_pred_preci_mat[,8]),
                          admin.sf,
                          precip,
                          min_limit =0 ,
                          max_limit =153, color= 'Blues',legend_name = "Precip\n(cm)")
