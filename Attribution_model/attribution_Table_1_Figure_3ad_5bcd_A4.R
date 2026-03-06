source('Functions/manuscript_functions_vs2.R')
packages = c("rlist","dplyr","tidyverse",
             "raster", "sf", "stringr", 
             "ggplot2", "ggpubr","RColorBrewer", "colorspace","MASS","plot3D","RobustGaSP","elevatr", 
             "lattice", "viridisLite", "latticeExtra", "grid","tseries", "ranger") 
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


##5 Table 1, Location Mean
baseline_NDVI <- matrix(NA, nrow=14000, ncol=18)


for (t in 1:18){
  for (i in 1:14000){
    baseline_NDVI[i,t] = mean(output_NDVI_all_mat_missing_filled[,Aug_index][i,1:t], na.rm=TRUE)
  }
}
#basline rmse for ndvi
sqrt(mean((output_NDVI_all_mat_missing_filled[,Aug_index][,11:18]-baseline_NDVI[,10:17])^2))


##6. Table 1 LM

lm_results <- lm_update_each_year_ndvi(testing_index=11:18,output_all=output_NDVI_all_mat_missing_filled[,Aug_index],
                                       input_1=list_eda_ndvi$VPD_loc_JulyAug_mean_record,
                                       input_2=list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record,
                                       pred_vpd_JulyAug_testing_est_par =list_eda_ndvi$VPD_loc_JulyAug_mean_record[,11:18],
                                       pred_precip = list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[,11:18])

#rmse
sqrt(mean((output_NDVI_all_mat_missing_filled[,Aug_index][,11:18]-lm_results$pred_lm)^2))
#L 95
mean(lm_results$interval_length_lm_model_record)
#P 95
sum(lm_results$prop95_lm_model_record)/(14000*8)

##7. Table 1, RF
NDVI <- output_NDVI_all_mat_missing_filled[, Aug_index]
VPD <- list_eda_ndvi$VPD_loc_JulyAug_mean_record
Preci <- list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record

n_cols <- ncol(NDVI)
n_pix <- nrow(NDVI)

pred_ndvi_attribution <- matrix(NA, nrow = n_pix, ncol = n_cols)
rmse_results_attribution <- rep(NA, n_cols)

for (t in 11:n_cols) {
  
  preds_t <- rep(NA, n_pix)
  
  for (i in 1:n_pix) {
    
    y_train <- NDVI[i, 1:(t-1)]
    vpd_train <- VPD[i, 1:(t-1)]
    pre_train <- Preci[i, 1:(t-1)]
    
    train_data <- data.frame(
      y = as.numeric(y_train),
      VPD = as.numeric(vpd_train),
      Precip = as.numeric(pre_train)
    )
    
    train_data <- na.omit(train_data)
    
    rf_model <- ranger(
      y ~ .,
      data = train_data
    )
    
    test_data <- data.frame(
      VPD = VPD[i, t],
      Precip = Preci[i, t]
    )
    
    preds_t[i] <- predict(rf_model, test_data)$predictions
  }
  
  pred_ndvi_attribution[, t] <- preds_t
  
  rmse <- sqrt(mean((NDVI[, t] - preds_t)^2, na.rm = TRUE))
  rmse_results_attribution[t] <- rmse
  
  cat(t, " : ", signif(rmse, 3), "\n") 
}

sqrt(mean((pred_ndvi_attribution[,11:18]-output_NDVI_all_mat_missing_filled[,Aug_index][,11:18])^2))

##8. Table 1, G-PPGP

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

set.seed(1234) #tried various seeds, all good
###subsampling some grids to predict
num_samples=500 ##number of subsamples
sample_grid_indices=sample(num_loc,size=num_samples) ###need to use a more principal way for estimation
#sample_grid_indices = 1 # to check
record_prediction_NDVI_attribution=matrix(NA,nrow=num_loc,ncol = num_yrs_all)
lower95_attribution = matrix(NA,nrow=num_loc,ncol = num_yrs_all)
upper95_attribution = matrix(NA,nrow=num_loc,ncol = num_yrs_all)

#G-PPGP for attribution
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
    
    pred_rgasp_i_grid_attribution=predict(m_rgasp_i_grid,testing_input_attribution)

    record_prediction_NDVI_attribution[i_grid,time_t]=pred_rgasp_i_grid_attribution$mean
    lower95_attribution[i_grid,time_t] = pred_rgasp_i_grid_attribution$lower95
    upper95_attribution[i_grid,time_t] = pred_rgasp_i_grid_attribution$upper95
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


##9. Attribution

#Figure 3a
loc = 5346
test_output <- record_prediction_NDVI_attribution[,11:num_yrs_all][loc,]
plot(2003:2020, output_NDVI_all_mat_missing_filled[,Aug_index][loc,],xlab='Year',ylab='August NDVI', cex.main = 2.5, ylim= c(0.1,0.35), type = 'p', pch = 16,cex.lab=2.5, cex.axis=2.5, cex = 1.5)
polygon( c(2013:2020,rev(2013:2020)),c(lower95_attribution[,11:num_yrs_all][loc,],rev(upper95_attribution[,11:num_yrs_all][loc,])),col = "darkseagreen2", border = F)
lines(2013:2020, test_output,type='l', col = 'darkgreen' , lwd = 4)
lines(2013:2020, output_NDVI_all_mat_missing_filled[,Aug_index][loc,11:18],type='p', col="black", pch = 16, cex =1.5)
abline(v = 2013 , lty = 2, lwd =2)
legend("topleft", legend= c("G-PPGP"), col ='darkgreen' , lty = 1, lwd =4 , cex=2.5)


# Figure 5 b-d
#the 2020(index 18 par$est):
#-4.394635 -4.424772 -4.111560
betas <- exp(c(-4.394635, -4.424772))
nugget <- exp(c(-4.111560))


#pick a location

#Fig 5b
#location = 11105

#Fig 5c
#location= 7340

#Fig 5d
#location = 4996

#Fig 3d
location= 5346

#floor and ceiling of covariates at location
vpd_lower <- floor(min(list_eda_ndvi$VPD_loc_JulyAug_mean_record[location,])) -5
vpd_upper <- ceiling(max(list_eda_ndvi$VPD_loc_JulyAug_mean_record[location,])) +5
precip_lower <- floor(min(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[location,])) -10
precip_upper <- ceiling(max(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[location,])) +10

vpd_grid <- seq(vpd_lower,vpd_upper, by = 0.5)
precip_grid <- seq(precip_lower,precip_upper, by = 0.5)
data_grid <- expand.grid(VPD=vpd_grid, Precipitation=precip_grid)
data_grid$NDVI <- NA


GP_fitting_one_step_grid_validation <-function(i){
  vpd_grid <- seq(vpd_lower,vpd_upper, by = 0.5)
  precip_grid <- seq(precip_lower,precip_upper, by = 0.5)
  data_grid <- expand.grid(VPD=vpd_grid, Precipitation=precip_grid)
  data_grid$NDVI <- NA
  
  
  testing_index=11:18
  output_all=output_NDVI_all_mat_missing_filled[,Aug_index]
  input_1=list_eda_ndvi$VPD_loc_JulyAug_mean_record
  input_2=list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record
  
  loc= location
  year = 18
  
  
  #ndvi output at training years at specific location
  output=as.vector(output_all[loc,1:(year-1)])
  
  #training data of first input (true vpd)
  input=cbind(input_1[loc,1:(year-1)])
  
  #training input of both vpd and precip
  input=cbind(input,input_2[loc,1:(year-1)])
  
  #predicted VPD at test location (pred_vpd_JulyAug_testing_est_par is only length of testing_index)
  testing_input= as.matrix(data_grid[i, 1:2])
  
  m_i_loc=rgasp(design=input,response =output,nugget.est=F,
                nugget=nugget,range.par=1/betas) #zero.mean="Yes"
  
  
  m_i_loc_pred=predict(m_i_loc,t(as.matrix(as.numeric(testing_input) )))
  
  
  NDVI_grid_output = c(m_i_loc_pred$mean, m_i_loc_pred$sd)
  NDVI_grid_output
  
  
}

# year 2020, possible grids 
ndvi_sapply <- sapply(1:dim(data_grid)[1],GP_fitting_one_step_grid_validation )
data_grid$NDVI <- ndvi_sapply[1,]
data_grid$SD <- ndvi_sapply[2,]

coul <- colorRampPalette(brewer.pal(9, "Greens"))(dim(data_grid)[1]) #have to do 9 cause there is only 9 in greens
coul_true <- coul[which.min(rowSums(abs(data_grid[,1:2]-c(list_eda_ndvi$VPD_loc_JulyAug_mean_record[location, 18], list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[location, 18] ))))]


plot_level <-levelplot(NDVI ~ VPD*Precipitation, data=data_grid  ,xlab=list("July-Aug VPD (hPa)", cex=2.5), ylab=list(" Jan-Aug Precipitation (cm)", cex = 2.5),
                       main="", col.regions=coul, colorkey = list(title=  list("NDVI", cex = 2.5),title.control = list(side = "top"),labels=list(cex=2.5)),
                       scales = list(x=list(cex=2.5),y=list(cex=2.5))) 


plot_level <- plot_level + xyplot(Precipitation ~ VPD, data = data_grid,
                                  panel = function(x, y, ...) {
                                    
                                    grid.circle(x = unit(list_eda_ndvi$VPD_loc_JulyAug_mean_record[location, 18], "native"),
                                                y = unit(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[location, 18], "native"),
                                                r = unit(0.5, "cm"),
                                                gp = gpar(fill = coul_true, lwd=3))
                                  })

print(plot_level)

##9, Suppement Figures A4 a-f

#NDVI 2019 true
admin.sf <- st_crop(admin, extent(modis))
plot_column_as_image_ndvi(output_NDVI_all_mat_missing_filled[,Aug_index][,17],
                          admin.sf,
                          modis,
                          min_limit =0 ,
                          max_limit =1, color= 'Greens', legend_name = "NDVI")

#2019 attribution prediction
plot_column_as_image_ndvi(record_prediction_NDVI_attribution[,17],
                          admin.sf,
                          modis,
                          min_limit =0 ,
                          max_limit =1, color= 'Greens', legend_name = "NDVI")
#2019 attribution abs difference

plot_column_as_image_ndvi(abs(output_NDVI_all_mat_missing_filled[,Aug_index][,17]-record_prediction_NDVI_attribution[,17]),
                          admin.sf,
                          modis,
                          min_limit =0 ,
                          max_limit =1, color= 'Greens', legend_name = "NDVI")


#NDVI 2020 true
admin.sf <- st_crop(admin, extent(modis))
plot_column_as_image_ndvi(output_NDVI_all_mat_missing_filled[,Aug_index][,18],
                          admin.sf,
                          modis,
                          min_limit =0 ,
                          max_limit =1, color= 'Greens', legend_name = "NDVI")

#2020 attribution prediction
plot_column_as_image_ndvi(record_prediction_NDVI_attribution[,18],
                          admin.sf,
                          modis,
                          min_limit =0 ,
                          max_limit =1, color= 'Greens', legend_name = "NDVI")

#2020 attribution abs diff
plot_column_as_image_ndvi(abs(output_NDVI_all_mat_missing_filled[,Aug_index][,18]-record_prediction_NDVI_attribution[,18]),
                          admin.sf,
                          modis,
                          min_limit =0 ,
                          max_limit =1, color= 'Greens', legend_name = "NDVI")

