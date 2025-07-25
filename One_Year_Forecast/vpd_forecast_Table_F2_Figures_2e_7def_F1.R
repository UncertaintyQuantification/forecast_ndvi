source('Functions/manuscript_functions.R')

packages = c("rlist","dplyr","tidyverse",
             "raster", "sf", "stringr", 
             "ggplot2", "ggpubr","RColorBrewer", "colorspace","MASS","plot3D","RobustGaSP","elevatr") 
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

##5. Table S14, VPD forecasting

baseline_vpd <-matrix(NA, nrow=14000, ncol=18)

#baseline for vpd
for (t in 1:18){
  for (i in 1:14000){
    baseline_vpd[i,t] = mean(list_eda_ndvi$VPD_loc_JulyAug_mean_record[i,1:t], na.rm=TRUE)
  }
}

#basline rmse for vpd
sqrt(mean((list_eda_ndvi$VPD_loc_JulyAug_mean_record[,11:18]-baseline_vpd[,10:17])^2))

#just using last year as prediction, vpd
sqrt(mean((list_eda_ndvi$VPD_loc_JulyAug_mean_record[,11:18]-list_eda_ndvi$VPD_loc_JulyAug_mean_record[,10:17])^2))

#lm
vpd_lin_predict <- lm_update_each_year_covariate(testing_index=11:18,output_all= list_eda_ndvi$VPD_loc_JulyAug_mean_record)
#rmse
sqrt(mean((list_eda_ndvi$VPD_loc_JulyAug_mean_record[,11:18]-vpd_lin_predict$pred_lm)^2))
#L95
mean(vpd_lin_predict$interval_length_lm_model_record)
#P95
sum(vpd_lin_predict$prop95_lm_model_record)/(14000*8)


#PPGP

vpd_julaug_grid = list_eda_ndvi$VPD_loc_JulyAug_mean_record
delta_vpd_julaug_grid = vpd_julaug_grid -  rowMeans(vpd_julaug_grid[, 1:10]) %*% matrix(1, ncol=num_yrs_all)

set.seed(1234)
sample_loc = sample(1:14000, size=10, replace = F)

training_index=as.numeric(1:10)
testing_index=as.numeric(11:18)
kernel_type='matern_5_2'

# matrix for nrmse storage
NRMSE_vpd_record=matrix(NA,20,15)

#standard deviation for VPD
sd_vpd=sd(list_eda_ndvi$VPD_loc_JulyAug_mean_record[,11:18])

#matrices to store metrics
pred_vpd_JulyAug_testing_est_par=matrix(NA,num_loc,length(testing_index))
lower95_vpd_JulyAug_testing_est_par=matrix(NA,num_loc,length(testing_index))
upper95_vpd_JulyAug_testing_est_par=matrix(NA,num_loc,length(testing_index))
l2_pred_error_rgasp_model_record=matrix(NA,num_loc,length(testing_index))
interval_length_rgasp_model_record=matrix(NA,num_loc,length(testing_index))
prop95_rgasp_model_record=matrix(NA,num_loc,length(testing_index))


param_ini_vpd=c(-2,-6)

for(i_test in testing_index){ 
  # optimize parameters for VPD using the optim function
  # The training data goes up to the index right before the current test index
  #validation set size is d
  m_vpd=optim(param_ini_vpd,one_step_pred_nrmse_vpd,method='BFGS',i_train_last=i_test-1,d=3)
  m_vpd
  
  # set the parameters obtained by optim as nugget and range param
  est_all_par_vpd=m_vpd$par

  # training indices up to the point right before the current test index
  training_index_aug=as.numeric(1:(i_test-1))
  
  # testing input is just the current test index
  testing_input=as.matrix(as.numeric(i_test))
  
  # train ppgasp model for vpd using parameters obtained in optim
  ppgasp.model.vpd=ppgasp(design=as.matrix(as.numeric(training_index_aug)),
                          response=t(list_eda_ndvi$VPD_loc_JulyAug_mean_record[,training_index_aug]),
                          nugget.est=F,nugget=exp(est_all_par_vpd[2]),range.par=1/exp(est_all_par_vpd[1]) )
  
  print(i_test)
  print(c(1/exp(est_all_par_vpd[1]), exp(est_all_par_vpd[2])))
  #predict for the test set using the model
  pred.model.vpd=predict(ppgasp.model.vpd,testing_input=testing_input )
  
  # save metrics
  pred_vpd_JulyAug_testing_est_par[,i_test-testing_index[1]+1]= pred.model.vpd$mean
  lower95_vpd_JulyAug_testing_est_par[,i_test-testing_index[1]+1]= pred.model.vpd$lower95
  upper95_vpd_JulyAug_testing_est_par[,i_test-testing_index[1]+1]= pred.model.vpd$upper95
  l2_pred_error_rgasp_model_record[,i_test-testing_index[1]+1] = (pred.model.vpd$mean-list_eda_ndvi$VPD_loc_JulyAug_mean_record[,i_test])^2
  interval_length_rgasp_model_record[,i_test-testing_index[1]+1] = pred.model.vpd$upper95 - pred.model.vpd$lower95
  prop95_rgasp_model_record[,i_test-testing_index[1]+1] = (pred.model.vpd$upper95>list_eda_ndvi$VPD_loc_JulyAug_mean_record[,i_test])*(pred.model.vpd$lower95<list_eda_ndvi$VPD_loc_JulyAug_mean_record[,i_test])
  
}

#predicted vpd using ppgp
pred_vpd_JulyAug_testing_est_par[which(pred_vpd_JulyAug_testing_est_par<0)]=0

#rmse
sqrt(mean((list_eda_ndvi$VPD_loc_JulyAug_mean_record[,11:18]-pred_vpd_JulyAug_testing_est_par)^2))
#L95
mean(interval_length_rgasp_model_record)
#P95
sum(prop95_rgasp_model_record)/(14000*8)


#Figure 2e
loc=5346
plot(2003:2020, list_eda_ndvi$VPD_loc_JulyAug_mean_record[loc,],xlab='Year',ylab='July-Aug VPD (hPa)', cex.main = 2.5, ylim= c(5,60), type = 'p', pch = 16,cex.lab=2.5, cex.axis=2.5, cex =1.5)
polygon( c(2013:2020,rev(2013:2020)),c(lower95_vpd_JulyAug_testing_est_par[loc,],rev(upper95_vpd_JulyAug_testing_est_par[loc,])),col = "orangered", border = F)
lines(2013:2020, pred_vpd_JulyAug_testing_est_par[loc,],type='l', col = 'orangered4' , lwd = 4)
lines(2013:2020, list_eda_ndvi$VPD_loc_JulyAug_mean_record[loc,11:18],type='p', col="black", pch = 16, cex=1.5)
abline(v = 2013 , lty = 2, lwd =2)
legend("topleft", legend= c("PPGP"), col ='orangered4' , lty = 1, lwd =4 , cex=2.5)


#Figure 7d
#true vpd 2020
admin.sf <- st_crop(admin, extent(vpdmax))
plot_column_as_image_ndvi(list_eda_ndvi$VPD_loc_JulyAug_mean_record[,18],
                          admin.sf,
                          vpdmax,
                          min_limit =0 ,
                          max_limit =72, color= 'OrRd',legend_name = "VPD\n(hPa)")

#Figure 7e
#forecast vpd 2020 ppgp
admin.sf <- st_crop(admin, extent(vpdmax))
plot_column_as_image_ndvi(pred_vpd_JulyAug_testing_est_par[,8],
                          admin.sf,
                          vpdmax,
                          min_limit =0 ,
                          max_limit =72, color= 'OrRd',legend_name = "VPD\n(hPa)")

#Figure 7f
#abs vpd difference
admin.sf <- st_crop(admin, extent(vpdmax))
plot_column_as_image_ndvi(abs(list_eda_ndvi$VPD_loc_JulyAug_mean_record[,18]-pred_vpd_JulyAug_testing_est_par[,8]),
                          admin.sf,
                          vpdmax,
                          min_limit =0 ,
                          max_limit =72, color= 'OrRd',legend_name = "VPD\n(hPa)")

#Figure S6

#rolling window vpd figure S6
start_year=6
rmse_last_year_pred=rep(NA,length(start_year:11))
rmse_average_pred=rep(NA,length(start_year:11))
rmse_ppgasp_pred=rep(NA,length(start_year:11))

for(i_year in start_year:11){
  print(i_year)
  
  testing_index=as.numeric(i_year:(i_year+7)) 
  
  rmse_last_year_pred[i_year-start_year+1]=sqrt(mean((list_eda_ndvi$VPD_loc_JulyAug_mean_record[,testing_index-1]-list_eda_ndvi$VPD_loc_JulyAug_mean_record[,testing_index])^2))
  
  pred_previous_years=matrix(NA,dim(list_eda_ndvi$VPD_loc_JulyAug_mean_record[,testing_index])[1],dim(list_eda_ndvi$VPD_loc_JulyAug_mean_record[,testing_index])[2] )
  for(i_previous in 1:length(testing_index)){
    pred_previous_years[,i_previous]=rowMeans(list_eda_ndvi$VPD_loc_JulyAug_mean_record[,1:(testing_index[i_previous]-1)])
  }
  rmse_average_pred[i_year-start_year+1]=sqrt(mean((pred_previous_years-list_eda_ndvi$VPD_loc_JulyAug_mean_record[,testing_index])^2))
  
  kernel_type='matern_5_2'
  alpha=c(1) #kernel power parameter for  power exp
  pred_vpd_JulyAug_testing_est_par_here=matrix(NA,num_loc,length(testing_index))
  
  for(i_test in testing_index){ ##i_test is the last training input
    
    m_vpd=optim(param_ini_vpd,one_step_pred_nrmse_vpd,method='BFGS',i_train_last=i_test-1,d=3,
                kernel_type= kernel_type,alpha=alpha)
    m_vpd
    
    # set the parameters obtained by optim as nugget and range param
    est_all_par_vpd=m_vpd$par
    
    # training indices up to the point right before the current test index
    training_index_aug=as.numeric(1:(i_test-1))
    
    # testing input is just the current test index
    testing_input=as.matrix(as.numeric(i_test))
    
    training_selected_index=training_index_aug
    
    #train model
    ppgasp.model.vpd=ppgasp(design=as.matrix(as.numeric(training_selected_index)),
                            response=t(list_eda_ndvi$VPD_loc_JulyAug_mean_record[,training_selected_index]),
                            nugget.est=F,nugget=exp(est_all_par_vpd[2]),range.par=1/exp(est_all_par_vpd[1]),
                            kernel_type=kernel_type,alpha=alpha)
    
    # predict for the test set using the model
    pred.model.vpd=predict(ppgasp.model.vpd,testing_input=testing_input )
    
    # save mean and lower and upper credible interval
    pred_vpd_JulyAug_testing_est_par_here[,i_test-testing_index[1]+1]= pred.model.vpd$mean
    
  }
  pred_vpd_JulyAug_testing_est_par_here[which(pred_vpd_JulyAug_testing_est_par_here<0)]=0
  rmse_ppgasp_pred[i_year-start_year+1]=sqrt(mean((list_eda_ndvi$VPD_loc_JulyAug_mean_record[,testing_index]-pred_vpd_JulyAug_testing_est_par_here)^2))
  
}

# save data in format easy for ggplot
last_year_method <- cbind.data.frame(rmse=rmse_last_year_pred, method = rep("Previous Yr", length(rmse_last_year_pred)), year = c(2007:2012))
ppgasp_method <- cbind.data.frame(rmse = rmse_ppgasp_pred, method = rep("PPGP", length(rmse_ppgasp_pred)), year = c(2007:2012))
average_method <- cbind.data.frame(rmse = rmse_average_pred, method = rep("Location Mean", length(rmse_average_pred)), year = c(2007:2012))
#code and results for fno can be found in One_Year_Forecast folder, in FNO_vpd_forecast.py 
fno_pred <- read.table('Data/rolling_window_FNO.txt')
fno_pred <- fno_pred[[1]]
fno_method <- cbind.data.frame(rmse = fno_pred, method = rep("FNO", length(fno_pred)), year = c(2007:2012))

rmse_data <- rbind.data.frame(last_year_method, ppgasp_method, average_method, fno_method)

rmse_plot <- ggplot(data = rmse_data, 
                    aes(x=year, y=rmse, color=method, shape=method))+
  geom_point(size=6)+ scale_y_continuous(limits = c(2.5, 4.5))+scale_color_manual(values = c("Previous Yr" = "darkred", "PPGP"= "darkgreen", "Location Mean"="purple", "FNO"="steelblue"))+
  scale_shape_manual(values = c( "Previous Yr" = 18,  "PPGP" = 15, "Location Mean" = 1, "FNO" = 17))

rolling_window_plot<- rmse_plot+ theme(legend.position="bottom", panel.background = element_blank(),axis.title = element_text(size = 20),axis.text=element_text(size=20), 
                 legend.text=element_text(size=20), legend.background = element_blank(), legend.title=element_text(size=20),plot.title = element_text(hjust = 0.5,size=20),
                 axis.line.x = element_line(colour = "black", size =1), axis.line.y = element_line(colour = "black", size = 1), legend.key = element_blank(),
                 axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(color = 'Methods:', shape = "Methods:",x = "Prediction Years", y = "RMSE" ) + ggtitle("Rolling 8 Year Window Prediction Performance of VPD") +
  scale_x_continuous(labels=c('2007-2015', '2008-2016', '2009-2017', '2010-2018', '2011-2019', '2012-2020'))


print(rolling_window_plot)
mean(rmse_last_year_pred)
mean(rmse_ppgasp_pred)
mean(rmse_average_pred)
mean(fno_pred)

