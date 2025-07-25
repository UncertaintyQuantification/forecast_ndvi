# This file contains functions that are called in the other R files to produce manuscript results and figures. 


## function that gets neighboring indices 
get_ne<-function(y){ 
  
  # get the interior square so we will have the same number of 9 index
  # these points are on the edge of the dataset and do not have 8 neighbors
  loc_exterior=c(1:dim(y)[1],(1:(dim(y)[2]-2))*dim(y)[1]+c(1),(1:(dim(y)[2]-2))*dim(y)[1]+dim(y)[1],(dim(y)[2]-1)*dim(y)[1]+1:dim(y)[1])
  loc_interior=(1:(dim(y)[1]*dim(y)[2]))[-loc_exterior]
  
  # 9 neighbor of the interior
  # calculate the indices of the 9 neighbors (including the point itself) for each interior point in the dataset
  loc_interior_ne=matrix(NA,length(loc_interior),9)
  for(loc_i in 1:length(loc_interior)){
    loc_interior_ne[loc_i,]=c(loc_interior[loc_i]-dim(y)[1]-1,loc_interior[loc_i]-dim(y)[1],loc_interior[loc_i]-dim(y)[1]+1,
                              loc_interior[loc_i]-1,loc_interior[loc_i],loc_interior[loc_i]+1,
                              loc_interior[loc_i]+dim(y)[1]-1,loc_interior[loc_i]+dim(y)[1],loc_interior[loc_i]+dim(y)[1]+1)
  }
  
  ##all neighbors just for smoothing
  #calculate the indices of the neighbors for the exterior points. 
  #exterior points do not have 8 neighbors so different rules are used depending on whether the point is a corner point or an edge point.
  loc_ne=matrix(NA,dim(y)[1]*dim(y)[2],9)
  loc_ne[loc_interior,]=loc_interior_ne
  
  ##neighbor of the 4 border
  for(i in 2:(dim(y)[1]-1) ){
    loc_ne[i,1:6]=c(i-1,i,i+1,i-1+dim(y)[1],i+dim(y)[1],i+1+dim(y)[1] )
  }
  upper_level_index=(1:(dim(y)[2]-2))*dim(y)[1]+1
  
  for(i in upper_level_index ){
    loc_ne[i,1:6]=c(i-dim(y)[1],i-dim(y)[1]+1,i,i+1,i+dim(y)[1],i+dim(y)[1]+1)
  }
  
  lower_level_index=(1:(dim(y)[2]-1))*dim(y)[1]
  
  for(i in lower_level_index ){
    loc_ne[i,1:6]=c(i-dim(y)[1]-1,i-dim(y)[1],i-1,i,i+dim(y)[1]-1,i+dim(y)[1])
  }
  
  for(i in ((dim(y)[2]-1)*dim(y)[1]+2):((dim(y)[2])*dim(y)[1]-1) ){
    loc_ne[i,1:6]=c(i-dim(y)[1]-1,i-dim(y)[1],i-dim(y)[1]+1,i-1,i,i+1)
  }
  ## neighbors of 4 corner points
  loc_ne[1,1:4]=c(1,2,dim(y)[1]+1,dim(y)[1]+2)
  loc_ne[dim(y)[1],1:4]=c(dim(y)[1]-1,dim(y)[1],2*dim(y)[1]-1,2*dim(y)[1])
  loc_ne[dim(y)[1]*(dim(y)[2]-1)+1,1:4]=c(dim(y)[1]*(dim(y)[2]-2)+1,dim(y)[1]*(dim(y)[2]-1)+2,dim(y)[1]*(dim(y)[2]-1)+1,dim(y)[1]*(dim(y)[2]-1)+2 )
  loc_ne[dim(y)[1]*(dim(y)[2]),1:4]=c(dim(y)[1]*(dim(y)[2]-1)-1,dim(y)[1]*(dim(y)[2]-1),dim(y)[1]*(dim(y)[2])-1,dim(y)[1]*(dim(y)[2]))
  
  return(loc_ne)
  
}

## function that fills NA values by smoothing
fill_NA<-function(output){
  output_all_mat=matrix(output,dim(output)[1]*dim(output)[2],dim(output)[3])
  
  # simply smooth
  output_all_mat_missing_filled=output_all_mat
  #around 1% missing
 
  for(i in 1:dim(output_all_mat)[1]){
    for(j in 1:dim(output_all_mat)[2]){
      if(is.na(output_all_mat_missing_filled[i,j])){
        output_all_mat_missing_filled[i,j]= mean(output_all_mat_missing_filled[loc_ne[i],j],na.rm=T)
      }
    }
  }
  return(output_all_mat_missing_filled)
}

## computes correlations between months and computes july-aug vpd average, jan-aug precip average 
eda_corr<-function(output_ndvi,input_vpd,input_precipitation,num_yrs_all=18){
  
  Jan_index=1+12*(0:(num_yrs_all-1))
  Feb_index=2+12*(0:(num_yrs_all-1))
  March_index=3+12*(0:(num_yrs_all-1))
  April_index=4+12*(0:(num_yrs_all-1))
  May_index=5+12*(0:(num_yrs_all-1))
  June_index=6+12*(0:(num_yrs_all-1))
  July_index=7+12*(0:(num_yrs_all-1))
  Aug_index=8+12*(0:(num_yrs_all-1))
  Aug_index=8+12*(0:(num_yrs_all-1))
  Sep_index=9+12*(0:(num_yrs_all-1))
  Oct_index=10+12*(0:(num_yrs_all-1))
  Nov_index=11+12*(0:(num_yrs_all-1))
  Dec_index=12+12*(0:(num_yrs_all-1))
  
  
  num_loc=dim(output_ndvi)[1]
  num_obs_all=dim(output_ndvi)[2]
  cor_VPD_prior_month_NDVI=rep(NA,num_loc)
  cor_VPD_July_NDVI_Aug=rep(NA,num_loc)
  cor_VPD_JulyAug_NDVI_Aug=rep(NA,num_loc)
  cor_VPD_Aug_NDVI_Aug=rep(NA,num_loc)
  cor_VPD_Jan_to_Aug_NDVI_Aug=rep(NA,num_loc)
  
  cor_preci_prior_month_NDVI=rep(NA,num_loc)
  cor_preci_July_NDVI_Aug=rep(NA,num_loc)
  cor_preci_JulyAug_NDVI_Aug=rep(NA,num_loc)
  cor_preci_Aug_NDVI_Aug=rep(NA,num_loc)
  cor_preci_Jan_to_Aug_NDVI_Aug=rep(NA,num_loc)
  
  cor_NDVI_June_NDVI_Aug=rep(NA,num_loc)
  cor_NDVI_July_NDVI_Aug=rep(NA,num_loc)
  cor_NDVI_JuneJuly_NDVI_Aug=rep(NA,num_loc)
  
  
  VPD_loc_JulyAug_mean_record=matrix(NA,num_loc,num_yrs_all)
  VPD_loc_JanAug_mean_record=matrix(NA,num_loc,num_yrs_all)
  
  Preci_loc_Jan_to_Aug_mean_record=matrix(NA,num_loc,num_yrs_all)
  
  
  for(i_loc in 1:num_loc){
    cor_VPD_prior_month_NDVI[i_loc]=cor(input_vpd[i_loc,1:(num_obs_all-1)],output_ndvi[i_loc,2:num_obs_all])
    cor_VPD_July_NDVI_Aug[i_loc]=cor(input_vpd[i_loc,July_index],output_ndvi[i_loc,Aug_index])
    VPD_i_loc_JulyAug_mean=rowMeans(cbind(input_vpd[i_loc,July_index],input_vpd[i_loc,Aug_index]))
    VPD_loc_JulyAug_mean_record[i_loc,]=VPD_i_loc_JulyAug_mean
    
    cor_VPD_JulyAug_NDVI_Aug[i_loc]=cor(VPD_i_loc_JulyAug_mean,output_ndvi[i_loc,Aug_index])
    cor_VPD_Aug_NDVI_Aug[i_loc]=cor(input_vpd[i_loc,Aug_index],output_ndvi[i_loc,Aug_index])
    VPD_i_loc_Jan_to_Aug_mean=rowMeans(cbind(input_vpd[i_loc,Jan_index],input_vpd[i_loc,Feb_index],
                                             input_vpd[i_loc,March_index],input_vpd[i_loc,April_index],
                                             input_vpd[i_loc,May_index],input_vpd[i_loc,June_index],
                                             input_vpd[i_loc,July_index],input_vpd[i_loc,Aug_index]))
    VPD_loc_JanAug_mean_record[i_loc,]=VPD_i_loc_Jan_to_Aug_mean
    cor_VPD_Jan_to_Aug_NDVI_Aug[i_loc]=cor(VPD_i_loc_Jan_to_Aug_mean,output_ndvi[i_loc,Aug_index])
    
    
    ##precipretation
    cor_preci_prior_month_NDVI[i_loc]=cor(input_precipitation[i_loc,1:(num_obs_all-1)],output_ndvi[i_loc,2:num_obs_all])
    cor_preci_July_NDVI_Aug[i_loc]=cor(input_precipitation[i_loc,July_index],output_ndvi[i_loc,Aug_index])
    Preci_i_loc_JulyAug_mean=rowMeans(cbind(input_precipitation[i_loc,July_index],input_precipitation[i_loc,Aug_index]))
    cor_preci_JulyAug_NDVI_Aug[i_loc]=cor(Preci_i_loc_JulyAug_mean,output_ndvi[i_loc,Aug_index])
    cor_preci_Aug_NDVI_Aug[i_loc]=cor(input_precipitation[i_loc,Aug_index],output_ndvi[i_loc,Aug_index])
    Preci_i_loc_Jan_to_Aug_mean=rowMeans(cbind(input_precipitation[i_loc,Jan_index],input_precipitation[i_loc,Feb_index],
                                               input_precipitation[i_loc,March_index],input_precipitation[i_loc,April_index],
                                               input_precipitation[i_loc,May_index],input_precipitation[i_loc,June_index],
                                               input_precipitation[i_loc,July_index],input_precipitation[i_loc,Aug_index]))
    Preci_loc_Jan_to_Aug_mean_record[i_loc,]=Preci_i_loc_Jan_to_Aug_mean
    
    cor_preci_Jan_to_Aug_NDVI_Aug[i_loc]=cor(Preci_i_loc_Jan_to_Aug_mean,output_ndvi[i_loc,Aug_index])
    
    ##self correlation
    cor_NDVI_June_NDVI_Aug[i_loc]=cor(output_ndvi[i_loc,Aug_index-2],output_ndvi[i_loc,Aug_index])
    cor_NDVI_July_NDVI_Aug[i_loc]=cor(output_ndvi[i_loc,Aug_index-1],output_ndvi[i_loc,Aug_index])
    
    mean_NDVI_JuneJuly=rowMeans(cbind(output_ndvi[i_loc,Aug_index-2],output_ndvi[i_loc,Aug_index-1]))
    cor_NDVI_JuneJuly_NDVI_Aug[i_loc]=cor(mean_NDVI_JuneJuly,output_ndvi[i_loc,Aug_index])
    
  }
  
  ###VPD
  return.list=list()
  return.list$cor_VPD_prior_month_NDVI=(cor_VPD_prior_month_NDVI)
  return.list$cor_VPD_July_NDVI_Aug=(cor_VPD_July_NDVI_Aug)
  return.list$cor_VPD_JulyAug_NDVI_Aug=(cor_VPD_JulyAug_NDVI_Aug) ##again the abs is very high
  return.list$cor_VPD_Aug_NDVI_Aug=(cor_VPD_Aug_NDVI_Aug)
  return.list$cor_VPD_Jan_to_Aug_NDVI_Aug=(cor_VPD_Jan_to_Aug_NDVI_Aug) ##this one seems also high, 
  ###preci
  return.list$cor_preci_prior_month_NDVI=(cor_preci_prior_month_NDVI)
  return.list$cor_preci_July_NDVI_Aug=(cor_preci_July_NDVI_Aug)
  return.list$cor_preci_Aug_NDVI_Aug=(cor_preci_Aug_NDVI_Aug) ##
  return.list$cor_preci_JulyAug_NDVI_Aug=(cor_preci_JulyAug_NDVI_Aug) ## 0.467
  return.list$cor_preci_Jan_to_Aug_NDVI_Aug=(cor_preci_Jan_to_Aug_NDVI_Aug) ##this one is quite high, 0.6455984
  
  ##NDVI
  return.list$cor_NDVI_June_NDVI_Aug=(cor_NDVI_June_NDVI_Aug)
  return.list$cor_NDVI_JuneJuly_NDVI_Aug=(cor_NDVI_JuneJuly_NDVI_Aug)
  return.list$cor_NDVI_July_NDVI_Aug=(cor_NDVI_July_NDVI_Aug)
  ##record
  return.list$NDVI_loc_Aug_record=(output_ndvi[,Aug_index])
  return.list$VPD_loc_JulyAug_mean_record=(VPD_loc_JulyAug_mean_record)
  return.list$VPD_loc_JanAug_mean_record=(VPD_loc_JanAug_mean_record)
  return.list$Preci_loc_Jan_to_Aug_mean_record=(Preci_loc_Jan_to_Aug_mean_record)
  
  return(return.list)
}

## function which plots a column structured data back into a map image (NDVI)
plot_column_as_image_ndvi <- function(column_data, admin.sf, base_raster, min_limit = 0, max_limit =1, color = 'Greens', legend_name ="") {
  # column to 100x140 matrix
  matrix_data <- matrix(column_data, nrow = 100, ncol = 140)
  
  # matrix to raster
  raster_data <- raster(matrix_data)
  extent(raster_data) <- extent(base_raster)
  
  # raster to dataframe
  df <- data.frame(rasterToPoints(raster_data, spatial = TRUE))
  colnames(df) <- c("layer","x","y")
  
  # image
  plot <- ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = layer)) + 
    geom_sf(data = admin.sf, fill = alpha("black", 0)) + 
    xlab("") + ylab("") + theme_bw() + 
    scale_fill_gradientn(name = legend_name, colours = brewer.pal(9,color),limits=c(min_limit,max_limit)) + 
    scale_x_continuous(breaks = seq(-112, -105, by = 3.5)) +
    scale_y_continuous(breaks = seq(39, 34, by = -2.5)) +
    theme(legend.text = element_text(size = 18))+
    theme(legend.title = element_text(size = 18))+
    theme(legend.box.just = "center")+
    theme(legend.key.height = unit(2, "cm"))+
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16))
  
  return(plot)
} 


#linear model where input of test year can be or doesn't have to be predicted values of vpd and precip
lm_update_each_year_ndvi<-function(testing_index=11:18,output_all,
                                   input_1=NA,input_2=NA,
                                   pred_vpd_JulyAug_testing_est_par = NULL, 
                                   pred_precip = NULL){
  
  pred_lm=matrix(NA,num_loc,length(testing_index))
  interval_length_lm_model_record=matrix(NA,num_loc,length(testing_index))
  prop95_lm_model_record=matrix(NA,num_loc,length(testing_index))
  
  for(i_test in testing_index){
    
    for(i_loc in 1:num_loc){
      
      #ndvi output at training years at specific location
      output=as.vector(output_all[i_loc,1:(i_test-1)])
      
      
      #get training input of both vpd and precip
      input_vpd=cbind(input_1[i_loc,1:(i_test-1)])
      input_precip=cbind(input_2[i_loc,1:(i_test-1)])
      
      #predicted VPD at test location (pred_vpd_JulyAug_testing_est_par is only length of testing_index)
      testing_input_vpd = pred_vpd_JulyAug_testing_est_par[i_loc, i_test - min(testing_index) + 1]
      
      #predicted precip at test location (pred_precip is only length of testing_index)
      testing_input_precip = pred_precip[i_loc, i_test - min(testing_index) + 1]
      
      df_pixel <- data.frame(August_NDVI = output,
                             avg_vpd_Jul_Aug = input_vpd,
                             avg_precip_Jan_Aug = input_precip)
      
      # using training years fit model
      model_pixel <- lm(August_NDVI ~ avg_vpd_Jul_Aug + avg_precip_Jan_Aug, data = df_pixel)
      
      df_pixel_test <- data.frame(avg_vpd_Jul_Aug = testing_input_vpd,
                                  avg_precip_Jan_Aug = testing_input_precip)
      
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


#linear model for covariate
lm_update_each_year_covariate<-function(testing_index=11:18,output_all){
  
  pred_lm_covariate=matrix(NA,num_loc,length(testing_index))
  interval_length_lm_model_record=matrix(NA,num_loc,length(testing_index))
  prop95_lm_model_record=matrix(NA,num_loc,length(testing_index))
  l2_pred_error_lm_model_record = matrix(NA,num_loc,length(testing_index))
  
  for(i_test in testing_index){
    
    for(i_loc in 1:num_loc){
      
      #vpd_output
      output=as.vector(output_all[i_loc,1:(i_test-1)])
      
      input= 1:(i_test - 1)
      
      testing_input = i_test
      
      df_pixel <- data.frame(covariate = output,
                             index = input)
      
      # train lm using training years
      model_pixel <- lm(covariate ~ index, data = df_pixel)
      
      df_pixel_test <- data.frame(index = testing_input)
      
      # predict using the linear model on the test years
      predictions <- predict(model_pixel, newdata = df_pixel_test, interval='predict')
      
      pred_lm_covariate[i_loc,i_test-testing_index[1]+1]=predictions[,1]
      interval_length_lm_model_record[i_loc,i_test-testing_index[1]+1] = predictions[,3]-predictions[,2]
      prop95_lm_model_record[i_loc,i_test-testing_index[1]+1] = (predictions[,3]>output_all[i_loc,i_test])*(predictions[,2]<output_all[i_loc,i_test])
      l2_pred_error_lm_model_record[i_loc, i_test-testing_index[1]+1]= (predictions[,1]-output_all[i_loc,i_test])^2
      
    }
  }
  
  return.list=list()
  return.list$pred_lm_covariate=pred_lm_covariate
  return.list$interval_length_lm_model_record =interval_length_lm_model_record
  return.list$prop95_lm_model_record =prop95_lm_model_record
  return.list$l2_pred_error_lm_model=l2_pred_error_lm_model_record
  
  return(return.list)
}


#calculate log posterior given samples from the data set
##this is for matern 2.5, a function to add up the log marginal post
##using TWO covariates in the past year, VPD_samples,precip_samples without past NDVI
log_post_samples_funct_two<-function(param,VPD_samples,precip_samples,NDVI_output){ ###one can add more argument
  if(is.null(dim(VPD_samples))){
    VPD_samples = matrix(VPD_samples, nrow=1)
    precip_samples = matrix(precip_samples, nrow=1)
    NDVI_output = matrix(NDVI_output, nrow=1)
  }
  beta_samples=exp(param[-length(param)])
  nugget_samples=exp(param[length(param)])
  log_post=0
  for(i_loc in 1:dim(VPD_samples)[1]){
    input=cbind(as.vector(VPD_samples[i_loc,]),as.vector(precip_samples[i_loc,]))
    output=as.vector(NDVI_output[i_loc,])
    
    m_rgasp_samples=rgasp(design=input,response=output,nugget.est =F,
                          range.par =1/beta_samples,nugget=nugget_samples )
    
    log_marginal_lik_samples=log_marginal_lik(param=param, nugget=m_rgasp_samples@nugget, nugget_est=F, 
                                              R0=m_rgasp_samples@R0, X=m_rgasp_samples@X, zero_mean=m_rgasp_samples@zero_mean,
                                              output=m_rgasp_samples@output, 
                                              kernel_type=as.integer(c(3,3)), ####be careful here, 3 indicates matern 2.5
                                              alpha=m_rgasp_samples@alpha)
    
    log_post=log_post+log_marginal_lik_samples
  }
  ##add prior
  log_post=log_post+log_approx_ref_prior(param=param,nugget=m_rgasp_samples@nugget,nugget_est=m_rgasp_samples@nugget,
                                         CL=m_rgasp_samples@CL,a=0.2,b=1/(length(output))^{
                                           1/dim(as.matrix(input))[2]
                                         } * (0.2 + dim(as.matrix(input))[2]));
  
  return(log_post)
}


#parameter tuning/optimization function for precipitation
preci_parameter_tuning = function(param, t){
  ## t is the number of training data NOT including the validation data
  ## with the given parameter -- return RMSE
  y_D = t(delta_preci_janaug_grid[sample_loc,1:t])
  R = ar1_cor(n=t, rho=param[1])
  RR = R + param[2]*diag(t)
  L = t(chol(RR))
  one = matrix(1, nrow=t, ncol=1)
  r1 = param[1]^matrix(t:1, ncol=1)
  r2 = param[1]^matrix((t+1):2, ncol=1)
  r3 = param[1]^matrix((t+2):3, ncol=1)
  R_inv_one = (backsolve(t(L), forwardsolve(L,one)))
  R_inv_r1 = (backsolve(t(L), forwardsolve(L,r1)))
  R_inv_r2 = (backsolve(t(L), forwardsolve(L,r2)))
  w1 = (1 - t(r1)%*%R_inv_one) %*% (1/ (t(one)%*%R_inv_one)) %*% t(R_inv_one)  + t(R_inv_r1)
  w2 = (1 - t(r2)%*%R_inv_one) %*% (1/ (t(one)%*%R_inv_one)) %*% t(R_inv_one)  + t(R_inv_r2)
  W = rbind(w1, w2)
  
  y_hat = W %*% y_D
  RMSE = sqrt(mean((t(delta_preci_janaug_grid[sample_loc,c(t+1, t+2)]) - y_hat)^2)) 
  return(RMSE)
}

#function that returns predicted precipitation values given year and parameters 
preci_predict = function(param, t){
  res_list = list()
  y_D = t(delta_preci_janaug_grid[,1:t])
  R = ar1_cor(n=t, rho=param[1])
  RR = R + param[2]*diag(t)
  L = t(chol(RR))
  one = matrix(1, nrow=t, ncol=1)
  r1 = param[1]^matrix(t:1, ncol=1)
  R_inv_one = (backsolve(t(L), forwardsolve(L,one)))
  R_inv_r1 = (backsolve(t(L), forwardsolve(L,r1)))
  theta_hat_mat = (1/ (t(one)%*%R_inv_one)) %*% t( t(y_D) %*% R_inv_one) 
  c_star = 1- t(r1)%*%R_inv_r1 + t(1-t(one)%*%R_inv_r1) %*% (1/ (t(one)%*%R_inv_one)) %*% (1-t(one)%*%R_inv_r1)
  sigma_hat_2_mat = diag(t(y_D - one%*%theta_hat_mat) %*% solve(RR) %*% (y_D - one%*%theta_hat_mat)) / (t-1) 
  w1 = (1 - t(r1)%*%R_inv_one) %*% (1/ (t(one)%*%R_inv_one)) %*% t(R_inv_one)  + t(R_inv_r1)
  W = rbind(w1)
  y_hat = W %*% y_D
  pred_var = as.numeric(c_star) * sigma_hat_2_mat
  res_list$pred = y_hat
  res_list$var = pred_var
  return(res_list)
}

#vpd nrmse forecast function for optimization
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
                          response=t(list_eda_ndvi$VPD_loc_JulyAug_mean_record[sample_loc,training_index_aug]),
                          nugget.est=F,nugget=exp(est_all_par_vpd[2]),range.par=1/exp(est_all_par_vpd[1]),
                          kernel_type=kernel_type,alpha=alpha )
  
  
  # predict vpd values
  pred.model.vpd=predict(ppgasp.model.vpd,testing_input=testing_input )
  
  #metrics
  sd_baseline=sd(list_eda_ndvi$VPD_loc_JulyAug_mean_record[,testing_input])
  NRMSE_vpd=sqrt(mean((list_eda_ndvi$VPD_loc_JulyAug_mean_record[sample_loc,testing_input]-t(pred.model.vpd$mean) )^2))/sd_baseline
  return(NRMSE_vpd)
  
}

#gppgp rmse validation function for one covariate instead of two
rmse_gppgp_validation_one<-function(param,input_samples,NDVI_output, d=2){ ###one can add more argument
  if(is.null(dim(input_samples))){
    #VPD_samples = matrix(VPD_samples, nrow=1)
    input_samples = matrix(input_samples, nrow=1)
    NDVI_output = matrix(NDVI_output, nrow=1)
  }
  beta_samples=exp(param[-length(param)])
  nugget_samples=exp(param[length(param)])
  n_years = ncol(input_samples)
  sse_gross = 0
  for(i_loc in 1:dim(input_samples)[1]){
    input=cbind(as.vector(input_samples[i_loc,1:(n_years-d)]))
    output=as.vector(NDVI_output[i_loc,1:(n_years-d)])
    
    m_rgasp_samples=rgasp(design=input,response=output,nugget.est =F,
                          range.par = 1/beta_samples,nugget=nugget_samples )
    
    testing_input=cbind(as.vector(input_samples[i_loc,(n_years-d+1):n_years]))
    pred_rgasp_i_grid=predict(m_rgasp_samples,testing_input)
    sse_loc = sum((NDVI_output[i_loc,(n_years-d+1):n_years] - pred_rgasp_i_grid$mean)^2)
    
    sse_gross=sse_gross+sse_loc
  }

  return(sqrt(sse_gross/(dim(input_samples)[1]*d)))
}


#linear model attribution model with only vpd input
lm_update_each_year_ndvi_input_VPD<-function(testing_index=11:18,output_all,
                                             input_1=NA,
                                             pred_vpd_JulyAug_testing_est_par = NULL){
  
  pred_lm=matrix(NA,num_loc,length(testing_index))
  interval_length_lm_model_record=matrix(NA,num_loc,length(testing_index))
  prop95_lm_model_record=matrix(NA,num_loc,length(testing_index))
  
  for(i_test in testing_index){
    
    for(i_loc in 1:num_loc){
      
      #ndvi output at training years at specific location
      output=as.vector(output_all[i_loc,1:(i_test-1)])
      
      
      #get training input of both vpd 
      input_vpd=cbind(input_1[i_loc,1:(i_test-1)])
      
      
      #gets predicted VPD at test location (pred_vpd_JulyAug_testing_est_par is only length of testing_index)
      testing_input_vpd = pred_vpd_JulyAug_testing_est_par[i_loc, i_test - min(testing_index) + 1]
      
      df_pixel <- data.frame(August_NDVI = output,
                             avg_vpd_Jul_Aug = input_vpd)
      
      # Eq 4 in Emily's paper, using training years
      model_pixel <- lm(August_NDVI ~ avg_vpd_Jul_Aug, data = df_pixel)
      
      df_pixel_test <- data.frame(avg_vpd_Jul_Aug = testing_input_vpd)
      
      # predict using the linear model on the test years
      predictions <- predict(model_pixel, newdata = df_pixel_test, interval='predict')
      
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

#linear model attribution model with input being only precipitation
lm_update_each_year_precip<-function(testing_index=11:18,output_all,
                                     input_1=NA,
                                     pred_precip = NULL){
  
  pred_lm=matrix(NA,num_loc,length(testing_index))
  interval_length_lm_model_record=matrix(NA,num_loc,length(testing_index))
  prop95_lm_model_record=matrix(NA,num_loc,length(testing_index))
  
  for(i_test in testing_index){
    
    for(i_loc in 1:num_loc){
      
      #ndvi output at training years at specific location
      output=as.vector(output_all[i_loc,1:(i_test-1)])
      
      input_precip=cbind(input_1[i_loc,1:(i_test-1)])
      
      #gets predicted precip at test location (pred_precip is only length of testing_index)
      testing_input_precip = pred_precip[i_loc, i_test - min(testing_index) + 1]
      
      df_pixel <- data.frame(August_NDVI = output,
                             avg_precip_Jan_Aug = input_precip)
      
      # Eq 4 in Emily's paper, using training years
      model_pixel <- lm(August_NDVI ~ avg_precip_Jan_Aug, data = df_pixel)
      
      df_pixel_test <- data.frame(avg_precip_Jan_Aug = testing_input_precip)
      
      # predict using the linear model on the test years
      predictions <- predict(model_pixel, newdata = df_pixel_test, interval='predict')
      
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