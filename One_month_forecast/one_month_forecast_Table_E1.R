source('One_Year_Forecast/precip_forecast_Table_F1_Figs2b_7abc.R')
source('One_Year_Forecast/vpd_forecast_Table_F2_Figures_2e_7def_F1.R')
library(RobustGaSP)

list_eda_ndvi=eda_corr(output_ndvi=output_NDVI_all_mat_missing_filled,
                       input_vpd=VPD_all_mat_missing_filled,
                       input_precipitation=precipitatian_all_mat_missing_filled,num_yrs_all=18)

names(list_eda_ndvi)
get_month_index<-function(num_yrs_training=10,num_yrs_testing=8){
  num_yrs_all=num_yrs_training+num_yrs_testing
  Jan_index=1+12*(0:(num_yrs_all-1))
  Feb_index=2+12*(0:(num_yrs_all-1))
  March_index=3+12*(0:(num_yrs_all-1))
  April_index=4+12*(0:(num_yrs_all-1))
  May_index=5+12*(0:(num_yrs_all-1))
  June_index=6+12*(0:(num_yrs_all-1))
  July_index=7+12*(0:(num_yrs_all-1))
  Aug_index=8+12*(0:(num_yrs_all-1))
  Sep_index=9+12*(0:(num_yrs_all-1))
  Oct_index=10+12*(0:(num_yrs_all-1))
  Nov_index=11+12*(0:(num_yrs_all-1))
  Dec_index=12+12*(0:(num_yrs_all-1))
  
  month_index_set=matrix(NA,num_yrs_all,12)
  month_index_set[,1]=Jan_index
  month_index_set[,2]=Feb_index
  month_index_set[,3]=March_index
  month_index_set[,4]=April_index
  month_index_set[,5]=May_index
  month_index_set[,6]=June_index
  month_index_set[,7]=July_index
  month_index_set[,8]=Aug_index
  month_index_set[,9]=Sep_index
  month_index_set[,10]=Oct_index
  month_index_set[,11]=Nov_index
  month_index_set[,12]=Dec_index
  
  month_index_set_training=matrix(NA,num_yrs_training,12)
  month_index_set_training[,1]=1+12*(0:(num_yrs_training-1))
  month_index_set_training[,2]=2+12*(0:(num_yrs_training-1))
  month_index_set_training[,3]=3+12*(0:(num_yrs_training-1))
  month_index_set_training[,4]=4+12*(0:(num_yrs_training-1))
  month_index_set_training[,5]=5+12*(0:(num_yrs_training-1))
  month_index_set_training[,6]=6+12*(0:(num_yrs_training-1))
  month_index_set_training[,7]=7+12*(0:(num_yrs_training-1))
  month_index_set_training[,8]=8+12*(0:(num_yrs_training-1))
  month_index_set_training[,9]=9+12*(0:(num_yrs_training-1))
  month_index_set_training[,10]=10+12*(0:(num_yrs_training-1))
  month_index_set_training[,11]=11+12*(0:(num_yrs_training-1))
  month_index_set_training[,12]=12+12*(0:(num_yrs_training-1))
  
  month_index_set_testing=matrix(NA,num_yrs_testing,12)
  month_index_set_testing[,1]=1+12*(0:(num_yrs_testing-1))
  month_index_set_testing[,2]=2+12*(0:(num_yrs_testing-1))
  month_index_set_testing[,3]=3+12*(0:(num_yrs_testing-1))
  month_index_set_testing[,4]=4+12*(0:(num_yrs_testing-1))
  month_index_set_testing[,5]=5+12*(0:(num_yrs_testing-1))
  month_index_set_testing[,6]=6+12*(0:(num_yrs_testing-1))
  month_index_set_testing[,7]=7+12*(0:(num_yrs_testing-1))
  month_index_set_testing[,8]=8+12*(0:(num_yrs_testing-1))
  month_index_set_testing[,9]=9+12*(0:(num_yrs_testing-1))
  month_index_set_testing[,10]=10+12*(0:(num_yrs_testing-1))
  month_index_set_testing[,11]=11+12*(0:(num_yrs_testing-1))
  month_index_set_testing[,12]=12+12*(0:(num_yrs_testing-1))
  
  return.list=list()
  return.list$month_index_set=month_index_set
  return.list$month_index_set_training=month_index_set_training
  return.list$month_index_set_testing=month_index_set_testing
  
  return(return.list)
}

list_month_index_set=get_month_index(num_yrs_training=10,num_yrs_testing=8)
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


num_loc=dim(output_NDVI_all_mat_missing_filled)[1] # Number of locations/pixels/grids
num_obs_all=dim(output_NDVI_all_mat_missing_filled)[2] # Number of months total 
num_yrs_all=num_obs_all/12 # Number of years total 


## ppgp
ppgp_fit_preci_mat = matrix(NA, nrow=num_loc, ncol=10)

#preci_predict(param = c(rho_vec[1], eta_vec[1]), t=1)$pred 
for(i in 1:10){
  
  ppgp_fit_preci_mat[,i] = preci_predict(param = c(rho_vec[1], eta_vec[1]), t=i)$pred 
}
sum(is.na(ppgp_fit_preci_mat))
dlm_fit_preci_mat = ppgp_fit_preci_mat + rowMeans(preci_janaug_grid[, 1:10]) %*% matrix(1, ncol=10)

up_to = 1
num_yrs_all=18
Aug_index=8+12*(0:(num_yrs_all-1))
lm_fit_preci_mat = matrix(NA, nrow=num_loc, ncol=10)
lm_preci_coef = matrix(NA, nrow=2, ncol=10)
for(i in 1:10){
  if(up_to == 1){
    fit_train = lm(preci_janaug_grid[,i]~precipitatian_all_mat_missing_filled[,list_month_index_set$month_index_set_training[i,1]])  
  }else{
    fit_train = lm(preci_janaug_grid[,i]~rowMeans(precipitatian_all_mat_missing_filled[,list_month_index_set$month_index_set_training[i,1:up_to]]))  
  }
  #print(summary(fit_train)$coefficients)
  lm_fit_preci = fitted(fit_train)
  lm_fit_preci_mat[,i] = lm_fit_preci
  lm_preci_coef[,i] = summary(fit_train)$coefficients[,1]
}
coef_mean = rowMeans(lm_preci_coef)
sum(is.na(lm_fit_preci_mat))

sum((preci_janaug_grid[,1:10]-lm_fit_preci_mat)*(dlm_fit_preci_mat-lm_fit_preci_mat))/sum((dlm_fit_preci_mat-lm_fit_preci_mat)^2)
optimize_weight_lagrange = function(dlm_mat, lm_mat){
  w_cand = sum((preci_janaug_grid[,1:10]-lm_mat)*(dlm_mat-lm_mat))/sum((dlm_mat-lm_mat)^2)
  if(w_cand<0) cat('negative weight') else if(w_cand>1) cat('weight>1') else return(w_cand)
}

(w = optimize_weight_lagrange(dlm_mat=dlm_fit_preci_mat, lm_mat=lm_fit_preci_mat))
# Jan - Jan 0.6165478 
# Jan - Feb 0.5837023
# Jan - Mar 0.4953567
# Jan - Apr 0.4089771
# Jan - May 0.3667995
# Jan - Jun 0.2989379
# Jan - Jul 0.1091241


test_preci_pred = matrix(NA, nrow=num_loc, ncol=8)
for(i in 1:8){
  if(up_to ==1 ){
    design_mat = cbind(rep(1,num_loc), precipitatian_all_mat_missing_filled[,120+list_month_index_set$month_index_set_testing[i,1]]) 
  }else{
    design_mat = cbind(rep(1,num_loc), rowMeans(precipitatian_all_mat_missing_filled[,120+list_month_index_set$month_index_set_testing[i,1:up_to]]))  
  }
  test_preci_pred[,i] = design_mat%*%matrix(coef_mean, ncol=1)
}
sum(is.na(test_preci_pred))

weigted_preci_pred = w*(dlm_pred_preci_mat) + (1-w)*test_preci_pred
sqrt(mean((preci_janaug_grid[,11:18]-weigted_preci_pred)^2))
# Jan - Jan 8.161729
# Jan - Feb 7.424498
# Jan - Mar 7.548514
# Jan - Apr 7.232915
# Jan - May 6.051593
# Jan - Jun 5.450503
# Jan - Jul 3.369651





set.seed(1234)
num_samples=500 ##number of subsamples
sample_grid_indices=sample(num_loc,size=num_samples) 
record_prediction_NDVI_weighted=matrix(NA,nrow=num_loc,ncol = num_yrs_all)
lower95_weighted = matrix(NA,nrow=num_loc,ncol = num_yrs_all)
upper95_weighted = matrix(NA,nrow=num_loc,ncol = num_yrs_all)


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
    
    testing_input_weighted=cbind(
      pred_vpd_JulyAug_testing_est_par[i_grid,(time_t-10)],
      weigted_preci_pred[i_grid,(time_t-10)]
    )
    
    
    pred_rgasp_i_grid_weighted=predict(m_rgasp_i_grid,testing_input_weighted)
    
    
    record_prediction_NDVI_weighted[i_grid,time_t]=pred_rgasp_i_grid_weighted$mean
  
    lower95_weighted[i_grid,time_t] = pred_rgasp_i_grid_weighted$lower95
    upper95_weighted[i_grid,time_t] = pred_rgasp_i_grid_weighted$upper95
  }
}


###by  PP-GP 
sqrt(mean((record_prediction_NDVI_weighted[,11:num_yrs_all]-output_NDVI_all_mat_missing_filled[,Aug_index][,11:18])^2))
# 0.02911757 with c(-7, -7, -9), size=500, seed=1234 this is the attribution
# Jan - Jan 0.03481296
# Jan - Feb 0.03515878
# Jan - Mar 0.03555282
# Jan - Apr 0.03527617
# Jan - May 0.03418972
# Jan - Jun 0.03401253
# Jan - Jul 0.03221293





#mean((output_NDVI_all_mat_missing_filled[,Aug_index][,11:18] > lower95_weighted[,11:18]) * (output_NDVI_all_mat_missing_filled[,Aug_index][,11:18] < upper95_weighted[,11:18]) ) 
# 95% coverage 0.9492768
#mean(upper95_weighted[,11:18]-lower95_weighted[,11:18])
# L95% 0.1098873


