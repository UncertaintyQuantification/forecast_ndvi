source('Functions/manuscript_functions.R')

##load some packages

packages = c("rlist","dplyr","tidyverse",
             "raster", "sf", "stringr", 
             "ggplot2", "ggpubr","RColorBrewer", "colorspace","MASS","plot3D","RobustGaSP","elevatr",
             "lattice", "viridisLite", "latticeExtra", "grid")
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
getwd()
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
  
  return(log_post)
}



dim(vpd.mean.list1$VPD_loc_JanAug_mean_record)
dim(precip.mean.list1$Precip_loc_JanAug_mean_record)
log_post_samples_funct_two(param = c(-5, -5, -7),
                           VPD_samples = vpd.mean.list1$VPD_loc_JanAug_mean_record[1:10,1:10],
                           precip_samples = precip.mean.list1$Precip_loc_JanAug_mean_record[1:10,1:10],
                           NDVI_output = output_NDVI_all_mat_missing_filled[1:10, Aug_index[1:10]])





combinations <- expand.grid(vpd = names(vpd.mean.list1), precip = names(precip.mean.list1))
set.seed(0)
###subsampling some grids to predict
num_samples=500 ##number of subsamples
sample_grid_indices=sample(num_loc,size=num_samples) 


results = c()

#counter start
count <- 0

for(comb_idx in 1:nrow(combinations)){
  vpd_comb_name = combinations[comb_idx, 1]
  preci_comb_name = combinations[comb_idx, 2] 
  VPD_samples_train = vpd.mean.list1[[vpd_comb_name]][sample_grid_indices,1:10]
  precip_samples_train = precip.mean.list1[[preci_comb_name]][sample_grid_indices,1:10]
  NDVI_output_train = output_NDVI_all_mat_missing_filled[sample_grid_indices, Aug_index[1:10]]
  
  param_ini=c(-7, -7, -9)
  
  m_est=optim(param_ini,log_post_samples_funct_two,VPD_samples=VPD_samples_train,
              precip_samples=precip_samples_train,
              NDVI_output=NDVI_output_train,
              control = list(fnscale = -1)
  ) ###if num of grids are large, do not use BFGS as the gradient is too large
  
  
  results = append(results, m_est$value)
  print(comb_idx)
}

# set seed 0 and 500 samples
cbind(combinations[order(-results),],results[order(-results)])




