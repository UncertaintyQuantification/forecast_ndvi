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


##5. Figure 1 
#1a true ndvi 2020
admin.sf <- st_crop(admin, extent(modis))
plot_column_as_image_ndvi(output_NDVI_all_mat_missing_filled[,Aug_index][,18],
                          admin.sf,
                          modis,
                          min_limit =0 ,
                          max_limit =1, color= 'Greens', legend_name = "NDVI")

#1b. true precip 2020
admin.sf <- st_crop(admin, extent(precip))
plot_column_as_image_ndvi(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[,18],
                          admin.sf,
                          precip,
                          min_limit =0 ,
                          max_limit =153, color= 'Blues', legend_name = "Precip\n(cm)")
#1c true vpd 2020
admin.sf <- st_crop(admin, extent(vpdmax))
plot_column_as_image_ndvi(list_eda_ndvi$VPD_loc_JulyAug_mean_record[,18],
                          admin.sf,
                          vpdmax,
                          min_limit =0 ,
                          max_limit =70, color= 'OrRd',legend_name = "VPD\n(hPa)")


#Figure 1d overlay plots. 
gross_ndvi <- colSums(output_NDVI_all_mat_missing_filled[,Aug_index])
min(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record)
max(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record)
min(gross_ndvi)
max(gross_ndvi)
years <- 2003:2020

y1 <- colMeans(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record)
y3 <- colMeans(list_eda_ndvi$VPD_loc_JulyAug_mean_record)
y2 <- gross_ndvi
  
# left y axis lim. Note, can use floor(min) and ceiling(max), but to get the axis to align on either side, I use these set values
y1_min <- 15
y1_max <- 40

# right y axis limot
y3_axis_min <- 25
y3_axis_max <- 40

# transformation such that left and right axis are scaled to one another
a_y3 <- (y1_max - y1_min)/(y3_axis_max - y3_axis_min)
b_y3 <- y1_min - a_y3*y3_axis_min
y3_scaled <- a_y3*y3 + b_y3

# scale ndvi histogram
min_height <- 0.1
ndvi_range <- y1_max - y1_min
y2_scaled <- ((y2 - min(y2))/(max(y2) - min(y2)))*(1 - min_height)*ndvi_range + y1_min + min_height * ndvi_range

# margins so the plot fits in my window
par(mar = c(5, 5, 2, 7))

# plot
plot(NA, xlim = range(years) + c(-0.5, 0.5), ylim = c(y1_min, y1_max),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")

#ndvi histogram
rect(xleft = years - 0.4,
     xright = years + 0.4,
     ybottom = y1_min,
     ytop = y2_scaled,
     col = "#74C476", border = NA)

# left axis, blue line for precip
lines(x = years, y = y1, type = "l", lwd = 3, col = "blue")
points(x = years, y = y1, pch = 16, col = "blue", cex = 1.5)

# right axis, orange line for vpd
lines(x = years, y = y3_scaled, type = "l", lwd = 3, col = "#FC8D59")
points(x = years, y = y3_scaled, pch = 17, col = "#FC8D59", cex = 1.5)

# gross ndvi labels, with 4 sig figs
text(x = years,
     y = y2_scaled + 0.5,
     labels = signif(y2, 4),
     col = "#238B45", cex = 1.5, font = 2)

# left axis details
axis(2, at = seq(15, y1_max, by = 5),
     labels = seq(15, y1_max, by = 5),
     las = 2, col = "blue", cex.axis = 2, lwd = 3)
mtext("Overall Jan–Aug Precip Avg (cm)", side = 2, line = 3.2, col = "blue", cex = 2)

# right axis details
ticks <- seq(25, 40, by = 5)
tick_position <- a_y3*ticks + b_y3
axis(4, at = tick_position,
     labels = ticks,
     las = 2, col = "#FC8D59", cex.axis = 2, lwd = 3)
mtext("Overall Jul–Aug VPD Avg (hPa)", side = 4, line = 5.2, col = "#FC8D59", cex = 2)

# x axis
axis(1, at = years, labels = years, las = 2, pos = y1_min, cex.axis = 2)
mtext("Year", side = 1, line = 4, col = "black", cex = 2)
abline(h = y1_min, col = "black")


##6. Elevation plot, and text location plot

# elevation data
#long and lat of region and resolution
extent_values <- c(-111.975, -104.975, 33.975, 38.975) 
resolution <- 0.05

x_coords <-seq(extent_values[1], extent_values[2], by =resolution)
y_coords <- seq(extent_values[3], extent_values[4], by =resolution)
grid_points <- expand.grid(x= x_coords, y= y_coords)

grid_points_sf <- st_as_sf(grid_points, coords= c("x", "y"), crs= 4326)

#elevation data 
elev_data <- get_elev_point(locations = grid_points_sf, prj = "+proj=longlat +datum=WGS84", src = "aws")

#elevation data to matrix
elevation_matrix <- matrix(elev_data$elevation, nrow =length(y_coords), ncol = length(x_coords), byrow =TRUE)
elevation_raster <- raster(elevation_matrix, xmn = extent_values[1], 
                           xmx = extent_values[2], ymn = extent_values[3], ymx = extent_values[4])

#seems this works if assigned like this
crs(elevation_raster) <-"+proj=longlat +datum=WGS84"

#100x140 dimensions
basic_raster <- raster(nrows = 100, ncols = 140, crs = crs(elevation_raster))
extent(basic_raster) <- extent(elevation_raster)
resampled_raster <- resample(elevation_raster, basic_raster, method = 'bilinear')

#raster to a matrix to begin matching with how our data is structured 
resampled_matrix <- as.matrix(resampled_raster)
flipped_matrix <- apply(resampled_matrix, 2, rev)

# matrix to a 100x140 matrix (column-wise)
elevation_column_data <- as.vector((flipped_matrix))

plot_column_as_image_elev2 <- function(column_data, admin.sf, base_raster, min_limit = 0, max_limit =1, color = 'Greens', legend_name ="") {
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
    scale_fill_gradientn(name = legend_name, colours = rev(brewer.pal(9, color)),limits=c(min_limit,max_limit)) + 
    scale_x_continuous(breaks = seq(-112, -105, by = 3.5)) +
    scale_y_continuous(breaks = seq(39, 34, by = -2.5)) +
    theme(legend.text = element_text(size = 28))+
    theme(legend.title = element_text(size = 28))+
    theme(legend.box.just = "center")+
    theme(legend.key.height = unit(2, "cm"))+
    theme(axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 24))
  
  return(plot)
}

plot_column_as_image_elev2(elevation_column_data,
                          admin.sf,
                          modis,
                          min_limit =600 ,
                          max_limit =4000, color="Spectral", legend_name = "Elevation\n(m)")


#temporarily set it to 0 so we can see where it is on the map, and make sure function index_to_df_row is correct
elevation_column_data_1234 <- elevation_column_data
elevation_column_data_1234[1234] <- 0

plot_column_as_image_ndvi_with_text <- function(column_data, admin.sf, base_raster, min_limit = 0, max_limit =1, color = 'Greens', legend_name ="", selected_loc = NA, assigned_labels =c(1,2,3,4)) {
  # column to 100x140 matrix
  matrix_data <- matrix(column_data, nrow = 100, ncol = 140)
  
  # matrix to raster
  raster_data <- raster(matrix_data)
  extent(raster_data) <- extent(base_raster)
  
  index_to_df_row <- function(index, raster_obj, df, nrow_m){
    row <- ((index - 1) %% nrow_m) + 1
    col <- ((index - 1) %/% nrow_m) + 1
    cell <- cellFromRowCol(raster_obj, row, col)
    which(cellFromXY(raster_obj, df[, c("x", "y")]) == cell)
  }
  
  # raster to dataframe
  df <- data.frame(rasterToPoints(raster_data, spatial = TRUE))
  colnames(df) <- c("layer","x","y")
  
  label_row <- c()
  for (loc in selected_loc){
    label_index <- index_to_df_row(loc, raster_data, df, 100)
    label_row <- c(label_row, label_index)
  }
  
  # image
  plot <- ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = layer)) + 
    geom_sf(data = admin.sf, fill = alpha("black", 0)) + 
    geom_text(data = df[label_row, ], aes(x = x, y = y), label = assigned_labels, size = 12, color = "black") +
    xlab("") + ylab("") + theme_bw() + 
    scale_fill_gradientn(name = legend_name, colours = rev(brewer.pal(9,color)),limits=c(min_limit,max_limit)) + 
    scale_x_continuous(breaks = seq(-112, -105, by = 3.5)) +
    scale_y_continuous(breaks = seq(39, 34, by = -2.5)) +
    theme(legend.text = element_text(size = 28))+
    theme(legend.title = element_text(size = 28))+
    theme(legend.box.just = "center")+
    theme(legend.key.height = unit(2, "cm"))+
    theme(axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 24))
  
  return(plot)
}

# Figure 2a
plot_column_as_image_ndvi_with_text(elevation_column_data,
                                    admin.sf,
                                    modis,
                                    min_limit =600 ,
                                    max_limit =4000, legend_name = "Elev.\n(m)", color= "Spectral",
                                    selected_loc = c(5346),
                                    assigned_labels = c(1))

# Figure 4a
plot_column_as_image_ndvi_with_text(elevation_column_data,
                          admin.sf,
                          modis,
                          min_limit =600 ,
                          max_limit =4000, legend_name = "Elevation\n(m)", color="Spectral",
                          #selected_loc = c(1234,11105,4996,7340),
                          selected_loc = c(11105,7340, 4996),
                          assigned_labels = c(1,2,3))


plot_column_as_image_ndvi_with_text_green <- function(column_data, admin.sf, base_raster, min_limit = 0, max_limit =1, color = 'Greens', legend_name ="", selected_loc = NA, assigned_labels =c(1,2,3,4)) {
  # column to 100x140 matrix
  matrix_data <- matrix(column_data, nrow = 100, ncol = 140)
  
  # matrix to raster
  raster_data <- raster(matrix_data)
  extent(raster_data) <- extent(base_raster)
  
  index_to_df_row <- function(index, raster_obj, df, nrow_m){
    row <- ((index - 1) %% nrow_m) + 1
    col <- ((index - 1) %/% nrow_m) + 1
    cell <- cellFromRowCol(raster_obj, row, col)
    which(cellFromXY(raster_obj, df[, c("x", "y")]) == cell)
  }
  
  # raster to dataframe
  df <- data.frame(rasterToPoints(raster_data, spatial = TRUE))
  colnames(df) <- c("layer","x","y")
  
  label_row <- c()
  for (loc in selected_loc){
    label_index <- index_to_df_row(loc, raster_data, df, 100)
    label_row <- c(label_row, label_index)
  }
  
  # image
  plot <- ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = layer)) + 
    geom_sf(data = admin.sf, fill = alpha("black", 0)) + 
    geom_text(data = df[label_row, ], aes(x = x, y = y), label = assigned_labels, size = 12, color = "black") +
    xlab("") + ylab("") + theme_bw() + 
    scale_fill_gradientn(name = legend_name, colours = brewer.pal(9,color),limits=c(min_limit,max_limit)) + 
    scale_x_continuous(breaks = seq(-112, -105, by = 3.5)) +
    scale_y_continuous(breaks = seq(39, 34, by = -2.5)) +
    theme(legend.text = element_text(size = 28))+
    theme(legend.title = element_text(size = 28))+
    theme(legend.box.just = "center")+
    theme(legend.key.height = unit(2, "cm"))+
    theme(axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 24))
  
  return(plot)
}

# Fig 6a
plot_column_as_image_ndvi_with_text_green(output_NDVI_all_mat_missing_filled[,Aug_index][,18],
                                    admin.sf,
                                    modis,
                                    min_limit =0 ,
                                    max_limit =1, color= 'Greens', legend_name = "NDVI",
                                    selected_loc = c(5346),
                                    assigned_labels = c(1))


##7. Covariate Correlation Plots

#3a. precip acf plot
preci_Jan_to_Aug_avg=colMeans(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record)
acf(preci_Jan_to_Aug_avg)$acf[-1][1]

acf_1=matrix(NA,14000,1)
for (i_loc in 1:14000) {
  acf_values_test <- acf(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[i_loc,], plot = FALSE)$acf[-1][1]
  acf_1[i_loc]=acf_values_test
}

admin.sf <- st_crop(admin, extent(precip))
plot_column_as_image_ndvi(acf_1, admin.sf, precip, min_limit = -1, max_limit =1, color = 'PRGn', legend_name = "ACF 1")

#3b. 
i_loc_elevation_large=which(elevation_column_data >= 2500)
i_loc_elevation_small=which(elevation_column_data < 2500)

num_years <- 18
par(mar = c(5, 6, 4, 2) + 0.1)
plot(2003:2020, colMeans(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[i_loc_elevation_large,]),
     ylab = 'Jan-Aug Precipitation (cm)',
     xlab = 'Year', type = 'l', col = 'red',
     ylim = range(c(15,70)),
     lwd =3,cex.lab=2, cex.axis=2, lty = "dashed")

# plot high precip areas
lines(2003:2020,colMeans(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record), type = 'l', col = 'black',lwd =3)
# plot low precip areas
lines(2003:2020,colMeans(list_eda_ndvi$Preci_loc_Jan_to_Aug_mean_record[-i_loc_elevation_large,]), type = 'l', col = 'blue',lwd =3, lty ="dotdash")
legend("topleft", legend = c("High Elevation Regions", "Low Elevation Regions", "Whole Region"), col = c("red", "blue", "black"), lty = c(2, 4, 1), cex = 1.7,lwd =3)
dev.off()

#3c. acf plot
VPD_july_aug_avg=colMeans(list_eda_ndvi$VPD_loc_JulyAug_mean_record)
acf(VPD_july_aug_avg)$acf[-1][1]

acf_1_vpd=matrix(NA,14000,1)
for (i_loc in 1:14000) {
  acf_values_test_vpd <- acf(list_eda_ndvi$VPD_loc_JulyAug_mean_record[i_loc,], plot = FALSE)$acf[-1][1]
  acf_1_vpd[i_loc]=acf_values_test_vpd
}

admin.sf <- st_crop(admin, extent(vpdmax))
plot_column_as_image_ndvi(acf_1_vpd, admin.sf, vpdmax, min_limit = -1, max_limit =1, color = 'PRGn', legend_name = "ACF 1")

#d. 
par(mar = c(5, 6, 4, 2) + 0.1)
plot(2003:2020, colMeans(list_eda_ndvi$VPD_loc_JulyAug_mean_record[i_loc_elevation_large,]),
     ylab = 'July-Aug VPD (hPa)',
     xlab = 'Year', type = 'l', col = 'red',
     ylim = range(c(15,46)),
     lwd =3,cex.lab=2, cex.axis=2, lty = "dashed")

# plot high vpd areas
lines(2003:2020,colMeans(list_eda_ndvi$VPD_loc_JulyAug_mean_record), type = 'l', col = 'black',lwd =3)
# plot low vpd areas
lines(2003:2020,colMeans(list_eda_ndvi$VPD_loc_JulyAug_mean_record[-i_loc_elevation_large,]), type = 'l', col = 'blue',lwd =3, lty ="dotdash")
legend("top", legend = c("High Elevation Regions", "Low Elevation Regions", "Whole Region"), col = c("red", "blue", "black"), lty = c(2, 4, 1), cex = 1.7,lwd =3)
