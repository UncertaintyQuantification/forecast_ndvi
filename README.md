Code for "Long-Term Probabilistic Forecast of Vegetation by Climate Attributes in the Four Corners Region" (Erika McPhillips, Hyeongseong Lee, Xiangyu Xie, Kathy Baylis, Chris Funk, and Mengyang Gu, 2025). 

The repository contains code and data to reproduce all numerical results and figures from the paper:
Data - contains VPD, NDVI, Precipitation data. Also contains output of various models that are called in other files. 
Data_figures - contains code to create the descriptive data figures and maps in the manuscript. 
Attribution_model - Contains the code of all attribution models (G-PPGP, FNO, DNN, Linear).
One_Year_Forecast - contains code of all one-year-ahead forecasts (G-PPGP, FNO, RNN, AR(1), Linear, Location Mean, Previous Year). 
functions - Required R functions. 
Variable_selection - contains code that provides results for the linear and non linear variable feature selection. 
One_month_forecast - contains code that produces NDVI forecasts closer than one year out.  

