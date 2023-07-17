#' @title Run the petchey extension partition on example data
#' 
#' @description Applies the Petchey extension partition to data from the
#' experiment presented in Gamfeldt and Hillebrand (2008, PLoSone)
#' 
#' note: How to deal with calculating percentage contributions from the
#' fox partition when there are negative numbers... that's something
#' for another day

# clear environment
rm(list = ls())

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# load relevant scripts
source("code/02-petchey-2003-extension.R")

# load the phytoplankton data list
dat_list <- readRDS("data/data_PLoSone_2011_formatted.rds")

# first check the fox partition
FP <- local_scale_part(data = dat_list$MYS, RYe = rep(1/5, 5), part = "fox_2005")
print(FP)

# run the petchy extension partition without stock data
PE <- petchey_2003_extension(stock_data = NA ,
                             fluxM_data = dat_list$MF, fluxY_data = dat_list$YF, 
                             RYe = rep(1/5, 5))
print(PE)

# run the petchy extension with stock data
PES <- 
  petchey_2003_extension(stock_data = dat_list$MYS ,
                         fluxM_data = dat_list$MF, fluxY_data = dat_list$YF, 
                         RYe = rep(1/5, 5)
                         )
print(PES)

