
# Petchey (2003) function for simple mono versus mixture expectations

# load the helper functions
source("code/helper-functions.R")
source("code/01-local-partitions.R")

# install and load libraries required for these functions
library(dplyr)
library(assertthat)

# set-up a function for the Petchey (2003) extension
petchey_2003_extension <- function(stock_data = NA,
                                   fluxM_data, fluxY_data, 
                                   RYe) {
  
  # sort out the fluxM_data
  fluxM <- 
    fluxM_data %>%
    arrange(sample, species)
  
  # define expected relative yields
  fluxM$RYe <- rep(RYe, n_unique(fluxM$sample))
  
  # get expected flux
  fluxM$Mflux_exp <- fluxM$Mflux * fluxM$RYe
  
  # calculate the net biodiversity effect
  Mflux_exp_tot_df <- aggregate(list(Yflux_exp_tot = fluxM[["Mflux_exp"]]), list(sample = fluxM$sample), sum)
  fluxY_tot <- full_join(Mflux_exp_tot_df, fluxY_data , by = "sample")
  
  # if there are no biomass stock data, we calculate the NBE and return
  if(!is.data.frame(stock_data)) {
    
    NBE_flux_tot <- 
      fluxY_tot %>%
      mutate(NBE_flux_tot = (Yflux-Yflux_exp_tot) ) %>%
      select(sample, NBE_flux_tot)
    
    return(NBE_flux_tot)
    
  } 
  
  # join the stock data to the monoculture flux data
  stock_fluxM <- full_join(stock_data, fluxM_data, by = c("sample", "species"))
  
  # calculate per biomass flux
  stock_fluxM$fluxM_per_stock <- (stock_fluxM$Mflux/stock_fluxM$M)
  
  # calculate the expected flux in the mixture
  stock_fluxM$Yflux_exp <- (stock_fluxM$Y*stock_fluxM$fluxM_per_stock)
  
  # calculate the expected mixture flux
  fluxY_no_abun_df <- aggregate(list(Yflux_exp_no_abun = stock_fluxM[["Yflux_exp"]] ), list(sample = stock_fluxM$sample), sum)
  fluxY_summary <- full_join(fluxY_no_abun_df, fluxY_tot, by = "sample")
  
  # calculate the different net biodiversity effects
  NBE_flux <- 
    fluxY_summary %>%
    mutate(NBE_flux_tot = (Yflux - Yflux_exp_tot),
           NBE_flux_no_abun = (Yflux - Yflux_exp_no_abun) )
  
  # calculate the NBE_flux_abun i.e. NBE flux due to changes in abundance
  NBE_abun <- (NBE_flux$NBE_flux_tot - NBE_flux$NBE_flux_no_abun)
  
  # separate out the NBE_flux_abun term into contributions from Fox (2005)'s mechanisms
  
  # perform the fox partition
  fox_part <- local_scale_part(data = stock_data, RYe = RYe, part = "fox_2005")
  
  # calculate the proportion of each of three terms in Fox (2005)'s partition
  x <- apply(fox_part[,c("TI_CE", "TD_CE", "DOM")], 1, function(x) x/sum(x), simplify = FALSE )
  y <- mapply(function(x, y) x*y, NBE_abun, x, SIMPLIFY = FALSE )
  z <- bind_rows(y, .id = "sample")
  names(z) <- c("sample", "NBE_flux_abun_TI.CE", "NBE_flux_abun_TD.CE", "NBE_flux_abun_DOM")
  z$sample <- as.numeric(z$sample)
  
  # join to the NBE_flux data.frame
  NBE_flux <- full_join(NBE_flux, z, by = "sample")
  
  # reorder the columns
  NBE_flux <- 
    NBE_flux %>%
    select(sample, NBE_flux_tot, 
           NBE_flux_abun_TI.CE, NBE_flux_abun_TD.CE, NBE_flux_abun_DOM,
           NBE_flux_no_abun)
  
  return(NBE_flux)
  
  }

# define some random examples
stock_data <- data.frame(sample = rep(c(1, 2), each = 2),
                         species = rep(c(1, 2), 2),
                         M = rep(c(500, 250), 2),
                         Y = c(300, 100, 330, 110))

# monoculture flux data
fluxM_data <- data.frame(sample = rep(c(1, 2), each = 2),
                         species = rep(c(1, 2), 2),
                         Mflux = c(80, 50, 100, 50))

# mixture flux data
fluxY_data <- data.frame(sample = c(1, 2),
                         Yflux = c(80, 90))

# test the function without stock data
petchey_2003_extension(stock_data = NA,
                       fluxM_data = fluxM_data, fluxY_data = fluxY_data, 
                       RYe = c(0.5, 0.5))

# test the function with stock data
petchey_2003_extension(stock_data = stock_data,
                       fluxM_data = fluxM_data, fluxY_data = fluxY_data, 
                       RYe = c(0.5, 0.5))

### END
