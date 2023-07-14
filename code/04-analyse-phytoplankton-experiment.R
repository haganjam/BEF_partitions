
# analyse an experiment using the Petchy extension

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# load the .csv file
dat <- read.csv(file = "data/data_PLoSone_2011.csv")

# get the lowest N_P treatment
dat_t1 <- 
  dat |>
  dplyr::filter(N_P == 2, tot_P == 0.13, tot_N == 0.26)
print(dat_t1)

# calculate P-uptake
dat_t1 <- dplyr::mutate(dat_t1, P_uptake = tot_P - P)

# calculate the mean for each species replicate
dat_t1 <- 
  dat_t1 |>
  dplyr::group_by(tot_P, tot_N, N_P, species) |>
  dplyr::summarise(biovol = mean(biovol, na.rm = TRUE),
                   P_uptake = mean(P_uptake, na.rm = TRUE), .groups = "drop")

# wrangle the data into the correct format
Y <- dplyr::filter(dat_t1, species == "All")
M <- dplyr::filter(dat_t1, species != "All")

# get the stock data
MS <- 
  M |>
  dplyr::mutate(sample = 1) |>
  dplyr::rename(M = biovol) |>
  dplyr::select(sample, species, M)





