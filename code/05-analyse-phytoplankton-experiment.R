#' @title Run the petchey extension partition on example data
#' 
#' @description Applies the Petchey extension partition to data from the
#' experiment presented in Gamfeldt and Hillebrand (2008, PLoSone)
#' 

# clear environment
rm(list = ls())

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(cowplot)

# load relevant scripts
source("code/02-petchey-2003-extension.R")
source("code/helper-plotting-theme.R")

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

# add the NP_ratio variable
PES$NP_ratio <- dat_list$id$N_P

# try to visualise these results
PES_long <- 
  PES |>
  dplyr::mutate(NBE_flux_abun = NBE_flux_abun_DOM + NBE_flux_abun_TD.CE + NBE_flux_abun_TI.CE) |>
  tidyr::pivot_longer(cols = dplyr::starts_with("NBE_"),
                      names_to = "NBE_type",
                      values_to = "NBE_mag") |>
  dplyr::mutate(sample = as.character(sample),
                NP_ratio = as.character(NP_ratio))

# without the fox partition
PES_long1 <- 
  PES_long |>
  dplyr::filter(NBE_type %in% c("NBE_flux_tot", "NBE_flux_abun", "NBE_flux_no_abun"))

# relevel the factor
PES_long1$NBE_type <- factor(PES_long1$NBE_type, 
                             levels = c("NBE_flux_tot", "NBE_flux_abun", "NBE_flux_no_abun"))

# rename the levels
levels(PES_long1$NBE_type) <- c("NBE total = ", " NBE abun +", " NBE no_abun")

# set the widths
widths1 <- rep(c(0.4, 0.2, 0.2), 3)

# set-up a nested colour palette
col_pal1 <- c("#AC0000", "#FF5D5D", "#FFC9C9")

# plot the results
p1 <- 
  ggplot(data = PES_long1,
       mapping = aes(x = NP_ratio, y = NBE_mag, fill = NBE_type, width = widths1)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_col(position = position_dodge2(width =0.5)) +
  xlab("N:P ratio") +
  ylab("Biodiversity effect (% P uptake)") +
  labs(fill = NULL) +
  scale_fill_manual(values = col_pal1) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key.size = unit(0.8,"line"))
plot(p1)

# plot the fox partition terms
PES_long2 <- 
  PES_long |>
  dplyr::filter( NBE_type %in% c("NBE_flux_abun", 
                                 "NBE_flux_abun_TI.CE", "NBE_flux_abun_TD.CE", "NBE_flux_abun_DOM") )

# relevel the factor
PES_long2$NBE_type <- factor(PES_long2$NBE_type)

# rename the levels
levels(PES_long2$NBE_type) <- c("NBE abun =", "DOM +", " TD CE +", " TI CE")

# set the widths
widths2 <- rep(c(0.2, 0.2, 0.2, 0.3), 3)
length(widths2)

# set-up a nested colour palette
col_pal2 <- c("#FF5D5D", "#203864", "#4674C6", "#B4C7E7")

# plot the results
p2 <- 
  ggplot(data = PES_long2,
       mapping = aes(x = NP_ratio, y = NBE_mag, fill = NBE_type, width = widths2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_col(position = position_dodge2(width = 0.6)) +
  xlab("N:P ratio") +
  ylab("") +
  labs(fill = NULL) +
  scale_fill_manual(values = col_pal2) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key.size = unit(0.8,"line"))
plot(p2)

# combine these plots
p12 <- 
  cowplot::plot_grid(p1, p2, align = "hv",
                     labels = c("a", "b"), 
                     label_size = 11,
                     label_fontface = "plain"
                     )

# export this figure as a .svg file
ggsave(filename = "figures-tables/fig_1.svg", plot = p12,
       units = "cm", width = 20, height = 9)

### END
