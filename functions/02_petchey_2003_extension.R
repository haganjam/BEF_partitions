
# Petchey (2003) function for simple mono versus mixture expectations

# load the helper functions
source("functions/helper_functions.R")
source("functions/01_local_scale_partitions.R")

# install and load libraries required for these functions
install_if("dplyr")
install_if("assertthat")

# define some random examples

# long-term change in biomass (Design 1a)
data <- data.frame(sample = rep(c(1, 2), each = 2),
                 species = rep(c(1, 2), 2),
                 M = rep(c(500, 250), 2),
                 Y = c(300, 100, 330, 110))

RYe <- c(0.5, 0.5)

# in combination with Design 1a or with Design 2

# calculate per biomass flux measurements then generate the expectation
# from the mixture biomass data
local_scale_part(data = f1.biomass, RYe = c(0.5, 0.5), part = "fox_2005")

# define some other function data

# monoculture data
f1.flux.M <- data.frame(sample = rep(c(1, 2), each = 2),
                        species = rep(c(1, 2), 2),
                        MF = c(80, 50, 100, 50))

# mixture
f1.flux.Y <- data.frame(sample = c(1, 2),
                        YF = c(80, 90))





# get the number of species in the data
n_sp <- n_unique(data$species)

# sort the data.frame
df <- 
  data %>%
  arrange(sample, species)

# define expected relative yields
df$RYe <- rep(RYe, n_unique(df$sample))

# define observed relative yields: Prevent dividing by zero
df$RYo <- ifelse(df$M == 0, 0, (df$Y/df$M))

# define the change in relative yield
df$dRY <- (df$RYo - df$RYe)

# calculate expected yield for each species
df$Ye <- (df$M*df$RYe)

# check the observed and expected yields
df %>%
  group_by(sample) %>%
  summarise(Ye = sum(Ye),
            Y = sum(Y))

# so, biomass increases by 25 and 65 units respectively

# how much increase in flux should this be as measured in the standard way

Exp.flux <- 
  f1.flux.M %>% 
  mutate(RYe = rep(RYe, 2)) %>%
  mutate(MF_Exp = (MF*RYe)) %>%
  group_by(sample) %>%
  summarise(MF_Exp = sum(MF_Exp))
print(Exp.flux)

# calculate the total NBE of the fluxes
NBE_flux <- 
  full_join(Exp.flux, f1.flux.Y, by = "sample") %>%
  mutate(NBE_flux_total = YF - MF_Exp)

# so that's the total difference in flux rates between mono and mixture
NBE_flux  

# we parse out any effects due to changing abundance
f1.flux.Y.Exp <- full_join(f1.biomass, f1.flux.M, by = c("sample", "species") )
print(f1.flux.Y.Exp)

# calculate per biomass flux
f1.flux.Y.Exp <- 
  f1.flux.Y.Exp %>%
  mutate(MF_per_bio = MF/M)
print(f1.flux.Y.Exp )

# calculate the expected flux in the mixture
f1.flux.Y.Exp <- 
  f1.flux.Y.Exp %>%
  mutate(Y_Exp = Y*MF_per_bio) %>%
  group_by(sample) %>%
  summarise(Y_Exp = sum(Y_Exp))
print(f1.flux.Y.Exp )

# join the expected flux to the true flux data
f1.flux <- full_join(f1.flux.Y.Exp, f1.flux.Y, by = "sample")
print(f1.flux)

f1.flux <- 
  f1.flux %>%
  mutate(NBE_flux_no_abun = YF-Y_Exp)

# join these together
full_join(NBE_flux, f1.flux, by = "sample") %>%
  mutate(NBE_flux_abun = (NBE_flux_total - NBE_flux_no_abun) ) %>%
  select(sample, NBE_flux_total, NBE_flux_abun, NBE_flux_no_abun)


local_scale_part(data = f1.biomass, RYe = c(0.5, 0.5), part = "fox_2005")





