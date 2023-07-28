#' @title Clean data from Gamfeldt and Hillebrand (2011, PLoSone) as an example
#' dataset for the flux-stock biodiversity effect partition
#' 
#' @description This scripts cleans the raw data of monoculture and mixture
#' biovolume and P-uptake into the format required to apply it to the 
#' petchey extension function.
#' 

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# load the .csv file
dat1 <- readr::read_csv(file = "data/data_PLoSone_2011.csv")
head(dat1)

# remove empty columns
dat1 <- dplyr::select(dat1, -dplyr::starts_with("..."))

# remove the empty rows
dat1 <- dat1[rowSums(is.na(dat1)) < ncol(dat1),]

# get the correct treatments
dat1 <- 
  dat1 |>
  dplyr::filter(tot_P == 5.02)

# give each treatment a unique id
dat1[["sample"]] <- as.integer(as.factor(with(dat1, paste(tot_P, N_P, sep = "_"))))

# reorder the columns and sort the data
dat1 <- 
  dat1 |>
  dplyr::select(sample, id, tot_P, tot_N, N_P, composition, biovol, P) |>
  dplyr::arrange(sample, composition)

# get the mixture stock data
dat2 <- readr::read_csv(file = "data/data_PLoSone_mixtures_2011.csv")
head(dat2)

# remove the empty column
dat2 <- dplyr::select(dat2, -dplyr::starts_with("..."))

# get the correct treatments
dat2 <- 
  dat2 |>
  dplyr::filter(tot_P == 5.02)


# give each treatment a unique id
dat2[["sample"]] <- as.integer(as.factor(with(dat2, paste(tot_P, N_P, sep = "_"))))

# reorder the columns and sort the data
dat2 <- 
  dat2 |>
  dplyr::select(sample, tot_P, N_P, species, biovol) |>
  dplyr::arrange(sample, species)

# loop over the different samples
all( unique(dat1$sample) == unique(dat2$sample) )

# list of samples
samples <- unique(dat1$sample)

# generate output lists
MYS_list <- vector("list", length = length(samples))
MF_list <- vector("list", length = length(samples))
YF_list <- vector("list", length = length(samples))

for(i in samples) {
  
  # get the lowest N_P treatment
  dat1_t1 <- 
    dat1 |>
    dplyr::filter(sample == i)
  
  # calculate P-uptake
  dat1_t1 <- dplyr::mutate(dat1_t1, P_uptake = ((tot_P - P)/tot_P)*100 )
  
  # calculate the mean for each species replicate
  dat1_t1 <- 
    dat1_t1 |>
    dplyr::group_by(sample, tot_P, N_P, composition) |>
    dplyr::summarise(biovol = mean(biovol, na.rm = TRUE),
                     P_uptake = mean(P_uptake, na.rm = TRUE), .groups = "drop")
  
  # wrangle the data into the correct format
  Y1 <- dplyr::filter(dat1_t1, composition == "All")
  M1 <- dplyr::filter(dat1_t1, composition != "All")
  
  # get the monoculture stock data
  MS <- 
    M1 |>
    dplyr::rename(species = composition) |>
    dplyr::rename(M = biovol) |>
    dplyr::select(sample, species, M)
  
  # subset the lowest N_P treatment
  dat2_t1 <- dplyr::filter(dat2, sample == i)
  
  # calculate the average of each species
  dat2_t1 <- 
    dat2_t1 |>
    dplyr::group_by(sample, N_P, tot_P, species) |>
    dplyr::summarise(biovol = mean(biovol, na.rm = TRUE), .groups = "drop")
  
  # get the mixture stock data
  YS <- 
    dat2_t1 |>
    dplyr::rename(Y = biovol) |>
    dplyr::select(sample, species, Y)
  
  # join the mixture stock data to the monoculture stock data
  MYS <- dplyr::full_join(MS, YS, by = c("sample", "species"))
  
  # get the flux data for the monocultures and the mixture
  
  # monocultures
  MF <-
    M1 |>
    dplyr::rename(species = composition, Mflux = P_uptake) |>
    dplyr::select(sample, species, Mflux)
  
  # mixture
  YF <- 
    Y1 |>
    dplyr::rename(species = composition, Yflux = P_uptake) |>
    dplyr::select(sample, species, Yflux)
  
  # write output into the relevant lists
  MYS_list[[i]] <- MYS 
  MF_list[[i]] <- MF 
  YF_list[[i]] <- YF 
  
}

# remove those treatments without complete monoculture data
rem <- lapply(MF_list, function(x) nrow(x) ) != 0

MYS_list <- MYS_list[rem]
MF_list <- MF_list[rem] 
YF_list <- YF_list[rem]

# bind these into data.frames
MYS_df <- dplyr::bind_rows(MYS_list)
MF_df <- dplyr::bind_rows(MF_list)
YF_df <- dplyr::bind_rows(YF_list)

# save this as a .rds file
id_vars <- 
  dat1 |>
  dplyr::select(sample, tot_P, tot_N, N_P) |>
  dplyr::distinct()

# remove the samples without the monoculture data
id_vars <- 
  id_vars |>
  dplyr::filter(sample %in% (unique(MF_df$sample)) )

# bind everything into a list
output <- list(id = id_vars,
               MYS = MYS_df,
               MF = MF_df,
               YF = YF_df)

# save as a .rds file
saveRDS(object = output, file = "data/data_PLoSone_2011_formatted.rds")

### END
