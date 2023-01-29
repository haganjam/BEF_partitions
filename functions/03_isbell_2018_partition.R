#'
#' @title Function to implement Isbell et al.'s (2018, Ecology Letters) biodiversity effect partition
#' 
#' @description This script contains a set of functions that are used to implement
#' Isbell et al.'s (2018, Ecology Letters) partition of biodiversity effects. The script includes
#' helper functions to calculate the number of unique elements, n_unique(), and
#' a function to calculate raw covariance, raw_cov(). The function to calculate the
#' different biodiversity effects is Isbell_2018_part(). In addition, we implement an
#' extension that uses samples from the Dirichlet distribution to estimate expected
#' relative yields ()
#'

# load the helper functions
source("functions/helper_functions.R")

# install and load libraries required for these functions
install_if("dplyr")
install_if("assertthat")

#'
#' @title Isbell_2018_part()
#' 
#' @description Function to calculate Isbell et al.'s (2018, Ecology Letters) 
#' biodiversity effect partition using species mixture and monoculture data
#' across multiple times and places.
#' 
#' @param data data.frame in the format defined by Isbell et al. (2018, Ecology Letters):
#' column 1 - sample: variable specifying the unique place-time combination
#' column 2 - place: variable specifying the place
#' column 3 - time: variable specifying the time-point
#' column 4 - species: variable specifying the species name (all times and places must have all species names present)
#' column 5 - M: monoculture functioning
#' column 6 - Y: mixture function
#' @param RYe expected relative yields for the species across in all times and places:
#' numeric vector of length = (N species) and which sums to one
#' 
#' @symbol Mi - monoculture of each species
#' @symbol Yoi - observed yield of each species in mixture
#' @symbol Yo - observed mixture yield - sum(Yoi)
#' @symbol RYei - expected relative yield of each species (usually 1/n species but can be anything)
#' @symbol RYoi - observed relative yield (Yoi/Mi) i.e. measures the extent to which species i overyields
#' @symbol dRY - RYoi - RYei
#' @symbol N - n spp in mixture
#' @symbol Poi - observed proportion of each species in mixture (i.e. RYoi/sum(RYoi))
#' 

Isbell_2018_part <- function(data, RYe) {
  
  # test if the input data is a data.frame
  test_1 <- function(x) {
    
    is.data.frame(x)
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(call$x, " is not a data.frame")
    
  }
  
  assertthat::assert_that(test_1(x = data))
  
  # test if all the required columns are present
  test_2 <- function(x) {
    
    all(names(x) %in% c("sample", "place", "time", "species", "M", "Y"))
    
  }
  
  assertthat::on_failure(test_2) <- function(call, env){
    
    paste0(deparse(call$x), " is missing one of the following columns: sample, place, time, species, M, Y")
    
  }
  
  assertthat::assert_that(test_2(x = data))
  
  # test if RYe is a number
  test_3 <- function(x) {
    
   is.vector(x) & is.numeric(x)
    
  }
  
  assertthat::on_failure(test_3) <- function(call, env){
    
    paste0(call$x, " is not a numeric vector")
    
  }
  
  assertthat::assert_that(test_3(x = RYe))
  
  # test if the length of the RYe vector equals the number of species
  test_4 <- function(x, y) {
    
    assertthat::are_equal(length(x), y)
    
  }
  
  assertthat::on_failure(test_4) <- function(call, env){
    
    paste0("x and y do not have equal length")
    
  }
  
  assertthat::assert_that(test_4(x = RYe, y = n_unique(data[["species"]])) )
  
  # test if the RYe vector sums to 1
  test_5 <- function(x) {
    
    sum(x) > 0.99
    
  }
  
  assertthat::on_failure(test_5) <- function(call, env){
    
    paste0("RYe does not sum to one")
    
  }
  
  assertthat::assert_that(test_5(x = RYe) )
  
  # get the number of species in the data
  n_sp <- n_unique(data$species)
  
  # get the number of times
  n_t <- n_unique(data$time)
  
  # get the number places
  n_p <- n_unique(data$place)
  
  # sort the data.frame
  df <- 
    data %>%
    arrange(sample, time, place, species)
  
  # define expected relative yields
  df$RYe <- rep(RYe, n_unique(df$sample))
  
  # define observed relative yields: Prevent dividing by zero
  df$RYo <- ifelse(df$M == 0, 0, (df$Y/df$M))
  
  # define the change in relative yield
  df$dRY <- (df$RYo - df$RYe)
  
  # calculate expected yield for each species
  df$Ye <- (df$M*df$RYe)
  
  # calculate observed proportion of each species in mixture (po,ijk, Isbell et al. 2018)
  df <- 
    df %>%
    group_by(sample) %>%
    mutate(Poi = if_else(is.na(Y/sum(Y)), 0, Y/sum(Y))  ) %>%
    ungroup()
  
  # calculate change in observed proportion relative to the expectation (d.po,ijk, Isbell et al. 2018)
  df$d.Poi <- (df$Poi - df$RYe)
  
  # calculate change in observed proportion (dRYo,ijk Isbell et al. 2018)
  df$d.RYoi <- (df$RYo - df$Poi)
  
  # add means from different hierarchical levels
  # species means for each time across places (pij and Mij single bar, Isbell et al. 2018)
  sm_t <- aggregate(df[, c("d.Poi", "M") ], list(df$time, df$species), mean)
  names(sm_t) <- c("time", "species", "d.Poi.t", "M.t")
  
  # species means for each place across times (pik and Mik single bar, Isbell et al. 2018)
  sm_p <- aggregate(df[, c("d.Poi", "M") ], list(df$place, df$species), mean)
  names(sm_p) <- c("place", "species", "d.Poi.p", "M.p")
  
  # overall species mean across all times and places
  sm_s <- aggregate(df[, c("d.Poi", "M") ], list(df$species), mean)
  names(sm_s) <- c("species", "d.Poi.s", "M.s")
  
  # calculate total number of species, times and places
  N <- n_sp*n_t*n_p
  
  # 1. Net biodiversity (E10)
  NBE <- sum(df$dRY*df$M)
  
  # 2. Total complementarity (E10)
  TC <- N * mean(df$dRY) * mean(df$M)
  
  # 3. Total selection effect (E2)
  TS <- N*raw_cov(df$dRY, df$M) 
  
  # 4. Non-random overyielding (E10)
  NO <- N * raw_cov(df$d.RYoi, df$M)
  
  # 5. Average selection (E9)
  AS <- n_t * n_p * sum((sm_s$d.Poi.s - mean(df$d.Poi)) * (sm_s$M.s - mean(df$M)))
  
  # 6. Temporal insurance (E9)
  TI_df <- merge(sm_t, sm_s)
  TI <- n_p * sum( (TI_df$d.Poi.t - TI_df$d.Poi.s)*(TI_df$M.t - TI_df$M.s) )
  
  # 7. Spatial insurance (E9)
  SI_df <- merge(sm_p, sm_s)
  SI <- n_t * sum( (SI_df$d.Poi.p - SI_df$d.Poi.s)*(SI_df$M.p - SI_df$M.s) )
  
  # 8. Spatio-temporal insurance (E9)
  ST_df <- 
    merge(
      merge(
        merge(df, sm_t, by = c("time", "species")), 
        sm_p, by = c("place", "species")), 
      sm_s, by = "species")
  
  # derive the terms in E9 separately      
  ST_T1 <- with(ST_df, (d.Poi - d.Poi.t - d.Poi.p + d.Poi.s + mean(d.Poi)) )
  ST_T2 <- with(ST_df, (M - M.t - M.p + M.s + mean(M)))
  ST <- sum(ST_T1*ST_T2)
  
  # 9. Total insurance effect (E9)
  IT <- sum( (df$d.Poi - mean(df$d.Poi))*(df$M - mean(df$M)) )
  
  # 10. Local complementarity and local selection
  LC_LS <- 
    df %>%
    group_by(sample) %>%
    summarise(dRY_m = mean(dRY),
              M_m = mean(M),
              cov_m = raw_cov(dRY, M),
              n = n()) %>%
    mutate(LC = n*dRY_m*M_m,
           LS = n*cov_m) %>%
    summarise(LC = sum(LC),
              LS = sum(LS))
  
  # 10. Local complementarity
  LC <- LC_LS$LC
  
  # 11. Local selection
  LS <- LC_LS$LS
  
  # check the internal consistency
  
  # NBE = TC + TS
  assertthat::are_equal(round(NBE, 1), round((TC + TS), 1) )
  
  # NBE = LC + LS
  assertthat::are_equal(round(NBE, 1), round((LC + LS), 1) )
  
  # TC - LC = LS - TS
  assertthat::are_equal(round((TC-LC), 1), round((LS - TS), 1) )
  
  # TS = NO + IT
  assertthat::are_equal(round(TS, 1), round((NO + IT), 1) )
  
  # IT = AS + TI + SI + ST
  assertthat::are_equal(round(IT, 1), round((AS + TI + SI + ST), 1) )

  # prepare output
  
  list(Beff = data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                         Value = c(NBE, TC, TS, NO, IT, AS, TI, SI, ST)) ,
       L.Beff = data.frame(L.Beff = c("LC", "LS"),
                           Value = c(LC, LS)),
       RYe = RYe)
  
}

#'
#' @title Test data presented in Isbell et al. (2018, Ecology Letters)
#' 
#' @description Digitise the examples presented in Isbell et al. (2018, Ecology Letters)
#' in order to test whether the partition functions calculate the correct effects. The
#' comments correspond to the location in the original paper
#' 

# table 1A
t1a <- data.frame(sample = c(1,1,2,2),
                  time = c(1,1,1,1), 
                  place = c(1,1,2,2), 
                  species = c(1,2,1,2), 
                  M = c(200,200,100,100),
                  Y = c(200,200,0,0))

t1a.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   0,    100, NA, NA, NA, NA, NA, NA))

# table 1B
t1b <- data.frame(sample = c(1,1,2,2),
                  time = c(1,1,1,1), 
                  place = c(1,1,2,2), 
                  species = c(1,2,1,2), 
                  M = c(50,350,0.44,1),
                  Y = c(8.15,291.7,0.88,0))

t1b.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   100,   0, NA, NA, NA, NA, NA, NA))

# case 1
cs1 <- data.frame(sample = rep(c(1:4), each = 2),
                  time=c(1,1,1,1,2,2,2,2), 
                  place=c(1,1,2,2,1,1,2,2), 
                  species=c(1,2,1,2,1,2,1,2),
                  M=c(100,50,100,50,100,50,100,50), 
                  Y=c(100,0,100,0,100,0,100,0))

cs1.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   0,  100,   0,   100,   100,   0,   0,   0))

# case 2
cs2 <- data.frame(sample = rep(c(1:4), each = 2), 
                  time=c(1,1,1,1,2,2,2,2), 
                  place=c(1,1,2,2,1,1,2,2), 
                  species=c(1,2,1,2,1,2,1,2),
                  M=c(100,50,100,50,50,100,50,100), 
                  Y=c(100,0,100,0,0,100,0,100))

cs2.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   0,    100,   0,   100,   0,   100,   0,   0))

# case 3
cs3 <- 
  data.frame(sample = rep(c(1:4), each = 2),
             time=c(1,1,1,1,2,2,2,2), 
             place=c(1,1,2,2,1,1,2,2), 
             species=c(1,2,1,2,1,2,1,2),
             M=c(100,50,50,100,100,50,50,100), 
             Y=c(100,0,0,100,100,0,0,100))

cs3.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   0,   100,   0,   100,   0,    0,   100,   0))

# case 4
cs4 <- 
  data.frame(sample = rep(c(1:4), each = 2),
             time=c(1,1,1,1,2,2,2,2), 
             place=c(1,1,2,2,1,1,2,2), 
             species=c(1,2,1,2,1,2,1,2),
             M=c(100,50,50,100,50,100,100,50), 
             Y=c(100,0,0,100,0,100,100,0))

cs4.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   0,   100,   0,   100,   0,    0,     0,   100))

# case 5
cs5 <- data.frame(sample = rep(c(1:4), each = 2),
                  time=c(1,1,1,1,2,2,2,2), 
                  place=c(1,1,2,2,1,1,2,2), 
                  species=c(1,2,1,2,1,2,1,2),
                  M=c(75,75,75,75,75,75,75,75), 
                  Y=c(50,50,50,50,50,50,50,50))

cs5.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   100,   0,   0,    0,   0,    0,     0,    0))

# case 6
cs6 <- data.frame(sample = rep(c(1:4), each = 2),
                  time=c(1,1,1,1,2,2,2,2), 
                  place=c(1,1,2,2,1,1,2,2), 
                  species=c(1,2,1,2,1,2,1,2),
                  M=c(100,50,100,50,100,50,100,50), 
                  Y=c(50,50,50,50,50,50,50,50))

cs6.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   150,   -50,   -50,    0,   0,    0,     0,    0))

# write the raw data into a list
test.data <- list(t1a, t1b, cs1, cs2, cs3, cs4, cs5, cs6)

# write the answer data into a list
ans.data <- list(t1a.ans, t1b.ans, cs1.ans, cs2.ans, cs3.ans, cs4.ans, cs5.ans, cs6.ans)

# run the test
results <- vector(length = length(test.data))
for (i in 1:length(test.data)) {
  
  u <- Isbell_2018_part(data = test.data[[i]], RYe = c(0.5, 0.5))
  v <- u$Beff
  
  w <- ans.data[[i]]
  x <- which(!is.na(w$Value))
  
  y <- dplyr::near(w$Value[x], v$Value[x], tol = 0.1)
  
  results[i] <- any(y != TRUE)
  
}

# check if any of the biodiversity effects were incorrectly calculated
if ( any(results) ) { 
  
  warning("Functions do not correctly calculate biodiversity effects of the test data") 
  
} else { 
  
  message("Functions correctly calculate biodiversity effects on the test data")
  
}

# remove the test data objects
rm(test.data, ans.data)

### END
