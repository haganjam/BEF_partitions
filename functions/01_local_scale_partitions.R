#'
#' @title local_scale_part()
#' 
#' @description Function to calculate Loreau and Hector's (2001, Nature) and
#' Fox (2005, Ecology Letters) biodiversity effect partition using species 
#' mixture and monoculture data
#' 
#' @param data data.frame in the following format
#' column 1 - sample: variable specifying the unique place-time combination
#' column 2 - species: variable specifying the species name (all times and places must have all species names present)
#' column 3 - M: monoculture functioning
#' column 4 - Y: mixture function
#' @param RYe expected relative yields for the species as a numeric vector of length = (N species) and which sums to one
#' @param part either "loreau_2001" or "fox_2005"
#' 
#' @symbol Mi - monoculture of each species
#' @symbol Yoi - observed yield of each species in mixture
#' @symbol Yo - observed mixture yield - sum(Yoi)
#' @symbol RYei - expected relative yield of each species (usually 1/n species but can be anything)
#' @symbol RYoi - observed relative yield (Yoi/Mi) i.e. measures the extent to which species i overyields
#' @symbol dRY - RYoi - RYei
#' @symbol N - n spp in mixture
#' 

# load the helper functions
source("functions/helper_functions.R")

# install and load libraries required for these functions
install_if("dplyr")
install_if("assertthat")

local_scale_part <- function(data, RYe, part = "loreau_2001") {
  
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
    
    all(names(x) %in% c("sample", "species", "M", "Y"))
    
  }
  
  assertthat::on_failure(test_2) <- function(call, env){
    
    paste0(deparse(call$x), " is missing one of the following columns: sample, species, M, Y")
    
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
  
  # calculate the net biodiversity effect
  nbe_df <- aggregate(df[, c("Y", "Ye") ], list(df$sample), sum)
  NBE <- (nbe_df$Y - nbe_df$Ye)
  
  if (part == "loreau_2001") {
    
    # calculate the complementarity effect (trait independent complementarity sensu Fox 2005)
    comp_df <- aggregate(df[, c("dRY", "M")], list(df$sample), mean)
    CE <- n_sp*(comp_df$dRY)*(comp_df$M)
    
    # calculate the selection effect
    SE <- sapply(split(df, df$sample), function(x) { n_sp*raw_cov(x$dRY, x$M) })
    
    # wrap this into a tibble()
    BEF_df <- tibble(sample = unique(df$sample),
                     NBE = NBE,
                     SE = SE,
                     CE = CE
                     )
    
  } else if (part == "fox_2005") {
    
    # calculate trait independent complementarity sensu Fox 2005
    comp_df <- aggregate(df[, c("dRY", "M")], list(df$sample), mean)
    TI_CE <- n_sp*(comp_df$dRY)*(comp_df$M)
    
    # calculate trait dependent complementarity
    TD_CE <- sapply(split(df, df$sample), 
                    function(x) { n_sp*raw_cov(x$M, (x$RYo - (x$RYo/sum(x$RYo)) ) ) })
    
    # calculate the dominance effect
    DOM <- sapply(split(df, df$sample), 
                  function(x) { n_sp*raw_cov(x$M, ((x$RYo/sum(x$RYo)) - x$RYe) ) })
    
    # wrap this into a data.frame
    BEF_df <- tibble(sample = unique(df$sample),
                     NBE = NBE,
                     TI_CE = TI_CE ,
                     TD_CE = TD_CE,
                     DOM = DOM
                     )
    
  } else {
    
    stop("Choose appropriate option for the `part` argument")
    
  }
  
  return(BEF_df)
  
}


# define examples from Fox (2005)
f1 <- data.frame(sample = rep(c(1, 2, 3, 4, 5, 6), each = 2),
                 species = rep(c(1, 2), 6),
                 M = rep(c(500, 250), 6),
                 Y = c(300, 100, 330, 110, 360, 120, 390, 130, 420, 140, 450, 150))

# define the answers from Fox (2005
f1.ans <- data.frame(SE = c(0, 2.5, 5, 7.5, 10, 12.5),
                     CE = c(0, 37.5, 75, 112.5, 150, 187.5))

# test if the function correctly calculates the biodiversity effects sensu Loreau and Hector (2001, Nature)

# use the partition to calculate the biodiversity effects
df.test <- local_scale_part(data = f1, RYe = c(0.60, 0.40), part = "loreau_2001")
SE.test <- all(near(df.test$SE, f1.ans$SE))
CE.test <- all(near(df.test$CE, f1.ans$CE))

# check if any of the biodiversity effects were incorrectly calculated
if ( any(c(SE.test, CE.test) == FALSE) ) { 
    
    warning("Functions do not correctly calculate biodiversity effects of the test data") 
    
  } else { 
    
    message("Functions correctly calculate biodiversity effects on the test data")
    
  }

# test if the function correctly calculates the biodiversity effects sensu Fox (2005, Ecology Letters)

# use the partition to calculate the biodiversity effects
df.test <- local_scale_part(data = f1, RYe = c(0.60, 0.40), part = "fox_2005")
CE.test <- near(df.test$TI_CE, f1.ans$CE)
SE.test <- near((df.test$TD_CE+df.test$DOM), f1.ans$SE)

# check if any of the biodiversity effects were incorrectly calculated
if ( any(c(SE.test, CE.test) == FALSE) ) { 
  
  warning("Functions do not correctly calculate biodiversity effects of the test data") 
  
} else { 
  
  message("Functions correctly calculate biodiversity effects on the test data")
  
}

### END
