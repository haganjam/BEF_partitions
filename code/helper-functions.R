
#'
#' @title n_unique()
#' 
#' @description Function that outputs the number of unique elements in a vector
#' 
#' @param x vector
#' 

n_unique <- function(x) { 
  
  # test if x is a numeric vector
  test_1 <- function(x) {
    
    is.vector(x)
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0("x is not a vector")
    
  }
  
  assertthat::assert_that(test_1(x = x))
  
  # return the number of unique elements
  return( length(unique(x)) ) 
  
}

#' 
#' @title raw_cov()
#' 
#' @description Function to calculate raw covariance from two equal length,
#' numeric vectors
#' 
#' @param x numeric vector of length N
#' @param y numeric vector of length N
#' 

raw_cov <- function(x, y) {
  
  test_1 <- function(x, y) {
    
    assertthat::are_equal(length(x), length(y))
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(call$x, " & ", call$y, " do not have equal length")
    
  }
  
  assertthat::assert_that(test_1(x = x, y = y))
  
  # test if x and y are numeric vectors
  test_2 <- function(x, y) {
    
    (is.vector(x) & is.numeric(x)) & (is.vector(y) & is.numeric(y))
    
  }
  
  assertthat::on_failure(test_2) <- function(call, env){
    
    paste0("either x or y are not a numeric vectors")
    
  }
  
  assertthat::assert_that(test_2(x = x, y = y))
  
  # calculate deviation from the mean of x and y
  c1 <- x - mean(x)
  c2 <- y - mean(y)
  
  # calculate the raw covariance
  c12 <- sum( (c1*c2) )/length(c1)
  
  return(c12)
  
}

### END
