# Code that implements the minimum interval quantile MVCI

#' Adjusted naive quantile based multivariate confidence interval (mvci)
#'
#' @description
#' Computes minimum width quantile CI separately for each column of matrix X i.e.
#' separately for each variable
#'
#' Adjusts the number of observations to exclude (i.e. the percentage a in function
#' f()) such that the total number of rows insice the confidence area mathces n-k
#'
#' @param X [N,M] numeric, Data matrix, vector valued M-dimensional observations on rows,
#' @param k [1,1] integer, Number of observations (rows) to remove
#' @param l [1,1] integer, Max number of outlier dimensions for a data row that is counted as within MVCI
#'
#' @return A list describing the multivariate confidence band.
#'
#' @export
mvci_min_width_quantile <- function( X, k, l ) {
  X <- as.matrix( X )
  n <- nrow(X)
  f <- function(a) {
    ci <- get_ci( X, a )
    total_contained( ci, X, l ) - (n-k)
  }
  aopt <- uniroot( f, c(0,0.5) )$root
  get_ci( X, aopt )
}


# A wrapper for backward combatibility
# TODO: remove this!
mvci_uniform_percentiles <- mvci_min_width_quantile


#' Naive quantile based multivariate confidence interval (mvci)
#'
#' @description
#' Mininum width quantile (1-a)*100% CI for a data matrix.
#' Applies smallest_ci() to each column of X.
#' Does not control FWER over multiple dimensions like mvci_min_width_quantile().
#'
#' @param X [N,M] numeric, Data, vector valued M-dimensional observations on rows,
#' @param alpha [1,1] numeric, Percentage of observations to leave outside the interval per dimension
#'
#' @return A list with upper and lower mvci boundary.
#' TODO: return value differs from that of other methods.
#'
#' @export
get_ci <- function( X, alpha ) {
  k <- floor( nrow(X)*alpha )
  foo <- apply( X, 2, function(x) smallest_ci( x, k ) )
  list( xu=foo[2,], xl=foo[1,] )
}


#' Mininum width confidence interval for a data vector
#'
#' @description
#' Returns (1 - k/length(x)) * 100% CI for a data vector using minimum width quantiles
#' TODO: same as min_width_quantile()
#'
#' @param c [M,1] numeric, A data vector
#' @param k [1,1] integer, Number of observations to leave outside the interval
#'
#' @return A list with upper and lower mvci boundary.
#' TODO: return value differs from that of other methods.
#'
#' @export
smallest_ci <- function(x,k) {
  n <- length(x)
  x <- sort(x)
  i <- which.min(x[(n-k):n]-x[1:(k+1)])
  c(x[i],x[n+i-k-1])
}


#' Number of ci violations of x
#'
#' @description
#' Returns the number of points in which a vector x lies outside ci
#'
#' @param ci list, A list with upper and lower mvci boundary.
#' @param x [M,1] numeric, A data vector
#'
#' @return integer, Number of values in x that are outside ci.
#'
#' @export
ci_violations <- function( ci, x ) {
  sum( ci$xu - x < 0 ) + sum( x - ci$xl < 0 )
}


#' Does ci contain x with respect to l?
#'
#' @description
#' Returns TRUE if ci contains x with respect to relaxation l and FALSE otherwise
#'
#' @param ci list, A list with upper and lower mvci boundary.
#' @param x [M,1] numeric, A data vector
#' @param l [1,1] integer, Max number of outlier dimensions for a data row that is counted as within MVCI
#'
#' @return logical, True if x violates ci <= l times. False if x violates ci more than l times.
#'
#' @export
ci_contains <- function( ci, x, l ) {
  ci_violations( ci, x ) <= l
}


#' How many rows of X are inside ci with respect to l?
#'
#' @description
#' Returns the number of rows in matrix X that are within ci with respect
#' to relaxation l.
#'
#' @param ci list, A list with upper and lower mvci boundary.
#' @param X [M,N] numeric, Data matrix with vector valued M-dimensional observations on rows
#' @param l [1,1] integer, Max number of outlier dimensions for a data row that is counted as within MVCI
#'
#' @return integer, Number of rows of X that are inside ci
#'
#' @export
total_contained <- function( ci, X, l ) {
  sum( apply( X, 1, function(x) ci_contains( ci, x, l ) ) )
}
