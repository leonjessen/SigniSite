# ---------------------------------------------------------------------------- #
# SigniSite, an RPackage for easy geno to phenotype analysis and visualisation #
#     Copyright (C) 2019 Leon Eyrich Jessen                                    #
# ---------------------------------------------------------------------------- #
#' Compute the tie correction factor
#'
#' @param x A numeric vector
#'
#' @return The computed tie correction factor as a numeric vector of length 1
#'
#' @examples
#' tie_c_factor(x = round(rnorm(10,0)))
#'
#' @export
tie_c_factor = function(x){

  # Check
  if( !is.numeric(x) ){
    stop("SigniSite: 'x' in tie correction factor is not numeric")
  }
  if( !is.vector(x) ){
    stop("SigniSite: 'x' in tie correction factor is not a vector")
  }

  # Compute tie correction factor
  N = length(x)
  T_vec = unname(table(x))
  tcf = 1 - (1 / ( N ** 3 - N )) * sum( T_vec ** 3 - T_vec )

  # Done
  return( tcf )
}
