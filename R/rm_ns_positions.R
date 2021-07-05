# ---------------------------------------------------------------------------- #
# SigniSite, an RPackage for easy geno to phenotype analysis and visualisation #
#     Copyright (C) 2019 Leon Eyrich Jessen                                    #
# ---------------------------------------------------------------------------- #
#' Remove positions where no significant association was identified
#'
#' @param z A matrix of z-scores
#' @inheritParams signisite_zmat
#'
#' @return A matrix of z-scores, with non-significant positions removed
#'
#' @examples
#' z_mat = get_signisite_zscores(ALIGNMENT$sequence, get_values(ALIGNMENT$fasta_header))
#' rm_ns_positions(z_mat, alpha = 0.05)
#'
#' @export
rm_ns_positions = function(z, alpha = 0.05){
  z_lim = qnorm( 1 - alpha / 2 )
  keep = apply(z, 1,
               function(z_i){
                 z_i = z_i[!is.na(z_i)]
                 ifelse( any( abs( z_i ) >= z_lim ), TRUE, FALSE)
               })
  z = z[keep,]
  return(z)
}
