# ---------------------------------------------------------------------------- #
# SigniSite, an RPackage for easy geno to phenotype analysis and visualisation #
#     Copyright (C) 2019 Leon Eyrich Jessen                                    #
# ---------------------------------------------------------------------------- #
#' Correct a matrix of z-scores for multiple testing
#'
#' @param z A matrix of z-scores
#'
#' @return A matrix of corrected z-scores
#'
#' @examples
#' z_mat = get_signisite_zscores(ALIGNMENT$sequence, get_values(ALIGNMENT$fasta_header))
#' correct_z_matrix(z = z_mat, method = 'bonferroni')
#'
#' @export
correct_z_matrix = function(z, method = 'bonferroni'){

  # Inherited from p.adjust are the following methods for correcting:
  adj_methods = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                  "fdr", "none")
  if( !method %in% adj_methods ){
    stop(paste0("SigniSite: method must be on of: ",
               paste0(adj_methods, collapse = ', ')),
         ". See ?p.adjust for more")
  }

  cons_pos      = ( z == 0 )
  z[ cons_pos ] = NA
  n_tests       = sum(!is.na(z))
  s_mat         = sign(z)
  p_mat         = 2 * pnorm( -1 * abs( z ) )
  p_mat_adj     = p.adjust(p = p_mat, method = method, n = n_tests)
  z_mat_adj     = -1 * qnorm( p_mat_adj / 2 ) * s_mat
  z_mat_adj[ cons_pos ] = 0
  return(z_mat_adj)
}
