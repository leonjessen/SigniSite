# ---------------------------------------------------------------------------- #
# SigniSite, an RPackage for easy geno to phenotype analysis and visualisation #
#     Copyright (C) 2019 Leon Eyrich Jessen                                    #
# ---------------------------------------------------------------------------- #
#' Run the SigniSite workflow
#'
#' @param file A multiple sequence alignment in FASTA format
#' @param method method for correcting for multiple testing, one of c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), see ?p.adjust for more. Default: bonferroni
#' @param alpha the level of significance. Default: 0.05
#'
#' @return A matrix of SigniSite z-scores, corrected for multiple testing and truncated to significant positions
#'
#' @examples
#' signisite_zmat(file = 'path/to/my/msa_file.fsa')
#'
#' @export
signisite_zmat = function(file, method = 'bonferroni', alpha = 0.05){

  # Step 1 - Read the multiple sequence alignment file
  msa = read_fasta(file = file)

  # Step 2 - Retrieve the phenotype for each sequence
  msa$phenotype = get_values(fasta_header = msa$fasta_header)

  # Step 3 - Compute the SigniSite z-scores
  z_mat = get_signisite_zscores(sequences = msa$sequence, values = msa$phenotype)

  # Step 4 - Correct the matrix of z scores for multiple testing
  z_mat_adj = correct_z_matrix(z = z_mat, method = method)

  # Step 5 - Remove non-significant positions
  z_mat_adj_sp = rm_ns_positions(z = z_mat_adj, alpha = alpha)

  # Done
  return(z_mat_adj_sp)
}
