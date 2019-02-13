# ---------------------------------------------------------------------------- #
# SigniSite, an RPackage for easy geno to phenotype analysis and visualisation #
#     Copyright (C) 2019 Leon Eyrich Jessen                                    #
# ---------------------------------------------------------------------------- #
#' Get the consensus sequence from a multiple sequence alignment
#'
#' For each position, the most frequently observed amino acid residue forms the consensus sequence
#'
#' @param sequences A set of sequences forming a multiple sequence alignment
#'
#' @return A single sequence, representing the consensus sequence
#'
#' @examples
#' get_signisite_zscores(ALIGNMENT$sequence, get_values(ALIGNMENT$fasta_header))
#'
#' @export
get_consensus = function(sequences){
  seq_mat = do.call(rbind, strsplit(x = ALIGNMENT$sequence, split = ''))
  cons_seq = apply(seq_mat, 2, function( x_j ){
    x_j_counts = table(x_j)
    max_count  = max(x_j_counts)
    max_pos    = which(x_j_counts == max_count)
    cons_res   = names(max_pos)
    if( length( cons_res ) > 1 ){
      cons_res = paste0('[', paste0(cons_res, collapse = ''), ']', collapse = '')
    }
    return( cons_res )
  })
  return(cons_seq)
}
