# ---------------------------------------------------------------------------- #
# SigniSite, an RPackage for easy geno to phenotype analysis and visualisation #
#     Copyright (C) 2019 Leon Eyrich Jessen                                    #
# ---------------------------------------------------------------------------- #
#' Compute the SigniSite z-scores
#'
#' @param sequences A set of sequences forming a multiple sequence alignment
#' @param values The for each sequence associated numerical value
#'
#' @return A matrix of z-scores, where each score quantifies the strength of the association between a given amino acid residue at a given position and the phenotype
#'
#' @examples
#' get_signisite_zscores(ALIGNMENT$sequence, get_values(ALIGNMENT$fasta_header))
#'
#' @export
get_signisite_zscores = function(sequences, values){

  # Check if sequences are a character vector
  if( !is.character(sequences) | !is.vector(sequences) ){
    stop("SigniSite: Sequences are not a character vector")
  }
  # Check if sequences are a numeric vector
  if( !is.numeric(values) | !is.vector(values) ){
    stop("SigniSite: Sequences are not a character vector")
  }
  # Check if all sequences have equal length
  if( length(unique(nchar(sequences))) > 1 ){
    stop("SigniSite: Not all sequence have the same length. The sequences are required to be pre-aligned")
  }
  # Check if the number of sequences and associated values are equal
  if( length(sequences) != length(values) ){
    stop("SigniSite: The number of sequences must equal the number of values")
  }

  # Define amino acid residue alphabet
  res_chars = c(AMINOACIDS$one, 'X', '-')

  # Gaps are represented as '*' rather than '-'
  gsub(pattern = "-", replacement = "*", x = sequences)

  # Get the total number of sequences in the MSA
  n_seqs = length(sequences)

  # Get the consensus sequence for the MSA
  cons_seq = get_consensus(sequences)

  # Create sequence matrix:
  # Rows: Each of the sequences
  # Columns: MSA positions
  seq_mat = do.call(rbind, strsplit(x = ALIGNMENT$sequence, split = ''))

  # Transpose the sequence matrix to work per position
  seq_mat = t(seq_mat)

  # Create empty matrix for holding final results
  # Rows: Positions
  # Columns: Amino acid res_chars
  z_mat = matrix(NA, nrow = nrow(seq_mat), ncol = length(res_chars))
  rownames(z_mat) = paste(seq(from = 1, to = nrow(z_mat)), cons_seq, sep = '_')
  colnames(z_mat) = res_chars

  # Pre-calculate statistical parameters
  ranks  = rank(-values)    # This needs to be up to the user high/low==strong
  mu_exp = (n_seqs + 1) / 2 # == mean(1:n_seqs)
  tcf    = tie_c_factor(x = values)
  sigmas = sapply(1:(n_seqs-1),
                  function(n_aa){
                    sqrt( ((n_seqs-n_aa) * (n_seqs+1) * tcf) / (12*n_aa) )
                  })

  # Compute SigniSite z-scores
  for( i in 1:nrow(z_mat) ){
    for( j in 1:ncol(z_mat) ){
      indices = which(seq_mat[i,] == res_chars[j])
      n_aa    = length(indices)
      if( 0 < n_aa & n_aa < n_seqs ){
        mu_obs     = mean(ranks[indices])
        sigma_exp  = sigmas[n_aa]
        z_mat[i,j] = (mu_exp - mu_obs) / sigma_exp
      } else if( n_aa == n_seqs ){
        z_mat[i,j] = 0
      }
    }
  }

  # Done, return matrix of SigniSite z-scores
  return(z_mat)
}
