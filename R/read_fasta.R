# ---------------------------------------------------------------------------- #
# SigniSite, an RPackage for easy geno to phenotype analysis and visualisation #
#     Copyright (C) 2019 Leon Eyrich Jessen                                    #
# ---------------------------------------------------------------------------- #
#' Read a fasta file
#'
#' @param file A multiple sequence alignment in FASTA format
#'
#' @return A data frame of FASTA headers and sequences
#'
#' @examples
#' read_fasta(file = 'path/to/my/msa_file.fsa')
#'
#' @export
read_fasta = function(file){

  # Check input
  if( !is.character( file ) ){
    stop("'file' has to be a string specifying a file name")
  }
  if( !file.exists( file ) ){
    stop(paste("Unable to read file",file))
  }

  # Read lines from file
  lines = readLines(con = file)

  # Remove any comments
  comments = grep(pattern = "^#", x = lines)
  if( length( comments > 0 ) ){
    lines = lines[-comments]
  }

  # Remove any empty lines
  empty_lines = grep(pattern = "^$", x = lines)
  if( length( empty_lines ) > 0 ){
    lines = lines[-empty_lines]
  }

  # Set variables
  n_seqs    = length(grep("^>", lines)) # The number of lines starting with ">"
  headers   = rep(NA, n_seqs)
  sequences = rep('', n_seqs)

  # Iterate over lines and fill containers
  entry_no = 1
  for( line in lines ){
    if( grepl(pattern = "^>", x = line) ){
      headers[entry_no] = line
      entry_no = entry_no + 1
    } else {
      sequences[entry_no - 1] = paste0(sequences[entry_no - 1], line)
    }
  }

  # Check sequences for any non-standard characters
  sequences = toupper(x = sequences)
  aa_pattern = "[^ARNDCQEGHILKMFPSTWYVX-]"
  if( any(grepl(pattern = aa_pattern, x = sequences)) ){
    msg = paste("Non standard character in sequences.",
                "Only 'ARNDCQEGHILKMFPSTWYVX-' allowed!",
                "Check your alignment!",
                "Proceeding by replacing with 'X'")
    warning(msg)
    sequences = gsub(pattern = "[^ARNDCQEGHILKMFPSTWYVX-]",
                     replacement = 'X', x = sequences)
  }

  # Check if the sequences are aligned
  if( length(unique(nchar(ALIGNMENT$sequence))) > 1 ){
    warning("SigniSite: The read FASTA files contain sequences of different length. To perform a SigniSite analysis, sequences must be pre-aligned")
  }

  # Return output data frame
  return(data.frame(fasta_header = headers, sequence = sequences,
                    stringsAsFactors = FALSE))
}
