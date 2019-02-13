# ---------------------------------------------------------------------------- #
# SigniSite, an RPackage for easy geno to phenotype analysis and visualisation #
#     Copyright (C) 2019 Leon Eyrich Jessen                                    #
# ---------------------------------------------------------------------------- #
#' Extract white space separated end-placed values from a set of FASTA headers
#'
#' @param fasta_header A vector of FASTA headers formatted like so: ">id info VALUE"
#'
#' @return A numerical vector of the extracted values
#'
#' @examples
#' get_values(fasta_header = ALIGNMENT$fasta_header)
#'
#' @export
get_values = function(fasta_header){

  split_lst = strsplit(x = fasta_header, split = "\\s+")
  values = unlist(lapply(split_lst, function(x){ return( x[length(x)] ) }))
  return( as.numeric(values) )

}
