#' Identification of residue-level genotype-phenotype correlations in protein multiple sequence alignments
#'
#' SigniSite performs residue level genotype phenotype correlation in 
#' protein multiple sequence alignments by identifying amino acid residues 
#' significantly associated with the phenotype of the data set. 
#' Input is a protein multiple sequence alignment in FASTA format. The 
#' phenotype is represented by a real-valued numerical value placed last in 
#' the identifier of each sequence (white-space separated).
#'
#' @name SigniSite-package
#' @aliases SigniSite
#' @docType package
#' @author
#' Maintainer: Leon Eyrich Jessen <ljess@dtu.dk>
#' @examples
#' library("SigniSite")
#' 
#' # Create the matrix of SigniSite z-scores
#' z = signisite_zmat(file = system.file(package = "SigniSite", "/signisite_alignment.fsa"), 
#'                    method = 'bonferroni', alpha = 0.05)
#' 
#' # Visualise z-score matrix as a logo plot
#' plot_signisite_logo(z)
#' @references Leon Eyrich Jessen, Ilka Hoof, Ole Lund and Morten Nielsen
##   (2013). SigniSite: Identification of residue-level
##   genotype-phenotype correlations in protein multiple sequence
##   alignments. Nucleic Acids Research, 41(W1), W286-W291. URL
##   \url{https://doi.org/10.1093/nar/gkt497}
#' @keywords package
#' @importFrom stats p.adjust pnorm qnorm rnorm
#' @importFrom ggplot2 ggplot ggsave geom_hline element_blank element_rect element_text labs scale_x_continuous theme
#' @importFrom ggseqlogo geom_logo theme_logo
NULL
