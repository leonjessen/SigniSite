# ---------------------------------------------------------------------------- #
# SigniSite, an RPackage for easy geno to phenotype analysis and visualisation #
#     Copyright (C) 2019 Leon Eyrich Jessen                                    #
# ---------------------------------------------------------------------------- #
#' Plot the results of a SigniSite analysis as a sequence logo
#'
#' @param z A matrix of z-scores
#'
#' @return A ggplot object
#'
#' @examples
#' z_mat = get_signisite_zscores(ALIGNMENT$sequence, get_values(ALIGNMENT$fasta_header))
#' z_mat_sp = rm_ns_positions(z_mat, alpha = 0.05)
#' plot_signisite_logo(z_mat_sp)
#'
#' @export
plot_signisite_logo = function(z){
  logo_plot = ggplot() +
    geom_logo(data = t(z), method = 'custom', seq_type = 'aa') +
    geom_hline(yintercept = 0) +
    labs(x = "MSA position and consensus sequence", y = "z-score") +
    scale_x_continuous(breaks = 1:nrow(z), labels = rownames(z)) +
    theme_logo() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  return(logo_plot)
}
