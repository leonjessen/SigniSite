# ---------------------------------------------------------------------------- #
# SigniSite, an RPackage for easy geno to phenotype analysis and visualisation #
#     Copyright (C) 2019 Leon Eyrich Jessen                                    #
# ---------------------------------------------------------------------------- #
# Set a data frame with names and abbreviations of the standard 20
# proteogenic amino acid residues
# Wed Feb 13 09:51:07 2019 ------------------------------
# Run as .set_aminoacids()
.set_aminoacids = function(){
  AMINOACIDS = data.frame(
    full      = c('Alanine', 'Arginine', 'Asparagine', 'Aspartate', 'Cysteine',
                  'Glutamine', 'Glutamate', 'Glycine', 'Histidine',
                  'Isoleucine', 'Leucine', 'Lysine', 'Methionine',
                  'Phenylalanine', 'Proline', 'Serine', 'Threonine',
                  'Tryptophan', 'Tyrosine', 'Valine'),
    three     = c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His',
                  'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp',
                  'Tyr', 'Val'),
    one       = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                  'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
    chemistry = c('Hydrophobic', 'Basic', 'Neutral', 'Acidic', 'Polar',
                  'Neutral', 'Acidic', 'Polar', 'Basic', 'Hydrophobic',
                  'Hydrophobic', 'Basic', 'Hydrophobic', 'Hydrophobic',
                  'Hydrophobic', 'Polar', 'Polar', 'Hydrophobic', 'Polar',
                  'Hydrophobic'),
    stringsAsFactors = FALSE
  )
  save(AMINOACIDS, file = "data/AMINOACIDS.RData")
  return(0)
}

# Create the main logo for SigniSite
# Thu Feb 14 08:55:56 2019 ------------------------------
# Run as .create_main_logo_png()
.create_main_logo_png = function(){

  # Set seed for reproducible logo generation
  set.seed(16020)

  # Define chars in logo
  res_chars = AMINOACIDS$one

  # Set matrix of negative residue scores
  m = 20
  n = 9
  res_mat = matrix(data = -abs(rnorm(m*n))/5, nrow = m, ncol = n)
  rownames(res_mat) = res_chars
  colnames(res_mat) = paste0('p_', 1:9)

  # Set scores for residues forming SIGNISITE
  res_mat['S','p_1'] = 11
  res_mat['I','p_2'] = 7
  res_mat['G','p_3'] = 7
  res_mat['N','p_4'] = 7
  res_mat['I','p_5'] = 7
  res_mat['S','p_6'] = 11
  res_mat['I','p_7'] = 7
  res_mat['T','p_8'] = 7
  res_mat['E','p_9'] = 7

  # Create logo
  main_logo = ggplot() +
    geom_logo(data = res_mat, method = 'custom') +
    theme_logo() +
    theme(legend.position  = 'none',
          axis.text.x      = element_blank(),
          axis.text.y      = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'transparent', colour = NA),
          plot.background  = element_rect(fill = 'transparent', colour = NA))

  # Save png
  ggsave(filename = 'man/figures/signisite_logo.png',
         plot = main_logo, device = 'png', width = 16, height = 16,
         units = 'cm', dpi = 144, bg = 'transparent')

  # Done
  return(0)
}
# Wed Feb 13 17:53:39 2019 ------------------------------
.set_example_data = function(){
  ALIGNMENT = read_fasta(file = 'data/signisite_alignment.fsa')
  save(ALIGNMENT, file = "data/ALIGNMENT.RData")
  return(0)
}

rm_zero_pos = function(pssm){
  pssm = pssm[-which(rowSums(abs(pssm))==0),]
  return(pssm)
}
rm_ns_res = function(pssm, alpha = 0.05){
  pssm[qnorm(alpha/2) < pssm & pssm < qnorm(1-alpha/2)] = 0
  return(pssm)
}
manhattan = function(x){
  s_mat = sign(x)
  p_mat = 2 * pnorm( -1 * abs( x ) )
  m_mat = -log10(p_mat) * s_mat
  return(m_mat)
}
