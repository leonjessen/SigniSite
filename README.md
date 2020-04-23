
# SigniSite <a href="http://www.cbs.dtu.dk/services/SigniSite/"><img src="man/figures/signisite_hex_logo.png" align="right" height="200" /></a>

By popular demand, this is an `R` implementation of the [SigniSite 2.1
Server: Residue level genotype phenotype correlation in protein multiple
sequence alignments](http://www.cbs.dtu.dk/services/SigniSite/). For
details on the method, please see the original research paper
[SigniSite: Identification of residue-level genotype-phenotype
correlations in protein multiple sequence
alignments](https://academic.oup.com/nar/article/41/W1/W286/1111259)
(See [Citation](#Citation))

## Overview

SigniSite performs residue level genotype phenotype correlation in
protein multiple sequence alignments by identifying amino acid residues
significantly associated with the phenotype of the data set. Input is a
protein multiple sequence alignment in FASTA format. The phenotype is
represented by a real-valued numerical value placed last in the
identifier of each sequence (white-space separated).

## License

SigniSite is free to use for academic users as [described in the
license](LICENSE); other users are requested to contact Software Package
Manager at <health-software@dtu.dk>.

## Installation

``` r
# The development version from GitHub:
install.packages("devtools")
devtools::install_github("leonjessen/SigniSite")
```

## Usage

Once installed, `SigniSite` is loaded like so

``` r
library("SigniSite")
```

There are two options for using `SigniSite`, either simply call
`SigniSite` with default options like so:

``` r
# Create the matrix of SigniSite z-scores
z = signisite_zmat(file = 'path/to/my/msa_file.fsa', method = 'bonferroni', alpha = 0.05)

# Visualise z-score matrix as a logo plot
plot_signisite_logo(z)
```

Or go through the following workflow:

1.  `read_fasta()`: Reads a multiple sequence alignment in FASTA format
2.  `get_values()`: Retrieve the sequence associated values from FASTA
    headers
3.  `get_signisite_zscores()`: Compute the SigniSite z-scores
4.  `correct_z_matrix()`: Correct the z-scores for multiple testing
5.  `rm_ns_positions()`: Remove any positions, where no residues are
    significantly associated with the sequence associated values
6.  `plot_signisite_logo()`: Visualise the results

The following is a walk through of the above workflow:

#### Step 1 - Read the multiple sequence alignment file

``` r
# Replace 'data/signisite_alignment.fsa' with your path to your file
msa = read_fasta(file = 'data/signisite_alignment.fsa')
```

If you do not have your own alignment file, but simply want to test the
`SigniSite` package, then you can use the build in example alignment
like so:

``` r
msa = ALIGNMENT
```

The multiple sequence alignment looks like so:

``` r
head(msa)
```

    ##       fasta_header
    ## 1  >144985_01 58.0
    ## 2 >144986_01 123.3
    ## 3 >144987_01 117.6
    ## 4 >144988_01 104.5
    ## 5  >144989_01 13.3
    ## 6  >144990_01 67.1
    ##                                                                                              sequence
    ## 1 PQITLWQRPFVTVKIGGQLKEALIDTGADDTVFEEINLPGRWKPKIIGGIGGFMKVREYDQIPVEICGHKAIGTVLVGPTPVNVIGRNLLTQIGCTLNF
    ## 2 PQITLWQRPFITVKIEGQLKEALIDTGADDTIFEEINLPGRWKPKIVGGIGGFMKVREYDQIPVEICGHKAIGTVLVGPTPVDVIGRNLLTQIGCTLNF
    ## 3 PQITLWQRPIVTVKIGGQLKEALLDTGADDTIFEELNLPGRWKPKIIGGIGGFIKVRQYDQVLVEICGHKAIGTVLVGPTPVDVIGRNLMTQIGCTLNF
    ## 4 PQITLWQRPIVTVKIGGQLKEALLDTGADDTIFEELNLPGRWKPKIIGGIGGFIKVRQYDQVLVEICGHKAIGTVVVGPTPVDVIGRNLMTQIGCTLNF
    ## 5 PQITLWQRPFITVKIGGQLKEALLDTGADDTVFEEMNLPGRWKPKMIGGIGGFLKVREYDQIPIEICGHKAIGTVLVGPTPVNVIGRNLLTQIGCTLNF
    ## 6 PQITLWQRPFITVKIGGQLKEALLDTGADDTIFEEMNLPGRWKPKMIGGIGGFLKVREYDQIPIEICGHKAIGPVLVGPTPVNVIGRNLLTQIGCTLNF

#### Step 2 - Retrieve the phenotype for each sequence

``` r
msa$phenotype = get_values(fasta_header = msa$fasta_header)
```

The multiple sequence alignment now looks like so:

``` r
head(msa)
```

    ##       fasta_header
    ## 1  >144985_01 58.0
    ## 2 >144986_01 123.3
    ## 3 >144987_01 117.6
    ## 4 >144988_01 104.5
    ## 5  >144989_01 13.3
    ## 6  >144990_01 67.1
    ##                                                                                              sequence
    ## 1 PQITLWQRPFVTVKIGGQLKEALIDTGADDTVFEEINLPGRWKPKIIGGIGGFMKVREYDQIPVEICGHKAIGTVLVGPTPVNVIGRNLLTQIGCTLNF
    ## 2 PQITLWQRPFITVKIEGQLKEALIDTGADDTIFEEINLPGRWKPKIVGGIGGFMKVREYDQIPVEICGHKAIGTVLVGPTPVDVIGRNLLTQIGCTLNF
    ## 3 PQITLWQRPIVTVKIGGQLKEALLDTGADDTIFEELNLPGRWKPKIIGGIGGFIKVRQYDQVLVEICGHKAIGTVLVGPTPVDVIGRNLMTQIGCTLNF
    ## 4 PQITLWQRPIVTVKIGGQLKEALLDTGADDTIFEELNLPGRWKPKIIGGIGGFIKVRQYDQVLVEICGHKAIGTVVVGPTPVDVIGRNLMTQIGCTLNF
    ## 5 PQITLWQRPFITVKIGGQLKEALLDTGADDTVFEEMNLPGRWKPKMIGGIGGFLKVREYDQIPIEICGHKAIGTVLVGPTPVNVIGRNLLTQIGCTLNF
    ## 6 PQITLWQRPFITVKIGGQLKEALLDTGADDTIFEEMNLPGRWKPKMIGGIGGFLKVREYDQIPIEICGHKAIGPVLVGPTPVNVIGRNLLTQIGCTLNF
    ##   phenotype
    ## 1      58.0
    ## 2     123.3
    ## 3     117.6
    ## 4     104.5
    ## 5      13.3
    ## 6      67.1

#### Step 3 - Compute the SigniSite z-scores

``` r
z_mat = get_signisite_zscores(sequences = msa$sequence, values = msa$phenotype)
head(z_mat)
```

    ##      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X  -
    ## 1_P NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0 NA NA NA NA NA NA NA
    ## 2_Q NA NA NA NA NA  0 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
    ## 3_I NA NA NA NA NA NA NA NA NA  0 NA NA NA NA NA NA NA NA NA NA NA NA
    ## 4_T NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0 NA NA NA NA NA
    ## 5_L NA NA NA NA NA NA NA NA NA NA  0 NA NA NA NA NA NA NA NA NA NA NA
    ## 6_W NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0 NA NA NA NA

#### Step 4 - Correct the matrix of z scores for multiple testing

``` r
z_mat_adj = correct_z_matrix(z = z_mat)
head(z_mat_adj)
```

    ##      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X  -
    ## 1_P NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0 NA NA NA NA NA NA NA
    ## 2_Q NA NA NA NA NA  0 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
    ## 3_I NA NA NA NA NA NA NA NA NA  0 NA NA NA NA NA NA NA NA NA NA NA NA
    ## 4_T NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0 NA NA NA NA NA
    ## 5_L NA NA NA NA NA NA NA NA NA NA  0 NA NA NA NA NA NA NA NA NA NA NA
    ## 6_W NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0 NA NA NA NA

#### Step 5 - Remove non-significant positions

``` r
z_mat_adj_sp = rm_ns_positions(z = z_mat_adj)
head(z_mat_adj_sp)
```

    ##       A  R  N  D  C  Q  E  G  H        I         L        K         M
    ## 10_L NA NA NA NA NA NA NA NA NA 2.055025 -5.318753       NA        NA
    ## 24_L NA NA NA NA NA NA NA NA NA 3.003257 -3.003257       NA        NA
    ## 32_V NA NA NA NA NA NA NA NA NA 3.803664        NA       NA        NA
    ## 33_L NA NA NA NA NA NA NA NA NA       NA -3.645199       NA        NA
    ## 43_K NA NA NA NA NA NA NA NA NA       NA        NA -1.96658        NA
    ## 46_M NA NA NA NA NA NA NA NA NA 1.295326  2.006375       NA -3.869101
    ##             F  P  S       T  W  Y         V  X  -
    ## 10_L 1.018376 NA NA      NA NA NA  0.000000 NA NA
    ## 24_L       NA NA NA      NA NA NA        NA NA NA
    ## 32_V       NA NA NA      NA NA NA -3.803664 NA NA
    ## 33_L 3.645199 NA NA      NA NA NA        NA NA NA
    ## 43_K       NA NA NA 1.96658 NA NA        NA NA NA
    ## 46_M       NA NA NA      NA NA NA        NA NA NA

#### Step 6 - Visualise the result

``` r
library('ggplot2')
library('ggseqlogo')
plot_signisite_logo(z = z_mat_adj_sp)
```

<img src="README_files/figure-gfm/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

\#\#Citation When using SigniSite for publications, please cite the
original reserach paper:

``` r
citation("SigniSite")
```

    ## 
    ## To cite SigniSite in publications use:
    ## 
    ##   Leon Eyrich Jessen, Ilka Hoof, Ole Lund and Morten Nielsen
    ##   (2013). SigniSite: Identification of residue-level
    ##   genotype-phenotype correlations in protein multiple sequence
    ##   alignments. Nucleic Acids Research, 41(W1), W286-W291. URL
    ##   https://doi.org/10.1093/nar/gkt497.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Article{,
    ##     title = {{SigniSite}: Identification of residue-level genotype-phenotype correlations in protein multiple sequence alignments},
    ##     author = {Leon Eyrich Jessen and Ilka Hoof and Ole Lund and Morten Nielsen},
    ##     journal = {Nucleic Acids Research},
    ##     year = {2013},
    ##     volume = {41},
    ##     number = {W1},
    ##     pages = {W286-W291},
    ##     url = {https://doi.org/10.1093/nar/gkt497},
    ##   }

Furthermore, `SigniSite` uses `ggseqlogo` for visualisation, so please
also cite:

[Wagih O. ggseqlogo: a versatile R package for drawing sequence logos.
Bioinformatics. 2017 Nov 15;33(22):3645-3647.
doi: 10.1093/bioinformatics/btx469.](https://academic.oup.com/bioinformatics/article/33/22/3645/3980251)
