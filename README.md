# AccuCor - Isotope Natural Abundance Correction

[![R-CMD-check](https://github.com/XiaoyangSu/AccuCor/actions/workflows/r_standard_check.yml/badge.svg)](https://github.com/XiaoyangSu/AccuCor/actions/workflows/r_standard_check.yml)

AccuCor is an isotope natural abundance correction algorithm that is needed
especially for high resolution mass spectrometers. AccuCor supports 13C, 2H and
15N isotope corretion.

AccuCor accepts Excel (`.xls` and `.xlsx`), comma-separated text (`.csv`), or
tab-separated test (`.tsv`) files.


## Installation
```R
install.packages("devtools")
library(devtools)
devtools::install_github("XiaoyangSu/AccuCor")
```

## Quickstart

```R
library(accucor)

# Input file (example file included)
# Or use your own: carbon_input_file <- "/path/to/my/datafile.csv"

carbon_input_file <- system.file("extdata", "C_Sample_Input_Simple.csv", package = "accucor")


# Output is written to [input_file]_corrected.xlsx by default
# Be sure to specify the appropriate resolution.
# For Exactive, the resolution is 100000, defined at 200 Mw

carbon_corrected <- natural_abundance_correction(
  path = carbon_input_file,
  resolution = 100000)


# The results are also returned as a named list of dataframes for further processing in R
# "Original", "Corrected", "Normalized", "PoolBeforeDF", "PoolAfterDF"

carbon_corrected


# Purity is set to 0.99 for C and N, and 0.98 for H
# Be sure to specify purity if your samples differ
carbon_corrected <- natural_abundance_correction(
  path = carbon_input_file,
  resolution = 100000,
  purity = 0.97)
```

## Introduction

### Input files

Two types of input are accepted:
*   *Simple table* including the compound, formula, isotope label, and values for
    the samples

*   *Classic Maven* copy/paste along with a separate compound database (e.g.
    `knowns.csv`).

#### Simple table input must contain the following columns:

*   **Compound** - The name of the compound.

*   **Formula** - The chemical formula for the compound.

*   **IsotopeLabel** - The isotope labeling. Currently `C12 PARENT` for the
    parent (unlabeled) peak group and `C13-label-#` for the labeled peak groups.
    `D-label-#` and `N15-label-` are also recognized for Deuterium and Nitrogen
    labeled data, respectively.

*   **Samples** - One column for the intensity data for each sample.

*   **metaGroupId** - *OPTIONAL* - Used to support multiple peak groups per
    compound. Exported from El-MAVEN (version 0.4.0 and greater).

Additional columns commonly exported from El-MAVEN will be automatically
ignored.

*NOTE* - Previous versions required the compound database (`KNOWNS.csv`) as
well as the output data. This is no longer necessary, as the formula is
included the exported data from El-MAVEN.


#### Classic Maven input

*   Data input file - one block per compound. The first row of a block is the
    compound name followed by sample names. The remaining rows in the block
    contain the isotople label (*e.g.* `C12 PARENT`, `C13-label-1`) in the first
    column followed by intensities for each sample.

*   Compound database file - A comma separated (`.csv`) file with columns for the
    compound name (compound) and molecular formula (formula).

#### Input file formats

Input and output files may be Excel (`.xls` and `.xlsx`), comma-separated text
(`.csv`), or tab-separated test (`.tsv`) files. The type is automatically
determined by file extension.


### Isotope support

The isotope for correction is automatically determined by the contents of the
`isotopeLabel` column in your input data. AccuCor supports 13C, 2H and
15N isotope corretion. Only single labeled experiments are supported. Be sure
to turn off isotope detection for other isotopes (see
[El-MAVEN isotope detection documentation](https://github.com/ElucidataInc/ElMaven/wiki/Labeled-LCMS-Workflow#isotope-detection)).


### Citation
If you use this software in your research, please cite the following paper
(also see `citation("accucor")`):

Su X, Lu W and Rabinowitz J (2017). "Metabolite Spectral Accuracy on Orbitraps." *Analytical Chemistry*, *89*(11), pp. 5940-5948. doi:
10.1021/acs.analchem.7b00396 (URL: [http://doi.org/10.1021/acs.analchem.7b0039](http://doi.org/10.1021/acs.analchem.7b0039)),
PMID: 28471646, R package version 0.2.3 (2018).

A BibTeX entry for LaTeX users is

```
  @Article{,
    title = {Metabolite Spectral Accuracy on Orbitraps},
    author = {Xiaoyang Su and Wenyun Lu and Joshua D. Rabinowitz},
    journal = {Analytical Chemistry},
    doi = {10.1021/acs.analchem.7b00396},
    volume = {89},
    number = {11},
    pages = {5940-5948},
    year = {2017},
    note = {PMID: 28471646, R package version 0.2.3 (2018)},
    url = {https://doi.org/10.1021/acs.analchem.7b00396},
    eprint = {https://doi.org/10.1021/acs.analchem.7b00396},
  }
```
