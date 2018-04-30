# AccuCor - Isotope Natural Abundance Correction

[![Travis build status](https://travis-ci.org/lparsons/accucor.svg?branch=master)](https://travis-ci.org/lparsons/accucor)

AccuCor is an isotope natural abundance correction algorithm that is needed 
especially for high resolution mass spectrometers. AccuCor supports 13C, 2H and
15N isotope corretion. 

AccuCor accepts Excel (`.xls` and `.xlsx`), comma-separated text (`.csv`), or
tab-separated test (`.tsv`) files.


## Installation
```R
install.packages("devtools")
devtools::install_github("lparsons/accucor")
```

## Quickstart

```R
library(accucor)

# Please make sure these parameters are accurate.
resolution <- 100000  # For Exactive, the Resolution is 100000, defined at Mw 200
resolution_defined_at <- 200

# Input file (example file included)
carbon_input_file <- system.file("extdata", "C_Sample_Input_Simple.csv", package = "accucor")

# Output is written to [input_file]_corrected.xlsx by default
carbon_corrected <- natural_abundance_correction(
  path = carbon_input_file,
  resolution = resolution, 
  resolution_defined_at = resolution_defined_at)

# The results are also returned as a named list of dataframes for further processing in R
# "Original", "Corrected", "Normalized", "PoolBeforeDF", "PoolAfterDF"
carbon_corrected
```

## Introduction

### Input files

Two types of input are accepted: 
* Simple table including the compound, formula, isotope label, and values for the samples
* Classic Maven copy/paste along with a separate compound database (e.g. `knowns.csv`).

#### Simple table input must contain the following columns:

* **Compound** - The name of the compound.

  *NOTE* - Multiple peaks per compound is not yet supported (see 
  [Issue #8](https://github.com/lparsons/Isotope-Natural-Abundance-Correction/issues/8)).
  
* **Formula** - The chemical formula for the compound.

   *NOTE* - Previous versions required the compound database (`KNOWNS.csv`) as
   well as the output data. This is no longer necessary, as the formula is
   included the exported data from El-MAVEN.
   
* **IsotopeLabel** - The isotope labeling. Currently `C12 PARENT` for the parent
  (unlabeled) peak group and `C13-label-#` for the labeled peak groups.
  `D-label-#` and `N15-label-` are also recognized for Deuterium and Nitrogen
  labeled data, respectively.
  
* **Samples** - One column for the intensity data for each sample.

Additional columns commonly exported from El-MAVEN will be automatically
ignored.


#### Classic Maven input

* Data input file - one block per compound. The first row of a block is the compound name followed by sample names. The remaining rows in the block contain the isotople label (*e.g.* `C12 PARENT`, `C13-label-1`) in the first column followed by intensities for each sample.

* Compound database file - A comma separated (`.csv`) file with columns for the compound name, molecular formula, and molecular weight.

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
