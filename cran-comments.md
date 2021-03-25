## Test environments
Testing is run using GitHub actions with following matrix:
* os: windows-latest, r: 'release'
* os: macOS-latest, r: 'release'
* os: ubuntu-20.04, r: 'release', rspm: "https://packagemanager.rstudio.com/cran/__linux__/focal/latest"
* os: ubuntu-20.04, r: 'devel', rspm: "https://packagemanager.rstudio.com/cran/__linux__/focal/latest"

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies
As this is a new package, there are no downstream dependencies yet.
