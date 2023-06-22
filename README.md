# Krill growth models intercomparison (Bahlburg et al. 2023)
This repository contains all scripts and analyses used to generate the results presented in

**Intercomparison of growth models for Antarctic krill (*Euphausia superba*): towards a generalised understanding (2023)**,
Dominik Bahlburg, Sally Thorpe, Bettina Meyer, Uta Berger, Eugene Murphy
https://www.biorxiv.org/content/10.1101/2023.05.09.540012v1

The required input- as well as output files (climatology datasets, simulation results) are available at:
https://zenodo.org/record/7560809

How to work with this repository:

Download/Clone the project, set the working directory to the project location.
The R environment (package versions + R version) used to develop the code is documented in the renv.lock-file contained in the repository.
Package versions can either be looked up in the renv.lock-file or the entire session can be restored by installing the renv-package and running
```
renv::restore()
```
after opening the project.

All scripts should then work as they are when the required input-files are in the correct locations (look at the file paths specified in the scripts)
