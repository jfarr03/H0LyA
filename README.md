# H0LyA
Repository documenting project work to constrain H0 from Lya BAO measurements.

MCMC chains are run using MontePython, the latest version of which can be downloaded from https://github.com/brinckmann/montepython_public/tree/3.0. 

Sample data, likelihoods and input parameter files are found in the relevant directories.

Some simple commands are:
 - To install, use `python setup.py install`.
 - To plot the BOSS chi2 scans use `python scripts/plot_scans.py`. Options can be modified in the script (function described in py/chi2_scan.py).
 - To use the module and MontePython to reproduce plots from Addison et al. 2018 run `source scripts/reproduce_addison18_plots.sh`.
