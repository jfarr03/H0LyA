# H0LyA
Repository documenting project work to constrain H0 from Lya BAO measurements. Rather than using a simple, Gaussian error approximation (as is common in the literature), the full BOSS chi2 surface is used, providing a more accurate representation of the errors.

MCMC chains are run using MontePython, the latest version of which can be downloaded from https://github.com/brinckmann/montepython_public/tree/3.0. 

Sample data, likelihoods and input parameter files are found in the relevant directories.

Some simple commands are:
 - To install, use `python setup.py install`.
 - To plot the BOSS chi2 scans use `python scripts/plot_scans.py`. Options can be modified in the script (function described in py/chi2_scan.py).
 - To use the module and MontePython to reproduce plots from Addison et al. 2018 run `source scripts/reproduce_addison18_plots.sh`.

To do:
 - Modify the plotting of `chi2_scan.py` to include covariance between the two distance measures when plotting the Gaussian approximation.
 - Explain the "wobbly edges" of the reproduced Addison et al. 2018 plots.
 - Add the likelihood to MontePython.
 
