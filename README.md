# H0LyA
Repository documenting project work to constrain H0 from Lya BAO measurements.

MCMC chains are run using MontePython, the latest version of which can be downloaded from https://github.com/brinckmann/montepython_public/tree/3.0. 

Sample data, likelihoods and input parameter files are found in the relevant directories.

To install, use ``python setup.py install''.

``python scripts/plot_scans.py'' plots the BOSS chi2 scans. Options can be modified in the script (function described in py/chi2_scan.py).
