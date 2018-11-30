#!/bin/bash -l

################################################################################

#Set the variables
MONTEPYTHONPATH="/Users/jfarr/Programs/montepython/"
N=1000000
OUTPUT="bao_chains_full_1000000"
PLOTDIR="reproduced_addison18_plots_1000000"

################################################################################

#Copy all of the relevant likelihoodis and data to montepython
cp -r likelihoods/* /$MONTEPYTHONPATH/montepython/likelihoods/
cp -r data/* /$MONTEPYTHONPATH/data

#Start python 2.7 environment
source activate py27

#Install the module 
python setup.py install --user

#Run the chains
python ${MONTEPYTHONPATH}/montepython/MontePython.py run -o ${OUTPUT}/bao_gal_lya_omegab -p input/bao.param -N $N --silent --conf ${MONTEPYTHONPATH}/default.conf
python ${MONTEPYTHONPATH}/montepython/MontePython.py run -o ${OUTPUT}/bao_lya_omegab -p input/bao_lya_omegab.param -N $N --silent --conf ${MONTEPYTHONPATH}/default.conf
python ${MONTEPYTHONPATH}/montepython/MontePython.py run -o ${OUTPUT}/bao_gal_omegab -p input/bao_gal_omegab.param -N $N --silent --conf ${MONTEPYTHONPATH}/default.conf
python ${MONTEPYTHONPATH}/montepython/MontePython.py run -o ${OUTPUT}/bao_gal_lya -p input/bao_gal_lya.param -N $N --silent --conf ${MONTEPYTHONPATH}/default.conf
python ${MONTEPYTHONPATH}/montepython/MontePython.py run -o ${OUTPUT}/bao_lya -p input/bao_lya.param -N $N --silent --conf ${MONTEPYTHONPATH}/default.conf
python ${MONTEPYTHONPATH}/montepython/MontePython.py run -o ${OUTPUT}/bao_gal -p input/bao_gal.param -N $N --silent --conf ${MONTEPYTHONPATH}/default.conf

#Make the plots
python ${MONTEPYTHONPATH}/montepython/MontePython.py info ${OUTPUT}/bao_gal/ ${OUTPUT}/bao_lya/ ${OUTPUT}/bao_gal_lya/ --no-mean --extra plot_files/plot_Om_vs_H0rd.txt 
python ${MONTEPYTHONPATH}/montepython/MontePython.py info ${OUTPUT}/bao_gal_omegab/ ${OUTPUT}/bao_lya_omegab/ ${OUTPUT}/bao_gal_lya_omegab/ --no-mean --extra plot_files/plot_Om_vs_H0.txt

#Move the plots to a more clear location
mkdir ${PLOTDIR}
cp ${OUTPUT}/bao_gal/plots/bao_gal-vs-bao_lya-vs-bao_gal_lya_triangle.pdf ${PLOTDIR}/plot_3.pdf
cp ${OUTPUT}/bao_gal_omegab/plots/bao_gal_omegab-vs-bao_lya_omegab-vs-bao_gal_lya_omegab_triangle.pdf ${PLOTDIR}/plot_4.pdf
open ${PLOTDIR}/*.pdf
