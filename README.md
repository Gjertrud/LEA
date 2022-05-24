# LEA
Matlab code for Likelihood Estimation of Affinity

This repository contains code that enables application of Likelihood Estimation of Affinity (LEA), written by Gjertrud Louise Laurell, Neurobiology Research Unit, 2022. LEA is based in Likelihood Estimation of Occupancy, developed by Dr. Martin Schain, and available at  https://github.com/martinschain/LEO.
To execute the code, you need MATLAB. The current code has been written and tested using MATLAB R2021a (version 9.10). 
The LEA function is found in lea.m. You have the option of providing your own initial parameter guesses. If you choose not to, initial parameter guesses will be evaluated using lea_getStartVals.m.
For estimation of the covariance structure of your data, we recommend that you download the code that calculates a non-linear shrinkage of the covariance matrix using test-retest data. That code is downloadable from http://www.econ.uzh.ch/en/people/faculty/wolf/publications.html#9 under the heading: Ledoit O. and Wolf, M. (2017). Numerical implementation of the QuEST function. Computational Statistics & Data Analysis 115, 199-223. 
For comments or questions regarding the LEA method or this code, please contact Gjertrud Louise Laurell.
