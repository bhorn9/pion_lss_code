# pion_lss_code
******************************************************************

Welcome to Pion LAgrangian for STructure In Cosmology! 

******************************************************************

The latest version of the code is available at
https://github.com/bhorn9/pion_lss_code

This code is for simulating the evolution of large scale structure in a Lambda CDM universe, using the pion field formalism of Nicolis and Vernizzi, discussed in hep-th/1406.0842 and hep-th/2604.xxxxx, which corresponds to the velocity potential $\vec{v} = \vec{\nabla} \pi$ of the cosmic matter fluid.  

Authors: Bart Horn and Bhavya Mishra, with additional input and testing from Lara Celik, Eamon McShane, and David Muqattash.  This work was supported in part by the National Science Foundation's Research at Undergraduate Institutions Program under research grant PHY-2210475.

Copyright: 2026 by Bart Horn and Bhavya Mishra.  This code may be used free of charge for research and educational purposes; however, we kindly request that you cite the corresponding paper hep-th/2604.xxxxx (full reference below) in publications and scientific applications. 

For further inquiries please contact: Bart Horn at bhorn01@manhattan.edu

******************************************************************

The code is written in Python version 3.11.9 and in C++ std = 17.  C++ typicially runs much more quicky but depending on your system Python may be more versatile at handling memory allocation at high resolution.  

The repository consists of the following programs and the linked inputs file.

inputs.txt
Cosmic and simulation input parameters can be adjusted here.  The order of parameters can be changed but not the names.  Remember to include decimals when declaring doubles if using C++.
-Omega_m, H0 are the matter fraction and Hubble parameter in the present day.
-z_init, z_final are the initial and final redshifts.
-boxsize is the size of the periodic 3D box.
-meshpoints is the resolution of the discrete spatial grid.  Runtime scales like meshpoints^3 log meshpoints.
-sigma is the scale of Gaussian smoothing for both initial and evolved conditions.
-cs2_0 and cv2_0 are the present day values of the dimensionless sound speed and viscosity EFT parameters.

Gaussian_ICs.py
generates Gaussian initial conditions for the pion field and outputs them as a text file.  The desired linear power spectrum can be calculated from a CLASS output file (given as a .dat file) or can be fit using either the Eisentein-Hu or the BBKS analytic fitting formula.

PLASTIC.py/PLASTIC++.cpp
evolves the pion field from the initial conditions until nonlinear structure begins to form, and produces a snapshot of pi(x), pi'(x) and delta(x) as a text file.  This step can be run either in Python or in C++.  PLASTIC++ requires having the FFTW3 library installed.  To compile in Ubuntu, use the command g++ PLASTIC++.cpp -std=c++17 -lfftw3 -lm

Plots.py
generates graphs of the power spectra and 3D fields for pi(x), pi'(x), and delta(x).

******************************************************************

References:

[1] Alberto Nicolis and Filippo Vernizzi, unpublished notes.

[2] "Soft-Pion Theorems for Large Scale Structure," Bart Horn, Lam Hui and Xiao Xiao, JCAP 09 (2014) 044, https://arxiv.org/hep-th/1406.0842

[3] "Effective Field Theory of a single Scalar Pion Field for Large Scale Structure," Lara Celik, Bart Horn, Bhavya Mishra and David Muqattash, https://arxiv.org/astro-ph-CO/2604.xxxxx
