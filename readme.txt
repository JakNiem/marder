marder -- a simple, non-parallel molecular dynamics program. Written by ???.
Implementation of spherical corrections by Isabel Nitzke.
Modifications by Jakob Niemann 2024.



***** General Information:

- The program is not designed for mixtures.

- This version can only handle droplets. For bubbles, a different center of mass alignment is required and some inequalities during the determination of the density profile have to be changed.

- It was tested with compiler: gfortran -fdefault-integer-8 CMyQI_NPT_RDF.f90
The option -fdefault-integer-8 is always required.


***** Input: 

- Information about the state point and fluid are specified in the input file like '1CLJ.dat'.

- The initial configuration must be provided in a separate file, like '1clj_4842.dat', which is referred to in the input file. 

- The program is called with: ./a.out 1CLJ


***** Spherical LRC:

- The density profile is determined in 'RhoProfile.f90'. At the same time, the homogeneous corrections are computed.

- In each time step, the current density is saved and particles receive a correction term according to their current location. Every 1000 steps, tanh density profile and corrections are computed anew. This is all included in 'ShellInit.f90'.

- The corrections per shell are computed in 'ShellCorr.f90'.

- The decomposition of the virial is computed in 'Force_ij.f90'.

  
***** Output:

- The parameters obtained in the fit of the density profile are listed in *_tanh.csv.

- Averages of the normal and tangential virial and pressure are found in *_Pres.csv.

- The average density and current correction terms per shell are found in *_Shells.csv.

- The file *_POS.dat can be used as an input to restart the simulation.


- file *_Gamma.csv containes data concerning density, pressure and surface tension


notes:

execute for debug:
    gfortran -g -o marder -fdefault-integer-8 CMyQI_NPT_RDF.f90

execute for production:
    gfortran -o marder -fdefault-integer-8 CMyQI_NPT_RDF.f90
