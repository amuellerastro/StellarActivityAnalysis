The archive you have received contains severals files and are required to run the analysis routine.
-epstopdf.pl
-Line_Setup.fits
-README
-RVSPY_analyze_activity.sav
-setup.ini

Requirements
------------

-IDL V8.2 or higher or IDL Virtual Machine
-The Fortran version of GLS has to be installed. The executable binary needs to be named 'GLS'.

Setup
-----

Following setup and directory structure is required:

1.) Create a folder for each star using the ID of the star.

2.) Inside this folder needs to be a folder named 'ceres'. The final reduced spectra of the ceres pipeline '*_sp.fits' and the correspondig CCF ascii files need to be there. The CCF file has to have the identical identifier as the _sp.fits. E.g. HD1835_20181205_UT04:09:47.928_CCF.txt and HD1835_20181205_UT04:09:47.928_sp.fits. If one or the other file is missing, the routine will crash.

3.) setup.ini:
    -This file has to be inside the same folder as RVSPY_analyze_activity.sav.
    -line 4: provide the path to the directory, which contains all the stars, e.g. /path/stars/
    -line 7: provide the path to a text fiel containing the target(s) to processed (one target per line, name has to be identical to the name of the target folder), e.g. /path/to/targets.txt
    -line 10: provide the path to GLS, e.g. /path/GLS_v2.3.02/
    -line 13: provide the path containing the script 'epstopdf.pl' (provided in this package), e.g. /path/ (do not add epstopdf.pl)
    -line 16: provide the path to the directory containing the file 'Line_Setup.fits', e.g. /path/Line_Setup.fits
    -line 19-23: Setup parameters for GLS, adjust as desired
    
4.) Line_Setup.fits contains the selected spectral lines, central wavelngth, and selected window size for the different values to be computed (e.g. equivalent width, ...).
    -This file has to be inside the same folder as RVSPY_analyze_activity.sav
    -The central wavelengths of the line are given in Angstrom
    -The wavelengths of the window size are given in Angstrom. The value provided is half the window size only. I.e. if you want a 15 Angstrom wide window a value of 7.5 has to be provided.
    -The (half) window size for the velocity correlation is in km/s.
    -The values for the window size have to be adjusted depending mainly on projected rotational velocity and activity level (e.g. broad Halpha emission vs narrow Halpha absorption).
    -Additional lines can be added if inside the spectral coverage of the data.
    

How to run the routine?
-----------------------

-Open a terminal and start IDL inside the directory containing RVSPY_analyze_activity.sav.
-type: restore, 'RVSPY_analyze_activity.sav'
-type: RVSPY_analyze_activity
-You will be asked if you want to prepare data. If the target is analysed for the first time or new data were added, answer with 'y'. It will create ascii versions of the spectra and a mean spectrum, which are needed for further processing. 

Results
-------

-The results are written in the 'Results' folder insside the folder of the processed target star.
-All values are stored in IDL SAVE files (can be restored in python) and therefore all plots can be reproduced if desired.


Good to know
------------

Depending on the number of epochs the analysis might take longer.
If stars have only one or two epochs the script will likely crash.
If the routine crashes, check for bad values in the CCF and spectra files and rerun ceres and/or remove the spectrum from folder.
