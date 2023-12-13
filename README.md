Supervised by Dr. Colin Goldblatt and Dr. Carsten Abraham


This repository contains the code used to analyze model output and make the figures for the undergraduate thesis titled:
'The Impact of the Meridional Temperature Gradient on the Spatial Distribution of Chemical Weathering Fluxes' by Keighan Gemmell.


Model output for the 3 CMIP6 CESM2 experiments used can be downloaded from: https://esgf-node.llnl.gov/search/cmip6/: 
You will need to select the CESM model; experiments: piControl, abrupt 2x CO2, abrupt 4x CO2; output variables: surface temperature, precipitation flux, and runoff flux. 


Model output for the S/S0 experiments can be downloaded from: https://doi.org/10.20383/101.0308.
	The output used for this analysis can be found under the 'CAM4' folder, with file names 'cam4_0x.nc' 
	where x is 8,9,1,11 representing solar constants used. 
	
	
Model output for the Eocene experiments can be downloaded from: https://www.deepmip.org/data-eocene/, following the steps outlined on that page and using the CESM model. 


As the data is found in NetCDF files, you will need to download and install the Climate Data Operator (CDO)  (https://code.mpimet.mpg.de/projects/cdo/wiki#Installation-and-Supported-Platforms), 
NetCDF (https://www.unidata.ucar.edu/downloads/netcdf/) in order to read in the files. 


Lithology data can be downloaded from: 10.1594/PANGAEA.788537. This data requires processing which can be found in CMIP6 folder weathering map scripts.

Soil shielding data can be downloaded from: https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1247. This data requires processing, 
which can be found in the soil shielding python script found in the CMIP6 folder. 


Scripts for processing and analyzing model output is written in python and make use of the climate data operators outlined above. 

Scripts are divided by model output source, and each contain scripts for global weathering fluxes, zonal weathering flux contributions, and cumulative fluxes.
The different model output source folders contain a few different additional scripts needed for a full analysis of the data. 
Each file name indicates what the script produces, and the directory in which it is located will indicate for which model output. 

Residence time scripts use the pyCO2sys package which can be downloaded from: https://pyco2sys.readthedocs.io/en/latest/. 


Many of the scripts use previously calculated variables (stored as .npy files), which are located in the directories in which they are used. 
The best way to ensure there are no errors with this is to download this entire repository, including .npy files. 

Paths to data files in all scripts are the paths used to where the data is stored in my directory on the anahim cluster.
If you download data to another location, or if you are not on the anahim cluster, the paths to data files will need to be changed in all scripts. 

Questions about the code in this repository can be sent to Keighan Gemmell (keighan.gemmell@gmail.com). 

