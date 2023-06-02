# Quick script to mask land values in 2 degree Munday config.
# Crudely done based on knowing where the land is. Ideally want to 
# do thisin MITGCM as part of the outputting

import netCDF4 as nc4
import numpy as np


#------------------------------------------------------------
# Add Land Mask based on knowing the locations where land is 
#------------------------------------------------------------
directory = '/data/hpcdata/users/racfur/MITgcm/verification/MundaySectorConfig_2degree/runs/100yrs/mnc_test_0002/'

# Open source file
ncfile = directory+'cat_tave.nc'
dset = nc4.Dataset(ncfile, 'r+')

# Get dimensions
Z_nc = dset.dimensions['Z'].name
Y_nc = dset.dimensions['Y'].name
X_nc = dset.dimensions['X'].name

# add the mask
Mask = dset.createVariable('Mask', 'i4', (Z_nc, Y_nc, X_nc))
dset.variables['Mask'][:,:,:]=1
dset.variables['Mask'][:,-2:,:]=0
dset.variables['Mask'][:,16:,-1:]=0
dset.variables['Mask'][32:,:,-1:]=0

# Save the file
dset.close()
