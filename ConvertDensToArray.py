# Script to read in MITGCM 'layers' output and save as npy array for quick and easy reading
# (note the mds routine is slow, so don't want to use this in my ReadRoutines each time I 
# change the dataset)
import sys
sys.path.append('.')
import mds as mds
import netCDF4 as nc4
import numpy as np

directory = '/data/hpcdata/users/racfur/MITgcm/verification/MundaySectorConfig_2degree/runs/100yrs/'
itrs = mds.scanforfiles(directory+'/layers_prho-tave')
density_full = mds.rdmds(directory+'layers_prho-tave',
                    itrs=itrs,machineformat='b',rec=None,fill_value=0,
                    returnmeta=False,astype=float,region=None,lev=(),
                    usememmap=False,mm=False,squeeze=True,verbose=False)
print(density_full.shape)  # (36000, 42, 78, 11)

#Discard first 50 years:
density = density_full[18000:,:,:,:]
print(density.shape)  # (36000, 42, 78, 11)

output_file = directory+'/DensityData.npy'
np.save(output_file, density)
