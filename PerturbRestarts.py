import netCDF4 as nc4
import numpy as np
import shutil
import glob
import gcm_filters as gcm_filters
import xarray as xr

random_noise=False
NN_restart=False
smoothing= True

data_dir = '/data/hpcdata/users/racfur/MITgcm/verification/MundayChannelConfig10km_LandSpits/runs/'

if smoothing:
   # Take standard restart and smooth with gcm filters 
   orig_file = data_dir+'/50yr_SpinUp/pickup.0002592000.nc'
   grid_file = data_dir+'/50yr_SpinUp/grid.nc'
   grid_ds = xr.open_dataset(grid_file)
   HFacC = grid_ds['HFacC'].values
   Tmask = np.ones(HFacC.shape)
   Tmask[:,:,:] = np.where( Tmask[:,:,:] > 0., 1, 0 )
   Tmask = xr.DataArray(Tmask, dims=['Z', 'Y', 'X'])
   Etamask = np.ones((HFacC.shape[1], HFacC.shape[2]))
   Etamask[:,:] = np.where( Etamask[:,:] > 0., 1, 0 )
   Etamask = xr.DataArray(Etamask, dims=['Y', 'X'])
   HFacW = grid_ds['HFacW'].values
   Umask = np.ones(HFacW.shape)
   Umask[:,:,:] = np.where( Umask[:,:,:] > 0., 1, 0 )
   Umask = xr.DataArray(Umask, dims=['Z', 'Y', 'Xp1'])
   HFacS = grid_ds['HFacS'].values
   Vmask = np.ones(HFacS.shape)
   Vmask[:,:,:] = np.where( Vmask[:,:,:] > 0., 1, 0 )
   Vmask = xr.DataArray(Vmask, dims=['Z', 'Yp1', 'X'])
   for smoothing_scale in [125, 150]:
      xr_dset = xr.open_dataset(orig_file)
      filter = gcm_filters.Filter( filter_scale=smoothing_scale, dx_min=10, filter_shape=gcm_filters.FilterShape.TAPER,
                                   grid_type=gcm_filters.GridType.REGULAR_WITH_LAND,
                                   grid_vars={'wet_mask': Tmask} )
      new_Temp = filter.apply( xr_dset['Temp'], dims=['Y', 'X']).values
      filter = gcm_filters.Filter( filter_scale=smoothing_scale, dx_min=10, filter_shape=gcm_filters.FilterShape.TAPER,
                                   grid_type=gcm_filters.GridType.REGULAR_WITH_LAND,
                                   grid_vars={'wet_mask': Etamask} )
      new_Eta  = filter.apply( xr_dset['Eta'], dims=['Y', 'X']).values
      filter = gcm_filters.Filter( filter_scale=smoothing_scale, dx_min=10, filter_shape=gcm_filters.FilterShape.TAPER,
                                   grid_type=gcm_filters.GridType.REGULAR_WITH_LAND,
                                   grid_vars={'wet_mask': Umask} )
      new_U    = filter.apply( xr_dset['U'], dims=['Y', 'Xp1']).values
      filter = gcm_filters.Filter( filter_scale=smoothing_scale, dx_min=10, filter_shape=gcm_filters.FilterShape.TAPER,
                                   grid_type=gcm_filters.GridType.REGULAR_WITH_LAND,
                                   grid_vars={'wet_mask': Vmask} )
      new_V    = filter.apply( xr_dset['V'], dims=['Yp1', 'X']).values
      proc=1
      for yproc in range(2):
         for xproc in range(6):
            print(xproc, yproc)
            perturbed_file = data_dir+'/50yr_smooth'+str(smoothing_scale)+'/pickup.0002592000.t'+str(proc).zfill(3)+'.nc'
            orig_proc_file = glob.glob(data_dir+'/50yr_SpinUp/mnc_output_00??/pickup.0002592000.t'+str(proc).zfill(3)+'.nc')
            print(orig_proc_file[0])
            shutil.copyfile(orig_proc_file[0], perturbed_file)
            nc_dset = nc4.Dataset(perturbed_file,'r+')
            nc_dset['Temp'][:,:,:,:] = new_Temp[:,:,yproc*52:(yproc+1)*52,xproc*40:(xproc+1)*40]
            nc_dset['U'][:,:,:,:] = new_U[:,:,yproc*52:(yproc+1)*52,xproc*40:(xproc+1)*40+1]
            nc_dset['V'][:,:,:,:] = new_V[:,:,yproc*52:(yproc+1)*52+1,xproc*40:(xproc+1)*40]
            nc_dset['Eta'][:,:,:] = new_Eta[:,yproc*52:(yproc+1)*52,xproc*40:(xproc+1)*40]
            nc_dset.close()
            proc=proc+1

if random_noise:
# Take standard restart and perturb by random noise
   for perturbation in range(5):
      for orig_file in glob.glob(data_dir+'/50yr_SpinUp/mnc_output_000?/pickup.0002592000.t???.nc'):
         perturbed_file = orig_file[:-3]+'_pert'+str(perturbation)+'.nc'
         shutil.copyfile(orig_file, perturbed_file)
         dset = nc4.Dataset(perturbed_file,'r+')
         dset['Temp'][:,:,:,:] = dset['Temp'][:,:,:,:]+np.random.random(dset['Temp'].shape) * 2 * (10**-10) - 10**-10
         dset['U'][:,:,:,:] = dset['U'][:,:,:,:]+np.random.random(dset['U'].shape) * 2 * (10**-10) - 10**-10
         dset['V'][:,:,:,:] = dset['V'][:,:,:,:]+np.random.random(dset['V'].shape) * 2 * (10**-10) - 10**-10
         dset['Eta'][:,:,:] = dset['Eta'][:,:,:]+np.random.random(dset['Eta'].shape) * 2 * (10**-10) - 10**-10
         dset.close()
      for orig_file in glob.glob(data_dir+'/50yr_SpinUp/mnc_output_001[012]/pickup.0002592000.t???.nc'):
         perturbed_file = orig_file[:-3]+'_pert'+str(perturbation)+'.nc'
         shutil.copyfile(orig_file, perturbed_file)
         dset = nc4.Dataset(perturbed_file,'r+')
         dset['Temp'][:,:,:,:] = dset['Temp'][:,:,:,:]+np.random.random(dset['Temp'].shape) * 2 * (10**-10) - 10**-10
         dset['U'][:,:,:,:] = dset['U'][:,:,:,:]+np.random.random(dset['U'].shape) * 2 * (10**-10) - 10**-10
         dset['V'][:,:,:,:] = dset['V'][:,:,:,:]+np.random.random(dset['V'].shape) * 2 * (10**-10) - 10**-10
         dset['Eta'][:,:,:] = dset['Eta'][:,:,:]+np.random.random(dset['Eta'].shape) * 2 * (10**-10) - 10**-10
         dset.close()

if NN_restart:
   # Take prediction from one iteration of neural network, and create restart from this
   NN_restart_dir = '/data/hpcdata/users/racfur/MITgcm/verification/MundayChannelConfig10km_LandSpits/runs/RestartFromNNPred'
   
   NN_prediction_file = NN_restart_dir+'/Spits12hrly_UNet2dtransp_histlen1_seed30475_200epochs_simple_Forecast120.nc'
   NN_dset = nc4.Dataset(NN_prediction_file,'r')
   
   NN_p1U = np.zeros((1,38,104,241))
   NN_p1U[0,:,:,:-1] = NN_dset['Pred_U'][0:1,:,:,:]
   NN_p1U[0,:,:,-1] = NN_dset['Pred_U'][0:1,:,:,0]
   
   NN_p1U_Inc = np.zeros((1,38,104,241))
   NN_p1U_Inc[0,:,:,:-1] = NN_dset['Pr_U_Inc'][0:1,:,:,:]
   NN_p1U_Inc[0,:,:,-1] = NN_dset['Pr_U_Inc'][0:1,:,:,0]
   
   NN_p1V = np.zeros((1,38,105,240))
   NN_p1V[0,:,:-1,:] = NN_dset['Pred_V'][0:1,:,:,:]
   NN_p1V[0,:,-1,:] = NN_dset['Pred_V'][0:1,:,0,:]
   
   NN_p1V_Inc = np.zeros((1,38,105,240))
   NN_p1V_Inc[0,:,:-1,:] = NN_dset['Pr_V_Inc'][0:1,:,:,:]
   NN_p1V_Inc[0,:,-1,:] = NN_dset['Pr_V_Inc'][0:1,:,0,:]
   
   S_array = np.ones((1,38,52,40))
   S_array[:,:,:,:] = 35. 
   
   proc=1
   for yproc in range(2):
      for xproc in range(6):
         print(xproc, yproc)
         NN_perturbed_file = NN_restart_dir+'/pickup.0002592001.t'+str(proc).zfill(3)+'.nc'
         print(NN_perturbed_file)
         nc_file = nc4.Dataset(NN_perturbed_file, 'w', format='NETCDF3_CLASSIC')
         nc_file.createDimension('T', None)
         nc_file.createDimension('Z', 38)
         nc_file.createDimension('Y', 52)
         nc_file.createDimension('X', 40)
         nc_file.createDimension('Yp1', 53)
         nc_file.createDimension('Xp1', 41)
         nc_time = nc_file.createVariable( 'T', 'i', ('T') )
         nc_Temp = nc_file.createVariable( 'Temp'  , 'f4', ('T', 'Z', 'Y', 'X') )
         nc_U    = nc_file.createVariable( 'U'     , 'f4', ('T', 'Z', 'Y', 'Xp1') )
         nc_V    = nc_file.createVariable( 'V'     , 'f4', ('T', 'Z', 'Yp1', 'X') )
         nc_Eta  = nc_file.createVariable( 'Eta'   , 'f4', ('T', 'Y', 'X') )
         nc_S    = nc_file.createVariable( 'S'     , 'f4', ('T', 'Z', 'Y', 'X') )
         nc_gTnm1= nc_file.createVariable( 'gTnm1' , 'f4', ('T', 'Z', 'Y', 'X') )
         nc_gUnm1= nc_file.createVariable( 'gUnm1' , 'f4', ('T', 'Z', 'Y', 'Xp1') )
         nc_gVnm1= nc_file.createVariable( 'gVnm1' , 'f4', ('T', 'Z', 'Yp1', 'X') )
         nc_dEtaHdt = nc_file.createVariable('dEtaHdt', 'f4', ('T', 'Y', 'X') )
         nc_gSnm1= nc_file.createVariable( 'gSnm1' , 'f4', ('T', 'Z', 'Y', 'X') )
         nc_time[:] = np.arange(1)
         nc_Temp[:,:,:,:] = NN_dset['Pred_Temp'][0:1,:,yproc*52:(yproc+1)*52,xproc*40:(xproc+1)*40]
         nc_U[:,:,:,:] = NN_p1U[:,:,yproc*52:(yproc+1)*52,xproc*40:(xproc+1)*40+1]
         nc_V[:,:,:,:] = NN_p1V[:,:,yproc*52:(yproc+1)*52+1,xproc*40:(xproc+1)*40]
         nc_Eta[:,:,:] = NN_dset['Pred_Eta'][0:1,yproc*52:(yproc+1)*52,xproc*40:(xproc+1)*40]
         nc_S[:,:,:] = S_array
         nc_gTnm1[:,:,:,:] = NN_dset['Pr_T_Inc'][0:1,:,yproc*52:(yproc+1)*52,xproc*40:(xproc+1)*40]
         nc_gUnm1[:,:,:,:] = NN_p1U_Inc[:,:,yproc*52:(yproc+1)*52,xproc*40:(xproc+1)*40+1]
         nc_gVnm1[:,:,:,:] = NN_p1V_Inc[:,:,yproc*52:(yproc+1)*52+1,xproc*40:(xproc+1)*40]
         nc_dEtaHdt[:,:,:] = NN_dset['Pr_E_Inc'][0:1,yproc*52:(yproc+1)*52,xproc*40:(xproc+1)*40]
         nc_gSnm1[:,:,:,:] = np.zeros((1,38,52,40))
         nc_file.tile_number = proc
         nc_file.bi = 1
         nc_file.bj = 1
         nc_file.sNx = 40
         nc_file.sNy = 52
         nc_file.Nx = 240
         nc_file.Ny = 104
         nc_file.Nr = 38
         nc_file.nSx = 1
         nc_file.nSy = 1
         nc_file.nPx = 6
         nc_file.nPy = 2
         nc_file.close()
         proc = proc+1
   NN_dset.close()
   
