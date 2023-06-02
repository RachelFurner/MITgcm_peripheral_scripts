# Script to contain nco commands to postprocess MITgcm output into more useable format

# String together the files along the time axis, taking only selected variables of interest
echo 'concatenating relevant variables from output files'
ncrcat -v ETAtave,Stave,T,Ttave,X,Xp1,Y,Yp1,Z,uVeltave,vVeltave,wVeltave tave.??????????.t???.nc cat_tave.nc
echo 'concatenating gm output files'
ncrcat gm_tave.??????????.t???.nc cat_gm_tave.nc
echo 'concatenating diag output files'
ncrcat diag.??????????.t???.nc cat_diag_tave.nc

# Discard first 50 years as spin up
echo 'ncks to remove first 50 years'
ncks -d T,18000,35999 -O cat_tave.nc cat_tave.nc
ncks -d T,18000,35999 -O cat_gm_tave.nc cat_gm_tave.nc
ncks -d T,18000,35999 -O cat_diag_tave.nc cat_diag_tave.nc

# Append bolus velocities and diag variables to the main file
echo 'adding bolus velocities to main file'
ncks -A cat_gm_tave.nc cat_tave.nc
echo 'adding diag variables to main file'
ncks -A cat_diag_tave.nc cat_tave.nc

# Deal with monitor output files - again concatenate and discard first 50 years
echo 'concatenating monitor files'
ncrcat monitor.??????????.t???.nc cat_monitor.nc
echo 'removing first 50 years from monitor file'
ncks -d T,600,1200 -O cat_monitor.nc cat_monitor.nc


