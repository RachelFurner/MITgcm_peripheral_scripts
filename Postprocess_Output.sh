# Script to contain nco commands to postprocess MITgcm output into more useable format

# String together all the files along the time axis, taking only delected variables of interest
ncrcat -v ETAtave,Stave,T,Ttave,X,Xp1,Y,Yp1,Z,uVeltave,vVeltave,wVeltave tave.??????????.t???.nc cat_tave.nc

# String together gm output files to get bolus velocities
ncrcat gm_tave.0000000000.t001.nc cat_gm_tave.nc

# Discard first 50 years as spin up
ncks -d T,18000,35999 -o cat_tave.nc cat_tave.nc
ncks -d T,18000,35999 -o cat_gm_tave.nc cat_gm_tave.nc

# Append bolus velocities to the main file
ncks -A cat_gm_tave.nc cat_tave.nc



# Deal with monitor output files - again concatenate and discard first 50 years
ncrcat monitor.??????????.t???.nc cat_monitor.nc
ncks -d T,600,1200 -o cat_monitor.nc cat_monitor.nc


