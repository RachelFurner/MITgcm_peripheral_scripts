# Quick script to rebuild all the MITgcm files from my run


for filename in 'daily_ave_2d' 'daily_ave_3d' 'monitor_grid' 'state' 'grid' 'phiHyd' 'monitor' 'phiHydLow' 
do
  ../../../../utils/python/MITgcmutils/scripts/gluemncbig -o $filename'.nc' 'mnc_output_'????/$filename'.??????????.t???.nc'
done

