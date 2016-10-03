undef ("cesmHeat")
procedure cesmHeat(figDir:string, ncSaveName:string, yearsToUse[*]:integer)
; Calculate seasonal cycle of heat content within CESM 
; uses simulation 002 from CESM1 Large Ensemble 
 
local rho0, cp

begin
  
   load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
   
   ; use same constants as argo analysis
   cp = 3990 ; J kg^-1 K^-1, water (Trenberth and Fasullo, 2008, JClim)
   a = 6.37122e6 ; m, radius of the earth
   pi = 4.*atan(1.0)
   secperday = 24*3600
   SA = 4*pi*a^2
   rho0 = 1026.5 ; kg m^-3 (Trenberth and Fasullo, 2008, JClim)

   ; pull down two cesm files POP files to combine
   ; these are available on yellowstone for those with access
   ; otherwise they can be accessed via ESGF: https://www.earthsystemgrid.org/dataset/ucar.cgd.ccsm4.CESM_CAM5_BGC_LE.html
   
   popDir = "/glade/p/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/TEMP"
   runNumber = "002"
   fnames = systemfunc("ls " + popDir + "/b.e11.*C5CNBDRD.f09_g16." + runNumber + ".pop.h.TEMP.*.nc")
   fnames := fnames(0:1) ; don't want future run

   f = addfiles(fnames, "r")
   time = f[:]->time
   TIME = cd_calendar(time, 0)
   ; pull out 2005-2015
   idxUse = ind((TIME(:, 0).ge.min(yearsToUse)).and.(TIME(:,0).le.max(yearsToUse)))
   T = f[:]->TEMP(idxUse, :, :, :)
   TLAT = f[0]->TLAT
   TLON = f[0]->TLONG
   TAREA = f[0]->TAREA ; cm^2
   dz = f[0]->dz ; (cm)
   z_t = f[0]->z_t ; depth from surface to midpoint of layer (cm)
   z_w_bot = f[0]->z_w_bot ; depth from surface to bottom of layer (cm)

   ; switch everything to meters
   cm2m = 100
   dz = dz/cm2m
   z_t = z_t/cm2m
   z_w_bot = z_w_bot/cm2m
   TAREA = TAREA/(cm2m^2)

   ; calculate heat budget down to 2000 meters
   idx2000 = closest_val(2000, z_w_bot)

   d = dimsizes(T)

   heatContent = cp*rho0*T ; Jm^-3

   heatContentAreaNorm = heatContent*conform_dims(dimsizes(heatContent), dz, 1) ; J m^-2
   
   ; add up heat vertically
   heatContentMap = dim_sum_n_Wrap(heatContentAreaNorm, 1)
   heatContentMap@units = "J m^-2"
   heatContentMap@lat2d = TLAT
   heatContentMap@lon2d = TLON
   heatContentMap!0 = "time"
   heatContentMap&time = T&time
   heatContentMap!1 = "nlat"
   heatContentMap!2 = "nlon"

   ; regrid
   ; Need weights to regrid from gx1v6 to 1x1 in a conservative manner
   ; Weights are created through ESMF_RegridWeightGen
   regridwgts = "/glade/u/home/mckinnon/gx1v6_1x1d/map_gx1v6_TO_1x1d_aave.160512.nc"
   heatContentMap1x1 = ESMF_regrid_with_weights(heatContentMap, regridwgts, False)
   
   ; check that total amount adds up across both estimates
   ; checks out!
   dA1x1 = area_global_rectilinear_grid(heatContentMap1x1&lat, heatContentMap1x1&lon, False)
   ; switch from km^2 to m^2
   dATemp = (10^3)^2 * dA1x1
   copy_VarMeta(dA1x1, dATemp)
   dATemp@units = "m^2"
   dA1x1 = dATemp
   delete(dATemp) 
   
   ; heatContentTotalOrig = heatContentMap*conform_dims(dimsizes(heatContentMap), TAREA, (/1,2/))
   heatContentTotal1x1 = heatContentMap1x1*conform_dims(dimsizes(heatContentMap1x1), dA1x1, (/1,2/))
   
   copy_VarMeta(heatContentMap1x1, heatContentTotal1x1)
   ; copy_VarMeta(heatContentMap, heatContentTotalOrig)
   
   delete( heatContentMap )
   delete( heatContentMap1x1 )
   
   ; add up heat from 2000 meters to surface
   heatContentMap = dim_sum_n_Wrap(heatContentAreaNorm(:, :idx2000, :, :), 1)
   heatContentMap@units = "J m^-2"
   heatContentMap@lat2d = TLAT
   heatContentMap@lon2d = TLON
   heatContentMap!0 = "time"
   heatContentMap&time = T&time
   heatContentMap!1 = "nlat"
   heatContentMap!2 = "nlon"

   ; regrid
   heatContentMap1x1 = ESMF_regrid_with_weights(heatContentMap, regridwgts, False)
   
   ; heatContent2000Orig = heatContentMap*conform_dims(dimsizes(heatContentMap), TAREA, (/1,2/))
   heatContent20001x1 = heatContentMap1x1*conform_dims(dimsizes(heatContentMap1x1), dA1x1, (/1,2/))
   copy_VarMeta(heatContentMap1x1, heatContent20001x1)
   ; copy_VarMeta(heatContentMap, heatContent2000Orig)
 
   ; heatContentTotalOrig@units = "J"
   heatContentTotal1x1@units = "J"
   heatContentTotal1x1@longname = "Heat content through full column"
   heatContent20001x1@units = "J"
   heatContent20001x1@longname = "Heat content down to 2000 meters"

   ; sumOrig = sum(ndtooned(heatContentTotalOrig))
   ; sum1x1 = sum(ndtooned(heatContentTotal1x1))
   ; percentDiff = (sum1x1 - sumOrig)/sumOrig
   ; print("Percent difference between maps is " + 100*percentDiff + "%")
   
   ; plot seasonal cycle
   heatContentTS = dim_sum_n_Wrap(dim_sum_n_Wrap(heatContent20001x1, 2), 1)
   seasonalCycle = plotSeasonalCycle(heatContentTS, figDir, "cesm.seasonal.cycle", "png", "CESM heat content")

   climatology = seasonalCycle@climo
   climatology!0 = "month"
   climatology&month = ispan(1,12,1)
   devFromClimo = ndtooned(seasonalCycle@deviation)
   devFromClimo!0 = "time"
   devFromClimo&time = heatContentTS&time

   ; save results into a netcdf
   system("/bin/rm -f " + ncSaveName)
   ncdf = addfile(ncSaveName, "c")

   fAtt = True
   fAtt@title = "Seasonal cycle in CESM heat content"
   fAtt@source_file = fnames
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )
   filedimdef(ncdf,"time",-1,True)

   ncdf->heatContentMapFull = heatContentTotal1x1
   ncdf->heatContentMap = heatContent20001x1
   ncdf->heatContentTS = heatContentTS
   ncdf->climatology = climatology
   ncdf->devFromClimo = devFromClimo

end