undef("infill0")
procedure infill0(argoHeatnc:string, maskname:string, ncSaveName:string, figDir:string)
; add zeros to empty argo regions
 
begin

   f = addfile(argoHeatnc, "r")
   heatContentMap = f->heatContentMap
     
   d = dimsizes(heatContentMap)
   ntime = d(0)
   nlat = d(1)
   nlon = d(2)
   
   maskVals := asciiread(maskname, -1, "integer")
   maskVals := reshape(maskVals, (/nlon, nlat/))
   maskVals := transpose(maskVals)
   maskVals!0 = "lat"
   maskVals!1 = "lon"
   maskVals&lat = heatContentMap&LATITUDE
   maskVals&lon = heatContentMap&LONGITUDE
   
   ; replace argo values with 0s (i.e. no anomaly) in missing regions
   mask2 = conform_dims(dimsizes(heatContentMap), maskVals, (/1,2/))
   locs = ind(ndtooned(mask2.and.ismissing(heatContentMap)))
   
   dummy = ndtooned(heatContentMap)
   dummy(ndtooned(locs)) = 0.
   
   heatContentMapFilled = reshape(dummy, d)
   delete(dummy)
   
   heatContentTS = dim_sum_Wrap(dim_sum_Wrap(heatContentMapFilled)) ; sum over lat and lon
   heatContentTS!0 = "TIME"
   heatContentTS&TIME = heatContentMap&TIME
   seasonalCycle = plotSeasonalCycle(heatContentTS, figDir, "argo.seasonal.cycle.filled.0", "png", "Argo heat content (zero anom)")
   climatology = seasonalCycle@climo
   climatology!0 = "month"
   climatology&month = ispan(1,12,1)
   devFromClimo = ndtooned(seasonalCycle@deviation)
   devFromClimo!0 = "TIME"
   devFromClimo&TIME = heatContentTS&TIME

   ; save results into a netcdf
   system("/bin/rm -f " + ncSaveName)
   ncdf = addfile(ncSaveName, "c")

   fAtt = True
   fAtt@title = "Seasonal cycle in Argo heat content (zero anom)"
   fAtt@source_file = argoHeatnc
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )
   filedimdef(ncdf,"TIME",-1,True)

   ; ncdf->heatContentMap = heatContentMapFilled
   ncdf->heatContentTS = heatContentTS
   ncdf->climatology = climatology
   ncdf->devFromClimo = devFromClimo

   
end