undef("infillGlobalAvg")
procedure infillGlobalAvg(ncDir:string, argoHeatnc:string, maskname:string, ncSaveName:string, figDir:string)
; add the global average to missing regions (volume normalized)

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
   
   ; replace argo values with global average of available values in missing regions
   ; should be normalized by volume
    
   mask2 = conform_dims(dimsizes(heatContentMap), maskVals, (/1,2/))
   noargo = mask2.and.ismissing(heatContentMap)
   locs = ind(ndtooned(noargo)) ; locations in the ocean without argo data
   
   globalAvgHeat = dim_sum_Wrap(dim_sum_Wrap(heatContentMap)) ; measured by argo in upper 2000 meters
   dA = area_global_rectilinear_grid(heatContentMap&LATITUDE, heatContentMap&LONGITUDE, False)
   ; switch from km^2 to m^2
   dATemp = (10^3)^2 * dA
   copy_VarMeta(dA, dATemp)
   dATemp@units = "m^2"
   dA = dATemp
   delete(dATemp)
   
   ; get ocean depth
   ncVol = "./data/ocean_depth.nc"
   fVol = addfile(ncVol , "r")
   zi = fVol->depth
   
   ; check lats/lons
   if (any(dA&lat.ne.zi&lat)).or.(any(dA&lon.ne.zi&lon))
      print("lats/lons not aligned! Abort!")
   end if
   
   ; make second zi capped at 2000 meters -- argo domain
   zi2000 = where(zi.gt.2000, 2000, zi)
   
   ; create volume 
   dV = zi2000*dA
   copy_VarMeta(dA, dV)
   dV@units = "m^3"
   dV@long_name = "Volume of each gridcell of ocean; depth capped at 2000 m"
   
   ; also for full ocean
   dVFull = zi*dA
   copy_VarMeta(dA, dVFull)
   dVFull@units = "m^3"
   dVFull@long_name = "Volume of each gridcell of ocean"
   
   ; calculate fraction of volume with argo data
   volVec = ndtooned(dV)
   locsWithData := .not.ismissing( heatContentMap(0 , :, :))
   ; total volume of locations that have data
   dVtotal = dim_sum(volVec(ind(ndtooned(locsWithData))))
   dVtotal@units = "m^3"
   
   ; calculate fraction of volume without argo data
   locsNoData = noargo(0, :, :)
   dVtotalMissing = dim_sum(volVec(ind(ndtooned(locsNoData))))
   dVtotalMissing@units = "m^3"
   
   ; print fraction of upper 2000 meters measured by argo
   fracMeasured2000 = (dVtotal)/(dVtotal+dVtotalMissing)
   print("Fraction of upper 2000 meters measured by argo: " + fracMeasured2000)
   
   ; assign this value to missing gridboxes x dV
   ; assumption: measured part of ocean is representative of unmeasured    
   globalAvgHeatNorm = globalAvgHeat/dVtotal
   globalAvgHeatInfill = conform_dims(dimsizes(heatContentMap), globalAvgHeatNorm, 0) * conform_dims(dimsizes(heatContentMap), dV, (/1, 2/))
   globalAvgHeatInfillVec = ndtooned(globalAvgHeatInfill)
    
   dummy = ndtooned(heatContentMap)
   dummy(ndtooned(locs)) = globalAvgHeatInfillVec(ndtooned(locs)) ; infill regions without argo
   
   heatContentMapFilled = reshape(dummy, d)
   delete(dummy)
    
   wks = gsn_open_wks("png", figDir + "/argodatafilledGlobal")
   res = True
   res@cnFillOn = True
   res@cnFillMode = "RasterFill"
   res@cnLinesOn = False
   plot = gsn_csm_contour(wks, heatContentMapFilled(0,:,:), res)
   
   heatContentTS = dim_sum_Wrap(dim_sum_Wrap(heatContentMapFilled)) ; sum over lat and lon
   heatContentTS!0 = "TIME"
   heatContentTS&TIME = heatContentMap&TIME
   seasonalCycle = plotSeasonalCycle(heatContentTS, figDir, "argo.seasonal.cycle.filled.global", "png", "Argo heat content (global anom)")
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
   fAtt@title = "Seasonal cycle in Argo heat content (global anom)"
   fAtt@source_file = argoHeatnc
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )
   filedimdef(ncdf,"TIME",-1,True)

   ; ncdf->heatContentMap = heatContentMapFilled
   ncdf->heatContentTS = heatContentTS
   ncdf->climatology = climatology
   ncdf->devFromClimo = devFromClimo
   
   ; also save dV
   dVname = ncDir + "/dV2000.nc"
   system("/bin/rm -f " + dVname)
   ncdf = addfile(dVname, "c")

   fAtt = True
   fAtt@title = "Volume of ocean down to 2000 meters"
   fAtt@source_file = ncVol
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )

   ncdf->dV = dV
   
   ; and dVFull
   dVFullname = ncDir + "/dV.nc"
   system("/bin/rm -f " + dVFullname)
   ncdf = addfile(dVFullname, "c")

   fAtt = True
   fAtt@title = "Volume of ocean"
   fAtt@source_file = ncVol
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )
   ncdf->dV = dVFull


   
end