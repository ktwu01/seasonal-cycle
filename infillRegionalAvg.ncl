undef("infillRegionalAvg")
procedure infillRegionalAvg(argoHeatnc:string, maskname:string, ncSaveName:string, figDir:string)
; add the regional average to missing regions

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

   ; replace argo values with regional averages (tropics, NH, SH)

   mask2 = conform_dims(dimsizes(heatContentMap), maskVals, (/1,2/))
   mask2!0 = "time"
   mask2!1 = "lat"
   mask2!2 = "lon"
   mask2&time = heatContentMap&TIME
   mask2&lat = maskVals&lat
   mask2&lon = maskVals&lon

   dA = area_global_rectilinear_grid(heatContentMap&LATITUDE, heatContentMap&LONGITUDE, False)
   ; switch from km^2 to m^2
   dATemp = (10^3)^2 * dA
   copy_VarMeta(dA, dATemp)
   dATemp@units = "m^2"
   dA = dATemp
   delete(dATemp)

   ; load in volume info
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

   lat1 = (/-90, -23, 23/)
   lat2 = (/-23, 23, 90/)

   dummy = ndtooned(heatContentMap)
   regionalAvgHeatInfill = new(dimsizes(heatContentMap), typeof(heatContentMap))

   do ct = 0, (dimsizes(lat1)-1)
      latIdx := ind((heatContentMap&LATITUDE.gt.lat1(ct)).and.(heatContentMap&LATITUDE.lt.lat2(ct)))
      regionalHeat = dim_sum_Wrap(dim_sum_Wrap(heatContentMap(:, latIdx, :))) ; dimensions of time
      volVec := ndtooned(dV(latIdx, :))
      ; get locations that have measurements
      locsWithData := .not.ismissing( heatContentMap(0 , latIdx, :))
      ; total volume of locations that have data
      dVtotal = dim_sum(volVec(ind(ndtooned(locsWithData))))
      dVtotal@units = "m^3"


      regionalAvgHeatNorm := regionalHeat/dVtotal ; time series, heat normalized by volume  
 
      ; get heat values for regions that need to be infilled
      ; units of J
      regionalAvgHeatInfill(:, latIdx, :)= conform_dims(dimsizes(heatContentMap(:, latIdx, :)), regionalAvgHeatNorm, 0) * conform_dims(dimsizes(heatContentMap(:, latIdx, :)), dV(latIdx, :), (/1, 2/))

   end do

   locs = ind(ndtooned(mask2.and.ismissing(heatContentMap)))
   regionalAvgHeatInfillVec = ndtooned(regionalAvgHeatInfill)

   dummy(locs) = regionalAvgHeatInfillVec(locs)
   heatContentMapFilled = reshape(dummy, d)
   copy_VarCoords(heatContentMap, heatContentMapFilled)

   wks = gsn_open_wks("png", figDir + "/argodatafilledRegional")
   res = True
   res@cnFillOn = True
   res@cnFillMode = "RasterFill"
   res@cnLinesOn = False
   plot = gsn_csm_contour(wks, heatContentMapFilled(0,:,:), res)

   heatContentTS = dim_sum_Wrap(dim_sum_Wrap(heatContentMapFilled)) ; sum over lat and lon
   heatContentTS!0 = "TIME"
   heatContentTS&TIME = heatContentMap&TIME
   seasonalCycle = plotSeasonalCycle(heatContentTS, figDir, "argo.seasonal.cycle.filled.regional", "png", "Argo heat content (regional anom)")
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
   fAtt@title = "Seasonal cycle in Argo heat content (regional anom)"
   fAtt@source_file = argoHeatnc
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )
   filedimdef(ncdf,"TIME",-1,True)

   ncdf->heatContentMap = heatContentMapFilled
   ncdf->heatContentTS = heatContentTS
   ncdf->climatology = climatology
   ncdf->devFromClimo = devFromClimo


end