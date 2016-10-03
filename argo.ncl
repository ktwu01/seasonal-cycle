undef ("argo")
procedure argo(figDir:string, ncSaveName:string, yearsToUse[*]:integer)
; Calculate ocean heat content based on Scripps Argo product

begin 

   ; Get Argo heat budget
 
   cp = 3990 ; J kg^-1 K^-1, water (Trenberth and Fasullo, 2008, JClim)
   a = 6.37122e6 ; m, radius of the earth
   pi = 4.*atan(1.0)
   secperday = 24*3600
   SA = 4*pi*a^2
   rho0 = 1026.5 ; kg m^-3 (Trenberth and Fasullo, 2008, JClim)

   ; (1) argo calculations
   argoFname = "/glade/scratch/mckinnon/argo/RG_ArgoClim_Temp.nc"
   f = addfile(argoFname, "r")

   ; chop to span 'yearsToUse'
   TIME = f->TIME
   TIME := cd_calendar(TIME, 0)
   idxUse = ind((TIME(:, 0).ge.min(yearsToUse)).and.(TIME(:,0).le.max(yearsToUse)))

   T = f->ARGO_TEMPERATURE_ANOMALY(idxUse, :, :, :)
   T0 = f->ARGO_TEMPERATURE_MEAN

   delete(TIME)
   delete(idxUse)

   d = dimsizes(T)

   ; alter Argo longitude
   T = lonPivot(T, 360.5)
 
   ; expand grid to spread -90 to 90 in longitude
   latExpand = latGlobeFo(180,"lat","latitude","degrees_north")
   Tnew = new((/d(0), d(1), 180, d(3)/), typeof(T), T@_FillValue)
   Tnew!2 = "LATITUDE"
   Tnew&LATITUDE = latExpand
   Tnew(:, :, {-64.5:64.5}, :) = T
 
   T := Tnew
   delete(Tnew)

   d = dimsizes(T)
   npres = d(1)
   nlat = d(2)

   ; get depths
   ; expand pressure and latitude
   pressureMat = conform_dims((/d(1), d(2)/), T&PRESSURE, 0)
   latMat = conform_dims((/d(1), d(2)/), T&LATITUDE, 1)
   z = pres2depth(pressureMat, latMat)

   midz = 0.5*(z(0:npres-2, :) + z(1:npres-1, : ))
   top = conform_dims((/1,nlat/),0.,0)
   bottom = conform_dims((/1,nlat/),2000.,0)
   total = array_append_record(array_append_record(top, midz, 0), bottom, 0)
 
   dz = total(1:npres, :) - total(0:npres-1,:)
   dz!0 = "PRESSURE"
   dz!1 = "LATITUDE"
   dz&PRESSURE = T&PRESSURE
   dz&LATITUDE = T&LATITUDE
   dz@units = "meters"
   dz@long_name = "depth associated with argo layer"

   if (any(dz.lt.0))
      print("Dz values less than zero!! Stop!")
   end if

   heatContent = cp*rho0*T
   copy_VarCoords(T, heatContent)
   heatContent@units = "J m^{-3}"
   heatContent@long_name = "Heat content based on 1x1 Argo measurements"
   heatContent@temp_data_source = argoFname

   printVarSummary(heatContent)
  
   ; multiply heat content by volume
   dA = area_global_rectilinear_grid(T&LATITUDE, T&LONGITUDE, False)
   ; switch from km^2 to m^2
   dATemp = (10^3)^2 * dA
   copy_VarMeta(dA, dATemp)
   dATemp@units = "m^2"
   dA = dATemp
   delete(dATemp)
   
   ; save dA for later
   system("/bin/rm -f /glade/u/home/mckinnon/seasonalCycleHeatContent/output/dAGlobal.nc")
   ncdf = addfile("/glade/u/home/mckinnon/seasonalCycleHeatContent/output/dAGlobal.nc", "c")
   fAtt = True
   fAtt@title = "dA for 1x1 grid"
   fAtt@source_fn = "area_global_rectilinear_grid"
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )
   ncdf->dA = dA 
 
   dV = conform_dims(d(1:3), dz, (/0, 1/)) * conform_dims(d(1:3), dA, (/1, 2/))

   heatContentTotal = heatContent*conform_dims(dimsizes(heatContent), dV, (/1, 2, 3/))
   copy_VarMeta(heatContent, heatContentTotal)
   heatContentTotal@units = "J"
   printVarSummary(heatContentTotal)

   ; add heat content vertically
   ; TESTING! ADD ONLY TO 1000 METERS
   ;idx = ind(T&PRESSURE.eq.1000)
   ; print(idx)
   ; heatContentMap =  dim_sum_n_Wrap(heatContentTotal(:, 0:idx, :, :), 1)
   heatContentMap =  dim_sum_n_Wrap(heatContentTotal, 1)
   
   printVarSummary(heatContentMap)

   ; add up all heat content and plot seasonal cycle
   heatContentTS = dim_sum_Wrap(dim_sum_Wrap(heatContentMap)) ; sum over lat and lon
   seasonalCycle = plotSeasonalCycle(heatContentTS, figDir, "argo.seasonal.cycle", "png", "Argo heat content")
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
   fAtt@title = "Seasonal cycle in Argo heat content"
   fAtt@source_file = argoFname
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )
   filedimdef(ncdf,"TIME",-1,True)

   ncdf->heatContentMap = heatContentMap
   ncdf->heatContentTS = heatContentTS
   ncdf->climatology = climatology
   ncdf->devFromClimo = devFromClimo


end