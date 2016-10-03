undef("seaiceHeat")
procedure seaiceHeat(figDir:string, ncSaveName:string, yearsToUse[*]:integer)
; Calculate seasonal cycle of heat content associated with sea ice changes

begin
   
   dataDir = "LOCATION WHERE YOU HAVE SAVED SEA ICE DATA"
   ; http://psc.apl.uw.edu/wordpress/wp-content/uploads/schweiger/ice_volume/PIOMAS.2sst.monthly.Current.v2.1.txt
   arcticFname = dataDir + "/PIOMAS.2sst.monthly.Current.v2.1.txt"
   
   ; https://www.ncdc.noaa.gov/snow-and-ice/extent/sea-ice/S/0.csv
   antarcticFname = dataDir + "/sh_seaiceextent.csv"
   
   ; constants
   rho = 917 ; kg m^-3
   Hfus = 334e3 ; J kg^-1, energy required to melt ice
   
   ; read in arctic sea ice
   nrow = numAsciiRow(arcticFname)
   ncol = numAsciiCol(arcticFname)
   ice1 = asciiread(arcticFname, (/nrow, ncol/), "float")
   
   data = ice1(:, 1:ncol-1) ; data is in 1000 km^3
   arcticSeaice = 1e3*(1e3)^3*ndtooned(data) ; now in m^3
   arcticSeaice@_FillValue = min(arcticSeaice) ; missing value is -1x10^12
   yrs = ice1(:, 0)
   
   ; make a time vector
   yyyymm = yyyymm_time(tointeger(min(yrs)), tointeger(max(yrs)), "integer")
   yyyy = yyyymm/100
   mm = yyyymm%100 ; this uses the modulus operator, month is now 1201 (http://www.ncl.ucar.edu/Document/Functions/Contributed/mod.shtml) 
   dd = mm
   hh = mm
   mn = mm
   sc = mm
   dd = 1 ; assign all the 1st of the month
   hh = 0
   mn = 0
   sc = 0
   time = cd_inv_calendar(yyyy, mm, dd, hh, mn, sc, "days since 1979-01-01 00:00:00", 0) 
   delete(yyyymm)
   delete(yyyy)
   delete(mm)
   delete(dd)
   delete(hh) 
   delete(mn)
   delete(sc)
 
   arcticSeaice!0 = "time"
   arcticSeaice&time = time
   arcticSeaice@units = "m^3"
   nTime = dimsizes(time)

   dIce = arcticSeaice(2:nTime-1) - arcticSeaice(0:nTime-3) 
   ; sign: more ice = more heat storage in system
   dTime = time(1:nTime-2)
   arcticHeat = dIce*rho*Hfus 
   arcticHeat!0 = "time"
   arcticHeat&time = dTime
   arcticHeat@units = "J"
   arcticHeat@source_file = arcticFname
   
   ; read in antarctic sea ice
   ice2 = asciiread(antarcticFname, (/445, 3/), "float") ; numAscii* functions not working    
   yyyymm = ice2(:, 0)
   yyyymm2 = floattointeger(yyyymm)
   yyyy = yyyymm2/100
   mm = yyyymm2%100 ; this uses the modulus operator, month is now 1201 (http://www.ncl.ucar.edu/Document/Functions/Contributed/mod.shtml) 
   dd = mm
   hh = mm 
   mn = mm
   sc = mm
   dd = 1 ; assign all the first of the month
   hh = 0
   mn = 0
   sc = 0
   delete(sc@_FillValue)
   timeAnt = cd_inv_calendar(yyyy, mm, dd, hh, mn, sc, "days since 1979-01-01 00:00:00", 0) 
 
   extent = ice2(:, 1) ; million km^2
   
   ; need to make an assumption about thickness, because not available
   antarctic_thickness = 0.9 ; m, http://onlinelibrary.wiley.com/doi/10.1029/2007JC004254/epdf
   antarcticSeaice = 1e6*(1e3)^2*extent*antarctic_thickness ; m^3
   
   delete(yyyymm)
   delete(yyyy)
   delete(mm)
   delete(dd)
   delete(hh)
   delete(mn)
   delete(sc)
   nTimeAnt = dimsizes(timeAnt)
   
   dIceAnt = antarcticSeaice(2:nTimeAnt-1) - antarcticSeaice(0:nTimeAnt-3) 
   
   dTimeAnt = timeAnt(1:nTimeAnt-2)
   antarcticHeat = dIceAnt*rho*Hfus 
   antarcticHeat!0 = "time"
   antarcticHeat&time = dTimeAnt
   antarcticHeat@units = "J"
   antarcticHeat@source_file = antarcticFname
    
   printVarSummary(antarcticHeat)
   printVarSummary(arcticHeat)
   
   ; really want to consider these together
   arcticTIME = cd_calendar(arcticHeat&time, 0)
   antarcticTIME = cd_calendar(antarcticHeat&time, 0)
   antarcticHeat := antarcticHeat(ind( (antarcticTIME(:, 0).ge.min(yearsToUse)).and.(antarcticTIME(:, 0).le.max(yearsToUse)) ))
   arcticHeat := arcticHeat(ind( (arcticTIME(:, 0).ge.min(yearsToUse)).and.(arcticTIME(:, 0).le.max(yearsToUse)) ))
   totalHeat = antarcticHeat + arcticHeat
       
   totalHeat!0 = "time"
   totalHeat&time = dTimeAnt(ind( (arcticTIME(:, 0).ge.min(yearsToUse)).and.(arcticTIME(:, 0).le.max(yearsToUse)) ))
   totalHeat@units = "J"
   
   ; plot seasonal cycle
   seasonalCycle = plotSeasonalCycle(antarcticHeat, figDir, "seaice.antarctic.seasonal.cycle", "png", "Antarctic Sea ice heat content")
   delete(seasonalCycle)
   seasonalCycle = plotSeasonalCycle(arcticHeat, figDir, "seaice.arctic.seasonal.cycle", "png", "Arctic Sea ice heat content")
   delete(seasonalCycle)
   seasonalCycle = plotSeasonalCycle(totalHeat, figDir, "seaice.seasonal.cycle", "png", "Sea ice heat content")
   
   climatology = seasonalCycle@climo
   climatology!0 = "month"
   climatology&month = ispan(1,12,1)
   devFromClimo = ndtooned(seasonalCycle@deviation)
   devFromClimo!0 = "time"
   devFromClimo&time = totalHeat&time
   
   ; save results into a netcdf
   system("/bin/rm -f " + ncSaveName)
   ncdf = addfile(ncSaveName, "c")
   
   fAtt = True
   fAtt@title = "Seasonal cycle in sea ice heat content"
   fAtt@source_file = (/arcticFname, antarcticFname/)
   fAtt@creation_date = systemfunc ("date")
   fAtt@rho = rho
   fAtt@Hfus = Hfus
   fileattdef( ncdf, fAtt )
   filedimdef(ncdf,"time",-1,True) 
   
   ncdf->heatContentTS = totalHeat
   ncdf->arcticHeat = arcticHeat
   ncdf->antarcticHeat = antarcticHeat
   ncdf->climatology = climatology
   ncdf->devFromClimo = devFromClimo

end