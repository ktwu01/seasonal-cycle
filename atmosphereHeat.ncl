procedure atmosphereHeat(figDir:string, ncSaveName:string, yearsToUse[*]:integer)
; use analysis of ERA-interim to calculate vertically integrated heat budget

begin

   ; get this file through contacting the CAS data manager: 
   ; http://www.cgd.ucar.edu/cas/catalog/reanalysis/ecmwf/erai/index.html
   ; or it is available on yellowstone
   atmosFname = "/glade/p/work/mckinnon/ERAI/ERAI.TE.1979-2015.nc"
   f = addfile(atmosFname, "r")
   time = tointeger(f->time) ; yyyymmdd
   yyyy = time/10000
   mmdd = time%10000 ; this uses the modulus operator, month is now 1201 (http://www.ncl.ucar.edu/Document/Functions/Contributed/mod.shtml)
   dd = mmdd%100 ;day is an integer, 01
   mm = mmdd/100 ;month just has the month now i.e. 12
   hh = dd
   mn = dd
   sc = dd
   hh = 0
   mn = 0
   sc = 0
   time2 = cd_inv_calendar(yyyy, mm, dd, hh, mn, sc, "days since 1979-01-01 00:00:00", 0)
   TIME = cd_calendar(time2, 0)

   ; only read in time span within yearsToUse
   idxUse = ind((TIME(:, 0).ge.min(yearsToUse)).and.(TIME(:,0).le.max(yearsToUse)))
   TE = f->TE(idxUse, :, :)
   printVarSummary(TE)
   TE&time = time2(idxUse)

   ; get weight of atmosphere
   psFname = "/glade/p/work/mckinnon/ERAI/ERAI.PS.1979-2015.nc"
   f2 = addfile(psFname, "r")
   PS = f2->PS(idxUse, :, :)
   printVarSummary(PS)
   PS&time = time2(idxUse) ; Pa

   ; to get mass of atmosphere for each gridbox, divide by g and multiply by surface area
   dA = area_global_rectilinear_grid(PS&lat, PS&lon, False)
   ; switch from km^2 to m^2
   dATemp = (10^3)^2 * dA
   copy_VarMeta(dA, dATemp)
   dATemp@units = "m^2"
   dA = dATemp
   delete(dATemp)

   massOfAtm = PS*conform_dims(dimsizes(PS), dA, (/1, 2/))/g

   ; multiply TE by mass of atmosphere to get J
   atmosEnergy = TE*massOfAtm
   copy_VarCoords(TE, atmosEnergy)
   atmosEnergy@units = "J"
   atmosEnergy@source_file = (/atmosFname, psFname/)

   ; sum over space to make time series
   atmosTS = dim_sum_Wrap(dim_sum_Wrap(atmosEnergy))

   ; and make a plot
   seasonalCycle = plotSeasonalCycle(atmosTS, figDir, "atmos.seasonal.cycle", "png", "Atmospheric heat content")

   climatology = seasonalCycle@climo
   climatology!0 = "month"
   climatology&month = ispan(1,12,1)
   devFromClimo = ndtooned(seasonalCycle@deviation)
   devFromClimo!0 = "time"
   devFromClimo&time = atmosTS&time

   ; save results into a netcdf
   system("/bin/rm -f " + ncSaveName)
   ncdf = addfile(ncSaveName, "c")

   fAtt = True
   fAtt@title = "Seasonal cycle in atmospheric heat content"
   fAtt@source_file = (/atmosFname, psFname/)
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )
   filedimdef(ncdf,"time",-1,True)

   ncdf->heatContentMap = atmosEnergy
   ncdf->heatContentTS = atmosTS
   ncdf->climatology = climatology
   ncdf->devFromClimo = devFromClimo


end