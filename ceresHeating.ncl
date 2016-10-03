undef("ceresHeating")
procedure ceresHeating(figDir:string, ncSaveName:string, yearsToUse[*]:integer)
; Integrate the CERES heating data to get heat

begin
  
   ; constants
   secperday = 60*60*24

   ; CERES data can be accessed and subset at https://ceres-tool.larc.nasa.gov/ord-tool/jsp/SYN1degSelection.jsp
   ; The below contains daily, TOA, all sky values from Dec 2004 to Nov 2015
   ceresFname = "/glade/p/work/mckinnon/CERES/CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed3A_Subset_200412-201511.nc"
   f = addfile(ceresFname, "r")
   toa = f->toa_net_all_daily ; Wm^-2

   ; get global TOA (W)
   dA = area_global_rectilinear_grid(toa&lat, toa&lon, False)
   ; switch from km^2 to m^2
   dATemp = (10^3)^2 * dA
   copy_VarMeta(dA, dATemp)
   dATemp@units = "m^2"
   dA = dATemp
   delete(dATemp)

   toaGlobal = dim_sum_Wrap(dim_sum_Wrap(toa*conform_dims(dimsizes(toa), dA, (/1, 2/))))
   toaGlobal@units = "W"
   toaGlobal!0 = "time"
   toaGlobal&time = toa&time
   
   ; need to switch Watts to Joules via temporal integration
   ; remove the mean over the period under the assumption that there is minimal acceleration
   ; of heating
 
   ; chop to complete years
   TIME = cd_calendar(toaGlobal&time, 0)
   ; alas, stops at Nov 2015
   startIdx = ind( (TIME(:, 1).eq.1).and.(TIME(:, 2).eq.1))
   startIdx := startIdx(0)
   endIdx = ind( (TIME(:, 1).eq.12).and.(TIME(:, 2).eq.31))
   endIdx := endIdx(dimsizes(endIdx)-1)
   
    ; remove mean of toaGlobal over full period
   toaGlobalAnom = toaGlobal(startIdx:endIdx) - dim_avg(toaGlobal(startIdx:endIdx))
   toaGlobalAnom@long_name = "Anomalies in heating"
   toaGlobalAnom@units = "W"
   toaGlobalAnom!0 = "time"
   toaGlobalAnom&time = toa&time(startIdx:endIdx)
   
   ; normalize by total m^2 to get Wm^-2, and calculate amplitude
   toaGlobal2 = toaGlobalAnom/dim_sum_Wrap(ndtooned(dA))
   toaGlobal2!0 = "time"
   toaGlobal2&time = toaGlobalAnom&time
   toaGlobal2 := calculate_monthly_values(toaGlobal2, "avg", 0, False)
   ; make climo
   toaGlobal2Climo = dim_avg_n(reshape(toaGlobal2, (/10, 12/)), 0)
   toaA = 0.5*(max(toaGlobal2Climo) - min(toaGlobal2Climo))
   print("Amplitude is " + toaA + " Wm^-2")

   TIME := cd_calendar(toaGlobalAnom&time, 0)
   yrsComplete = get_unique_values(TIME(:, 0))
   nyrs = dimsizes(yrsComplete)
   heatMon = new((/nyrs, 12/), typeof(toaGlobal))
   do i = 0, nyrs-1
      idxUse = ind( TIME(:, 0).eq.yrsComplete(i) )
      heatingAnom = toaGlobalAnom(idxUse) ; anomaly for the year, since mean already removed
      ; heating = toaGlobal(idxUse) ; heating for the year
      ; heatingAnom = heating - dim_avg(heating) ; anomalous heating for the year (zero mean)
      heat = cumsum(heatingAnom, 2)*secperday ; no missing values, so opt irrelevent
      heat!0 = "time"
      heat&time = heatingAnom&time
      heatMon(i, :) = calculate_monthly_values(heat, "avg", 0, False)
      delete(idxUse)
      delete(heatingAnom)
      delete(heat)
   end do
  
   yyyymm = yyyymm_time(tointeger(yrsComplete(0)), tointeger(yrsComplete(nyrs-1)), "integer")
   yyyy = yyyymm/100
   mm = yyyymm%100 ; this uses the modulus operator, month is now 1201 (http://www.ncl.ucar.edu/Document/Functions/Contributed/mod.shtml)
   dd = mm
   hh = mm
   mn = mm
   sc = mm
   dd = 15
   hh = 0
   mn = 0
   sc = 0
   time = cd_inv_calendar(yyyy, mm, dd, hh, mn, sc, "days since 1979-01-01 00:00:00", 0)
   
   ; chop to span yearsToUse
   delete(TIME)
   TIME = cd_calendar(time, 0)
   idxUse = ind((TIME(:, 0).ge.min(yearsToUse)).and.(TIME(:,0).le.max(yearsToUse)))
   time := time(idxUse)

   heatContentTS = ndtooned(heatMon)
   heatContentTS := heatContentTS(idxUse)
   heatContentTS!0 = "time"
   heatContentTS&time = time

   seasonalCycle = plotSeasonalCycle(heatContentTS, figDir, "ceres.seasonal.cycle", "png", "CERES integrated heating")
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
   fAtt@title = "Seasonal cycle in CERES heating (integrated)"
   fAtt@source_file = ceresFname
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )
   filedimdef(ncdf,"time",-1,True)

   ncdf->heatContentTS = heatContentTS
   ncdf->climatology = climatology
   ncdf->devFromClimo = devFromClimo

end
