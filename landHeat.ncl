undef("landHeat")
procedure landHeat(figDir:string, ncSaveName:string, yearsToUse[*]:integer)
; Calculate seasonal cycle of heat content associated with heating of the land surface

begin

   ; Data is from http://berkeleyearth.lbl.gov/auto/Global/Complete_TAVG_complete.txt
   ; Subset used in available in repo under 'data' folder
   landFnameAnom = "./data/Tavg.BEST.monthly.txt"
   landFnameClim = "./data/Tavg.BEST.climo.txt"
 
   Tanom = asciiread(landFnameAnom, -1, "float")
   ncol = 12
   nrow = dimsizes(Tanom)/ncol
   TanomMat = onedtond(Tanom, (/nrow, ncol/))

   yyyy = tointeger(TanomMat(:, 0))
   mm = tointeger(TanomMat(:, 1))
   dd = mm
   hh = mm
   mn = mm
   sc = mm
   dd = 15
   hh = 0
   mn = 0
   sc = 0 
   delete(sc@_FillValue) 

   time = cd_inv_calendar(yyyy, mm, dd, hh, mn, sc, "days since 1979-01-01 00:00:00", 0)

   anomVals = TanomMat(:, 2)

   climo = asciiread(landFnameClim, -1, "float")

   ; cut off complete years
   startIdx = ind( mm.eq.1 )
   startIdx := startIdx(0)
   endIdx = ind( mm.eq.12 )
   endIdx := endIdx(dimsizes(endIdx)-1)

   anomVals := anomVals(startIdx:endIdx)
   time := time(startIdx:endIdx)
   TIME = cd_calendar(time, 0)
   
   yrs = get_unique_values(TIME(:,0))
   nyrs = dimsizes(yrs)
   yrStart = min(TIME(:, 0))
   yrEnd = max(TIME(:, 0))
      
   ; add climatology back in
   totalVals = ndtooned(reshape(anomVals, (/nyrs, 12/)) + conform_dims((/nyrs, 12/), climo, 1))
   
   ; switch interpolate to daily data   
   doy = day_of_year(tointeger(TIME(:, 0)), tointeger(TIME(:, 1)), tointeger(TIME(:, 2)))
   doy2 = doy + 365*ndtooned(conform_dims((/nyrs, 12/), ispan(0, 11, 1), 0))
   doyNew = ispan(doy2(0), doy2(dimsizes(doy2)-1), 1)
   Tdaily = linint1_n(doy2, totalVals, False, doyNew, 0, 0)
   timeDaily = doyNew - 1 + cd_inv_calendar(tointeger(yrStart), 1, 1, 0, 0, 0, "days since 1979-01-01 00:00:00", 0)
   timeDaily@units = "days since 1979-01-01 00:00:00"
   TIMEdaily = cd_calendar(timeDaily, 0)
   
   ; fit a heat model for each year (1D heat conduction equation)
   dz = 1; % 1 meter
   nz = 20/dz; % go to 20 meters
   dt = 24*3600; % time step is one day
   K =  5e-7; % Hartmann (1994)
   maxiter = 500
   Tvertical = new((/nz + 1, dimsizes(Tdaily)/), "float")

   do i = 0, nyrs-1
      vals := Tdaily(ind(TIMEdaily(:,0).eq.yrs(i)))
      t := timeDaily(ind(TIMEdaily(:,0).eq.yrs(i)))
      
      nt = dimsizes(vals) - 1
      Tmat := new((/nz + 1, nt + 1/), "float")
      Tmat = 10
      
      Tmat(0, :) = vals ; surface temperature through time
      
      tol = 1e-4
      errVal = 10.
      counter = 0
      do while(errVal.gt.tol)
         Tlast := Tmat ; last guess
         Tmat(:, 0) = Tlast(:, nt)
         
         do j = 1, nt
            
            d2Tdz2 = (Tmat(0:nz-2, j-1) - 2*Tmat(1:nz-1, j-1) + Tmat(2:nz, j-1))/dz^2
            dTdt = K*d2Tdz2
            Tmat(1:nz-1, j) = dTdt*dt + Tmat(1:nz-1, j-1)
            Tmat(nz, j) = Tmat(nz-1, j)
         
         end do 
         errVal = max(abs(ndtooned(Tmat) - ndtooned(Tlast)))
         ; print(errVal)
         counter = counter + 1
         if counter.gt.maxiter
            print("Equation not converging")
         end if
         
      end do
      delete(errVal)
      
      ; remove the annual mean
      Tvertical(:, ind(TIMEdaily(:,0).eq.yrs(i))) = Tmat - conform_dims(dimsizes(Tmat), dim_avg_n(Tmat, 1), 0)
      
   end do
 
   
   ; get global average volumetric heat capacity of soil, assuming organic
   cv = 2.5e6 ; J m^-3 K^-1
   saLand = 150e6 ; surface area of land (km^2)
   dV = dz*saLand*(1e3)^2 ; m^3
   
   heatContentDaily = cv*dV*dim_sum_n(Tvertical, 0)
   heatContentDaily!0 = "time"
   heatContentDaily&time = timeDaily
   heatContent = calculate_monthly_values(heatContentDaily, "avg", 0, False)

   ; remove the ends because they are not full months
   heatContent := heatContent(1:dimsizes(heatContent)-2)

   ; pull out complete years
   TIMEnew = cd_calendar(heatContent&time, 0)
   idxStart = ind((TIMEnew(:,1).eq.1).and.(TIMEnew(:,0).eq.min(yearsToUse)))

   idxEnd = ind((TIMEnew(:,1).eq.12).and.(TIMEnew(:,0).eq.max(yearsToUse)))
   heatContent := heatContent(idxStart:idxEnd)
   
   ; plot
   seasonalCycle = plotSeasonalCycle(heatContent, figDir, "land.seasonal.cycle", "png", "Land heat content")
   climatology = seasonalCycle@climo
   climatology!0 = "month"
   climatology&month = ispan(1,12,1)
   devFromClimo = ndtooned(seasonalCycle@deviation)
   devFromClimo!0 = "time"
   devFromClimo&time = heatContent&time

   ; save results into a netcdf
   system("/bin/rm -f " + ncSaveName)
   ncdf = addfile(ncSaveName, "c")

   fAtt = True
   fAtt@title = "Seasonal cycle in land heat content"
   fAtt@source_file = (/landFnameAnom, landFnameClim/)
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )
   filedimdef(ncdf,"time",-1,True)

   ncdf->heatContentTS = heatContent
   ncdf->climatology = climatology
   ncdf->devFromClimo = devFromClimo



end


