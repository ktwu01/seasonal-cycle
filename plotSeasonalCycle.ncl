undef ("plotSeasonalCycle")
function plotSeasonalCycle(ts[*]:numeric, figDir:string, figName:string, figExt:string, title:string)
;   
; Given a time series of monthly values, 
; (1) remove the mean for each year
; (2) calculate the climatology
; (3) plot the climatology and spread around climatology
;  
begin
  
; error checks 

   d = dimsizes(ts)
   if (dimsizes(d).gt.1)
      print("Time series should be one dimensions. Yours has " + dimsizes(d) + " dimensions")
      exit
   end if
    
   nMonths = 12.
   nYrs = d/nMonths
   
   ; make sure the data is divisible by 12
   if (round(nYrs, 0).ne.nYrs)
      print("The number of data points is not divisible by 12")
      exit
   end if
   
   ; now we can use integers
   nMonths := 12
   nYrs := d/nMonths
   
   ; (1) remove the mean for each year
   ; NOTE THAT NCL DOES THIS OPPOSITE FROM MATLAB! IS ROW MAJOR!!
   tsMat = reshape(ts, (/nYrs, nMonths/))
   tsMatAnom = tsMat - conform_dims(dimsizes(tsMat), dim_avg_n_Wrap(tsMat, 1), 0)
   ; wks = gsn_open_wks("png", "test")
   ; plot  = gsn_csm_y(wks, ndtooned(tsMatAnom), False)
   
   ; (2) calculate the climatology
   tsMatAnomClimo = dim_avg_n_Wrap(tsMatAnom, 0)
   
   ; (3) plot the climatology and spread around climatology
   tsDev = tsMatAnom - conform_dims(dimsizes(tsMatAnom), tsMatAnomClimo, 1)
   range = dim_stddev_n_Wrap(tsDev, 0)
   
   X = ndtooned((/ispan(1,12,1), ispan(12,1,1)/))
   Y = ndtooned((/tsMatAnomClimo + range, tsMatAnomClimo(::-1) - range/))
   
   system("rm -f " + figDir + "/" + figName + "." + figExt)
   wks = gsn_open_wks(figExt, figDir + "/" + figName)
   
   res = True
   res@gsnDraw = False
   res@gsnFrame = False
   res@xyLineThicknessF = 2
   
   res@trXMinF = 1.
   res@trXMaxF = 12.
   res@trYMinF = min(Y)
   res@trYMaxF = max(Y)
   
   res@tiYAxisString = "Seasonal anomaly"
   res@tiXAxisString = "Month"
   res@tiMainString = title
   
   plot  = gsn_csm_xy(wks, X(0:11), tsMatAnomClimo, res)   ; create plot
   gsres                   = True                        ; poly res
   gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
   gsres@gsFillColor       = "SlateBlue"                 ; color chosen
   dummy = gsn_add_polygon(wks, plot, X, Y, gsres)
   
   draw(plot)
   frame(wks)
   
   seasonalCycle = True
   seasonalCycle@climo = tsMatAnomClimo
   seasonalCycle@deviation = tsDev
   
   return(seasonalCycle)

end
