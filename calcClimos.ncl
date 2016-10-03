; Initial code for calculation of climatologies for the components of the seasonal heat budget in
; McKinnon and Huybers, Seasonal constraints on inferred planetary heat content

; Need to either put all scripts in current directory, or add the path to your scripts

; Data sources are listed in each function call

load "./area_rectlinear_grid.ncl"
load "./plotSeasonalCycle.ncl"
load "./ocean_funcs.ncl"
load "./argo.ncl"
load "./atmosphereHeat.ncl"
load "./seaiceHeat.ncl"
load "./cesmHeat.ncl"
load "./ceresHeating.ncl"
load "./landHeat.ncl"
load "./infill0.ncl"
load "./infillGlobalAvg.ncl"
load "./infillRegionalAvg.ncl"
 
begin
 
ncDir = "LOCATION WHERE OUTPUT NETCDF FILES SHOULD BE SAVED"
figDir = "LOCATION WHERE FIGS SHOULD BE SAVED"
; (0) constants
; cp = 3990 ; J kg^-1 K^-1, water (Trenberth and Fasullo, 2008, JClim)
; a = 6.37122e6 ; m, radius of the earth
; pi = 4.*atan(1.0)
; secperday = 24*3600
; SA = 4*pi*a^2
; rho0 = 1026.5 ; kg m^-3 (Trenberth and Fasullo, 2008, JClim)
; H_fus = 334e3 ; J kg^-1, energy required to melt ice
; g  = 9.80665  ;standard gravity

; time span: 2005-2014 based on data availability
yearsToUse = ispan(2005, 2014, 1)

argo(figDir, ncDir + "/argo.unfilled.nc", yearsToUse)
atmosphereHeat(figDir, ncDir + "/atmos.erai.nc", yearsToUse)
seaiceHeat(figDir, ncDir + "/seaice.nc", yearsToUse)
cesmHeat(figDir, ncDir + "/cesm.heat.nc", yearsToUse)
ceresHeating(figDir, ncDir + "/ceres.heat.nc", yearsToUse)
landHeat(figDir, ncDir + "/land.heat.nc", yearsToUse)

; get mask for where argo + cesm don't overlap
dblquote = integertochar(34)
cesmFname = ncDir + "/cesm.heat.nc"
argoFname = ncDir + "/argo.unfilled.nc"
maskSaveName = ncDir + "/argo.missing.mask.csv"
system("matlab -nosplash -nodesktop -r " + dblquote + "makeArgoMask('" + cesmFname + "', '" + argoFname + "', '" + maskSaveName + "')" + dblquote)

; various infilling techniques
; (1) no anomaly in unmeasured regions
infill0(ncDir + "/argo.unfilled.nc", ncDir + "/argo.missing.mask.csv", ncDir + "/argo.filled.zeros.nc", figDir)

; (2) unmeasured = global mean
infillGlobalAvg(ncDir, ncDir + "/argo.unfilled.nc", ncDir + "/argo.missing.mask.csv", ncDir + "/argo.filled.global.nc", figDir)

; (3) unmeasured = hemispheric / tropical mean
infillRegionalAvg(ncDir + "/argo.unfilled.nc", ncDir + "/argo.missing.mask.csv", ncDir + "/argo.filled.regional.nc", figDir)


end