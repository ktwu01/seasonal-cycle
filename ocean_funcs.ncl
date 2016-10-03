;library of oceanography functions for NCL
;many functions taken from the CSIRO SEAWATER matlab library
;see  Fofonoff, P. & Millard, R.C. Unesco 1983. Algorithms for computation of fundamental properties of seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44.

 ;constants
 pi = 3.14159
 Re = 6371220.0  ;mean radius of earth (m)
 g  = 9.80665    ;standard gravity



;list of functions
 
;fm2meter(z) - convert dept in fathoms to meters
;pres2depth(P,lat) - convert pressure (db) to depth (m)
;depth2pres(z,lat) - convert depth (m) to pressure (db)
;dT_adiab_sw - calculate adiabatic temperature gradient
;theta_sw(T,S,P,Pr) - calculate potential temp. for sea water
;cp_sw - calculate constant pressure specific heat
;sw_smow - density of Standard Mean Ocean Water
;sw_dens0 - calculate seawater density at atmos. surface pressure
;sw_seck - calculate Secant Bulk Modulus (K) of seawater
;sw_dens(T,S,P) - calculate seawater density
;sw_svan - calc specific volume anomaly (only use if you don't already have density)





;======================
;fm2meter
;=======================
 
;convert dept in fathoms to meters
 undef("fm2m")

 procedure fm2m (Z:numeric)

 begin
   
   Z = Z * 1.822
 end

;====================================
;pres2depth
;====================================

;convert pressure (db) to depth (m)
;(see Saunders, P.M., 1981. Practical conversion of pressure to depth. J. Phys. Ocean, 11, 573-574.

  undef("pres2depth")
  
  function pres2depth (P:numeric, lat:numeric)

;p   - depth (m)
;lat - latitude (deg. North); must be same size as p

  local dimd, diml,z, c1, c2, c3, c4, d2rad, gam, y, rad, gy, z, bline, tline

  begin

    ;check array sizes match
    
    dimd = dimsizes(P)
    diml = dimsizes(lat)


    if (dimsizes(dimd).ne.dimsizes(diml).or.any(dimd.ne.diml) ) then
      
      print("Error. ocean_funcs.ncl: pres2depth")
      print("P and lat arrays must be same size")
      exit
    end if

    d2rad = pi/180.
    c1    = 9.72659
    c2    = -2.2512e-5
    c3    = 2.279e-10
    c4    = -1.82e-15
    gam   = 2.184e-6

    
    y     = abs(lat)
    rad   = sin(d2rad*y)
    rad   = rad^2

    gy    = 9.780318 * (1. + (rad*5.2788e-3) + (rad^2*2.36e-5) )
   

    ;calc. depth
    bline =gy + (gam*.5*P)
    tline = (c1*P) + (c2*P^2) + (c3*P^3) + (c4*P^4)

    z = tline/bline
   

    return(z)


  end

;====================================
;depth2pres
;====================================

;convert depth (m) to pressure (db)
;(see Saunders, P.M., 1981. Practical conversion of pressure to depth. J. Phys. Ocean, 11, 573-574.

  undef("depth2pres")
  
  function depth2pres (d:numeric, lat:numeric)

;d   - depth (m)
;lat - latitude (deg. North); must be same size as d

  local dimd, diml, c1, c2, Y, deg2rad, p

  begin

    ;check array sizes match
    
    dimd = dimsizes(d)
    diml = dimsizes(lat)


    if (dimsizes(dimd).ne.dimsizes(diml).or.any(dimd.ne.diml) ) then
      
      print("Error. ocean_funcs.ncl: depth2pres")
      print("d and lat arrays must be same size")
      exit
    end if

    
    c2      = 2.21e-6

    ;calc constant c1
    deg2rad = pi/180.
    Y       = sin(abs(lat)*deg2rad)
    c1      = (5.92 + (5.25 * Y^2.)) * 1.e-3

    ;calc. pressure


    p       = ( (1-c1) - sqrt( (1-c1)^2 - (4*c2*d) )  )/(2 * c2)


    return(p)


  end


;=========================================
;dT_adiab_sw
;==========================================

;calculate adiabatic temperature gradient (viz Fofonoff & Millard, 1983)

  undef("dT_adiab_sw")

  function dT_adiab_sw(T:numeric, S:numeric, P:numeric)

;T - temperature (C)
;S - salinity (psu)
;P - pressure (db)
;all three arrays must have the same dimensions


  local dimp, dimt, dims, T68, a0, a1, a2, a3, b0, b1, c0, c1, c2, c3, e0, e1, e2, out

  begin

    ;check array sizes
    dimp = dimsizes(P)
    dimt = dimsizes(T)
    dims = dimsizes(S)
    

    if (dimsizes(dimp).ne.dimsizes(dimt).or.any(dimp.ne.dimt).or.dimsizes(dims).ne.dimsizes(dimt).or.any(dims.ne.dimt) ) then
      
      print("Error. ocean_funcs.ncl: dT_adiab_sw")
      print("T, S and P arrays must be same size")
      exit
    end if


    
    ;constants
    a0 =  3.5803E-5
    a1 =  8.5258E-6
    a2 = -6.836E-8
    a3 =  6.6228E-10

    b0 =  1.8932E-6
    b1 = -4.2393E-8

    c0 =  1.8741E-8
    c1 = -6.7795E-10
    c2 =  8.733E-12
    c3 = -5.4481E-14

    d0 = -1.1351E-10
    d1 =  2.7759E-12

    e0 = -4.6206E-13
    e1 =  1.8676E-14
    e2 = -2.1687E-16

    T68 = T * 1.00024  ;convert to 1968 temperature scale


    out = a0 + (a1 + (a2 + a3*T68)*T68)*T68 + (b0 + b1*T68)*(S-35) + ( (c0 + (c1 + (c2 + c3*T68)*T68)*T68) + (d0 + d1*T68)*(S-35) )*P + (  e0 + (e1 + e2*T68)*T68 )*P*P

    return(out)
         

  end


  

;==============================================
;theta_sw
;===============================================

;calculate potential temp. for sea water, from salinity/pressure/temp.



  undef("theta_sw")

  function theta_sw( T:numeric, S:numeric, P:numeric, Pr:numeric)

;T  - temperature (Celsius)
;S  - salinity (psu): must be same size as T and P
;P  - pressure (db) : must be same size as T and S
;Pr - reference pressure (db), scalar OR same size as T and P

  local c68, dP, dth, th, q, dimpr, dimt


  begin

    c68  = 1.00024   ;conversion constant to 1968 T scale

    ;array size check is done by function dT_adiab_sw
    dimpr = dimsizes(Pr)
    
    ;reference pressure dimensions
    if (dimsizes(dimpr).ne.1.and.dimpr(0).ne.1) then
      dimt = dimsizes(T)
      if (dimsizes(dimpr).ne.dimsizes(dimt).or.any(dimpr.ne.dimt) ) then
      
        print("Error. ocean_funcs.ncl: theta_sw")
        print("Pr array must scalar or same dimensions as T,S,P")
        exit
      end if
    end if
    

    dP = Pr - P  ;pressure difference

    
    ;1st iteration
    dth = dP * dT_adiab_sw(T, S, P)
    th  = (T * c68) + (0.5 * dth)       
    q   = dth


    ; 2nd interation
    dth = dP * dT_adiab_sw(th/c68, S, (P+ (0.5*dP)) )
    
    th  = th + ( (1- (1/sqrt(2)) )*(dth-q) )
    
    q   = ( (2-sqrt(2))*dth) + ( ( (3/sqrt(2)) - 2) * q )
    


    ;3rd iteration
    dth = dP * dT_adiab_sw(th/c68, S, (P + (0.5*dP)) )

    th  = th + ( (1 + (1/sqrt(2)) )*(dth-q) )
    
    q   = ( (2+ sqrt(2))*dth) + ( ((-3/sqrt(2)) - 2) * q )


    ;4th interation
    dth = dP * dT_adiab_sw(th/c68, S, (P + dP) )

    th  = ( th + (dth - (2*q))/6 )/ c68
    
    return(th)

  end

;********************************
;cp_sw
;*******************************


 ;calculate constant pressure specific heat (cp) for seawater, from T, P and S 

  undef("cp_sw")

  function cp_sw( T:numeric, S:numeric, P:numeric)

;T  - temperature (Celsius)
;S  - salinity (psu): must be same size as T and P
;P  - pressure (db) : must be scalar or same size as T and S


  local dimp, dimt, dims, T1,T2, T3, T4, c0,c1,c2,c3,c4,a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,d0,d1,d2,d3,d4,e0,e1,e2,f0,f1,f2,f3,g0,h0,h1,h2,j1, A, B, cp_st0, cp_0t0, d1_cp, d2_cp, cp, Pbar


  begin

 ;check array sizes

    dimp = dimsizes(P)
    dimt = dimsizes(T)
    dims = dimsizes(S)

    if (product(dimp).ne.1.and.( dimsizes(dimp).ne.dimsizes(dimt).or.any(dimp.ne.dimt))) then
      print("Error. ocean_funcs.ncl: cp_sw")
      print("P must be scalar or same size as T, S arrays")
      exit
    end if


    if (dimsizes(dims).ne.dimsizes(dimt).or.any(dims.ne.dimt)) then
      print("Error. ocean_funcs.ncl: cp_sw")
      print("T, S arrays must be same size")
      exit
    end if


;check valid ranges

    if (any(S.gt.40..or.S.lt.0.)) then
      print("Warning: ocean_funcs.ncl: cp_sw")
      print("S is outside valid range of 0-40")
    end if

    if (any(T.gt.35..or.T.lt.0.)) then
      print("Warning: ocean_funcs.ncl: cp_sw")
      print("T is outside valid range of 0-35C")
    end if

;convert P from dB to B, and convert T to 1968 T-scale

    Pbar = P/10.
    T1   = T * 1.00024


;specific heat at P = 0


    ;temperature powers
    T2 = T1^2
    T3 = T1^3
    T4 = T1^4


    ;empirical constants
    c0 = 4217.4
    c1 = -3.720283
    c2 = 0.1412855
    c3 = -2.654387e-3
    c4 = 2.093236e-5

    a0 = -7.643575
    a1 = 0.1072763
    a2 = -1.38385e-3

    b0 = 0.1770383
    b1 = -4.07718e-3
    b2 = 5.148e-5


    
    cp_0t0 = c0 +(c1*T1) + (c2*T2) + (c3*T3) + (c4*T4)
    
    A      = a0 + (a1*T1) + (a2*T2)

    B      = b0 + (b1*T1) + (b2*T2)

    cp_st0 = cp_0t0 + (A*S) + (B*S^1.5)




;pressure dependance

    a0 = -4.9592e-1
    a1 =  1.45747e-2
    a2 = -3.13885e-4
    a3 =  2.0357e-6
    a4 =  1.7168e-8

    b0 =  2.4931e-4
    b1 = -1.08645e-5
    b2 =  2.87533e-7
    b3 = -4.0027e-9
    b4 =  2.2956e-11

    c0 = -5.422e-8
    c1 =  2.6380e-9
    c2 = -6.5637e-11
    c3 =  6.136e-13


    d1_cp = (Pbar * (a0 + (a1*T1) + (a2*T2) + (a3*T3) + (a4*T4) ) ) + ( Pbar^2 * (b0 + (b1*T1) + (b2*T2) + (b3*T3) + (b4*T4) )) + (Pbar^3 * (c0 + (c1*T1) + (c2*T2) + (c3*T3) ))


    d0 =  4.9247e-3
    d1 = -1.28315e-4
    d2 =  9.802e-7
    d3 =  2.5941e-8
    d4 = -2.9179e-10

    e0 = -1.2331e-4
    e1 = -1.517e-6
    e2 =  3.122e-8

    f0 = -2.9558e-6
    f1 =  1.17054e-7
    f2 = -2.3905e-9
    f3 =  1.8448e-11

    g0 =  9.971e-8

    h0 =  5.540e-10
    h1 = -1.7682e-11
    h2 =  3.513e-13

    j1 = -1.4300e-12
    

    d2_cp = Pbar *( (S *(d0 +(d1*T1)+(d2*T2)+(d3*T3)+(d4*T4)))  + (S^1.5 *(e0 +(e1*T1) + (e2*T2))))  + \
 (Pbar^2 *( (S *(f0+(f1*T1)+(f2*T2)+(f3*T3))) + (g0*S^1.5))  ) + \
(Pbar^3 * ( (S *(h0+(h1*T1)+(h2*T2))) +(j1*T1*S^1.5)) )



    cp = cp_st0 + d1_cp + d2_cp


    return(cp)

   end

;*************************
;sw_smow
;**************************

 ;density of Standard Mean Ocean Water (pure water)

  undef("sw_smow")

  function sw_smow( T:numeric )

;T  - temperature (Celsius)


  local a0, a1, a2, a3, a4, a5, T68, dens

  begin
  
  ;coefficients
    a0 = 999.842594
    a1 = 6.793952e-2
    a2 = -9.095290e-3
    a3 = 1.001685e-4
    a4 = -1.120083e-6
    a5 = 6.536332e-9

    T68 = T*1.00024

    dens = a0 +(a1*T68) + (a2*T68^2) + (a3*T68^3) + (a4*T68^4) + (a5*T68^5)

    return(dens)

  end

;**************************
;sw_dens0
;**********************

 ;calculate seawater density at atmos. surface pressure


  undef("sw_dens0")

  function sw_dens0( T:numeric, S:numeric)

;T  - temperature (Celsius)
;S  - salinity (psu): must be same size as T 


  local dims,dimt, T68, b0, b1, b2, b3, b4, c0, c1, c2, d0, dens


;check dimension sizes

  begin

    dims = dimsizes(S)
    dimt = dimsizes(T)

    if (dimsizes(dims).ne.dimsizes(dimt).or.any(dims.ne.dimt)) then
      print("Error. ocean_funcs.ncl: sw_dens0")
      print("T, S arrays must be same size")
      exit
    end if


;constants

  b0 = 8.24493e-1
  b1 = -4.0899e-3
  b2 = 7.6438e-5
  b3 = -8.2467e-7
  b4 = 5.3875e-9

  c0 = -5.72466e-3
  c1 = 1.0227e-4
  c2 = -1.6546e-6

  d0 = 4.8314e-4


  T68 = T * 1.00024
   
  dens = S*(b0 + (b1*T68)+(b2*T68^2)+(b3*T68^3)+(b4*T68^4) ) + S^1.5*(c0 + (c1*T68)+(c2*T68^2) ) + (d0 * S^2)

  dens = dens + sw_smow(T68)

  return(dens)

  end


;**************************
;sw_seck
;************************

 ;calculate Secant Bulk Modulus (K) of seawater

  undef("sw_seck")

  function sw_seck( T:numeric, S:numeric, P:numeric)


;T  - temperature (Celsius)
;S  - salinity (psu): must be same size as T 
;P  - Pressure (db): must be same size as T

  local dims, dimt, dimp, h0, h1, h2, h3, T68, AW, k0, k1, k2, BW, e0, e1, e2, e3, e4, KW, j0, i0, i1, i2, A, m0, m1, m2, f0, f1, f2, f3, g0, g1, g2, B, K0, K, Patm

  begin


    dims = dimsizes(S)
    dimt = dimsizes(T)


    ;check dimension sizes

    if (dimsizes(dims).ne.dimsizes(dimt).or.any(dims.ne.dimt)) then
      print("Error. ocean_funcs.ncl: sw_seck")
      print("T, S arrays must be same size")
      exit
    end if


    dimp = dimsizes(P)

    if (dimsizes(dimp).ne.dimsizes(dimt).or.any(dimp.ne.dimt)) then
      print("Error. ocean_funcs.ncl: sw_seck")
      print("T, P arrays must be same size")
      exit
    end if


    ;compression terms

    T68  = T * 1.00024
    Patm = P/10. ;vonverty to mb

    h3 = -5.77905E-7
    h2 = 1.16092E-4
    h1 = 1.43713E-3
    h0 = 3.239908


    AW  = h0 + (h1*T68) + (h2*T68^2) + (h3*T68^3)


    k2 =  5.2787E-8
    k1 = -6.12293E-6
    k0 =  8.50935E-5

    BW  = k0 + (k1 + k2*T68)*T68

    e4 = -5.155288E-5
    e3 = 1.360477E-2
    e2 = -2.327105
    e1 = 148.4206
    e0 = 19652.21

    KW  = e0 + (e1 + (e2 + (e3 + e4*T68)*T68)*T68) *T68 


    ;K at stmos. pressure


    j0 = 1.91075E-4

    i2 = -1.6078E-6
    i1 = -1.0981E-5
    i0 =  2.2838E-3


    A  = AW + S*(i0 + (i1*T68) + (i2*T68^2) ) + (j0*S^1.5)


    m2 =  9.1697E-10
    m1 = 2.0816E-8
    m0 = -9.9348E-7

    B = BW + (m0 + (m1*T68) + (m2*T68^2) )*S  ;   eqn 18

    f3 =  -6.1670E-5
    f2 =  1.09987E-2
    f1 =  -0.603459
    f0 = 54.6746

    g2 = -5.3009E-4
    g1 = 1.6483E-2
    g0 = 7.944E-2

    K0 = KW + S*(f0 + (f1*T68)+ (f2*T68^2) + (f3*T68^3) )  + S^1.5*(g0 + (g1*T68) + (g2*T68^2) )   ;eqn 16



  ; K at T, S, P
    K = K0 + (A*Patm) + (B*Patm^2)   ; eqn 15
    return(K)

  end

;************************
;sw_dens
;**********************


 ;calculate seawater density from T, S, P

  undef("sw_dens")

  function sw_dens( T:numeric, S:numeric, P:numeric)


;T  - temperature (Celsius)
;S  - salinity (psu): must be same size as T 
;P  - Pressure (db): must be same size as T

  local dims, dimt, dimp, dens0, K, Patm, dens

  begin


    dims = dimsizes(S)
    dimt = dimsizes(T)


    ;check dimension sizes

    if (dimsizes(dims).ne.dimsizes(dimt).or.any(dims.ne.dimt)) then
      print("Error. ocean_funcs.ncl: sw_dens")
      print("T, S arrays must be same size")
      exit
    end if


    dimp = dimsizes(P)

    if (dimsizes(dimp).ne.dimsizes(dimt).or.any(dimp.ne.dimt)) then
      print("Error. ocean_funcs.ncl: sw_dens")
      print("T, P arrays must be same size")
      exit
    end if

;check valid ranges

    if (any(S.gt.42..or.S.lt.0.)) then
      print("Warning: ocean_funcs.ncl: sw_dens")
      print("S is outside valid range of 0-42")
    end if

    if (any(T.gt.40..or.T.lt.-2.)) then
      print("Warning: ocean_funcs.ncl: sw_dens")
      print("T is outside valid range of -2 to 40C")
    end if

    if (any(P.gt.10000..or.P.lt.0.)) then
      print("Warning: ocean_funcs.ncl: sw_dens")
      print("P is outside valid range of 0 to 10000db")
    end if



    dens0    = sw_dens0(T,S)
    
    K         = sw_seck(T,S,P)

    Patm     = P/10.  ;convert db to mb


    dens       = dens0/(1-Patm/K)

    return(dens)


    end 
;*******************************
;sw_svan
;****************************

    undef("sw_svan")
    function sw_svan(T:numeric,S:numeric,P:numeric)


;calc specific volume anomaly

;T  - temperature (Celsius)
;S  - salinity (psu): must be same size as T 
;P  - Pressure (db): must be same size as T

  local dims, dimt, dimp, dens0, dens, svan

  begin


    dims = dimsizes(S)
    dimt = dimsizes(T)


    ;check dimension sizes

    if (dimsizes(dims).ne.dimsizes(dimt).or.any(dims.ne.dimt)) then
      print("Error. ocean_funcs.ncl: sw_svan")
      print("T, S arrays must be same size")
      exit
    end if


    dimp = dimsizes(P)

    if (dimsizes(dimp).ne.dimsizes(dimt).or.any(dimp.ne.dimt)) then
      print("Error. ocean_funcs.ncl: sw_svan")
      print("T, P arrays must be same size")
      exit
    end if



    rho     = sw_dens(T, S, P)
    rho0   = sw_dens(conform(P,0.,-1), conform(P,35.,-1), P)


    svan   = (1/rho) - (1/rho0)


    return(svan)

  end


