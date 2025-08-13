;+
;
; PROJECT:
;
;        CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;        mantained by Giulio Del Zanna, to develop additional features and
;        provide them to the astrophysics community.
;        Contributions via github are welcomed. 
;
; NAME:
;      doppler_factor_theta()
;       
; PURPOSE:
;
;      A function to compute the Doppler dimming factor for a line resonantly
;      photo-excited assuming uniform disk brightness and either
;      isotropic or bi-Maxwellian ion velocities, Del Zanna, 2025, MNRAS
;      If an error occurs, the function returns -1
;
; EXPLANATION:
;
; 	This function calculates the integral over the solid angle
; 	subtended by the scattering point C the integral of the disk
; 	line profile radiance I (in photons cm-2 s-1 sr-1 A-1) with
; 	the ion absorption distribution of velocity, assumed to be
; 	Maxwellian and isotropic or bi-Maxwellian. 
;       The integral, considering only the total scattered intensity,
;       reduces to an integral over wavelength of the line profile
;       radiance I and an exponential.
;       The default assumption is that the input disk radiance is an
;       averaged over the solar disk and does not depend on
;       direction/angles. However, a limb-brightening function can be
;       input. A radial velocity can also be input, to calculate the
;       Doppler dimming factor.
;
; INPUTS:
;
;       lambda_0
;               the rest wavelength of the line in Angstroms.
;
;
;       RPE
;               An IDL structure containing information to be passed
;               to the routines. The required tags currently are:
;
;        .gname
;               e.g. 'fe_13'
;
;        .model 
;               1 for Maxwellian isotropic;
;               2 for bi-Maxwellian.
;
;        if RPE.model = 1:
;
;        .Tion
;               the  ion temperature [K], from which 
;               thermal_v  is calculated:
;               The thermal velocity in cm/s, i.e. the most probable
;               speed of the Maxwellian distribution of the ion
;               velocities: 
;
;               thermal_v=sqrt((2*kb*Tion)/(mp*get_atomic_weight(Z))) 
;
;               where kb=1.38062e-16  is Boltzmann's constant
;               (cgs units) 
;               Tion  is the ion temperature [K]
;               mp=1.672661e-24 is the proton mass in grams
;               get_atomic_weight(Z) is the average atomic weight for
;               the element Z (Z=26 for Iron).
;
;  if RPE.model = 2:
;
;        .Tpar
;        .Tperp
;              The temperatures [K] of the parallel and perpendicular
;              distributions of the bi-Maxwellian, from which thermal_v
;              the averaged thermal velocity between the parallel and
;              perpendicular directions is calculated.
;
;        .r
;              the distance of the scattering point C from Sun centre
;              in solar radii units.
;
;        .u
;              The outflow velocity [cm/s] assumed to be in the radial
;              direction from Sun centre (default is 0.)
;
;        .psi 
;              The angle between  the plane of the sky and u (default
;              is 0.)
;
;        .a, .b
;              
;             The scattering factor is the angular dependence which
;             varies with the type of transition. It can be written as
;
;             a+b*(n x n')^2
;
;             where a, b are coefficients. By default, no angular
;             dependence is included, i.e. a=1 and b=0.
;
;        .radiance_ergs
;
;             The averaged disk radiance in ergs cm-2 s-1 sr-1 
;
;        .fwhm_a
;
;             The FWHM of the line disk profile in Angstroms.
;             Note: the thermal FWHM in Angstroms is
;
;             2*sqrt(alog(2))*thermal_v* lambda_a/c 
;
;             where lambda_a is the wavelength in Angstroms and c the
;             speed of light in cm/s. To this, the non-thermal FWHM should be
;             added in quadrature. The non-thermal FWHM for coronal
;             lines from Full-disk irradiance measurements is about 34
;             km/s (Feldman & Behring 1974) which in Angstroms is:
;
;             2*sqrt(alog(2))*34e5* lambda_a/c
;
;         ALTERNATIVELY, instead of radiance_ergs, fwhm_a:
;
;            .disk_lambda
;               a wavelength array (Angstroms)
;
;           .disk_spectrum
;               a disk spectrum, in photons cm-2 s-1 sr-1 Angstroms-1
;
; OPTIONAL:
;
;        RPE.lb_helio_angle, RPE.lb_values
;
;             the heliocentric angle (radians, between 0 and !pi/2.)
;             and the variation of the radiance across the disk, appropriately scaled.
;
;        RPE.radius
;
;             this is the radius of a circular source on the solar disk
;             at Sun centre, in solar radius units. By default assume the
;             whole disk. This option can be used to estimate the the RPE
;             from a region on the Sun, e.g. an active region or a flare.
;
; OPTIONAL INPUTS:
; 
;         n_theta:
;
;             number of bins in theta (default=300)
;
; OUTPUT:
;
;       The Doppler Dimming factor
;
; PROGRAMMING NOTES:
;
;       The input disk line profile is either constructed assuming a
;       Gaussian with a FWHM specified and total radiance, or taken as
;       an input spectrum. In the first instance, the bin in
;       wavelength is assumed to be fwhm_a/10  and the Gaussian is
;       constructed over a +/- 3*fwhm_a wavelength range.
;
;       In the second case, a 6 Angstrom region centred on the line is
;       taken. If the bin of the input spectrum is greater than 0.002 Angstroms, 
;       the spectrum is spline interpolated on a grid of 0.002 Angstroms. 
;
;       n is the direction toward the observer.
;       n' is the direction of the disk radiance from P' toward the
;       scattering point P.
;       u is the radial direction 
;
;       The main angles are:
;        theta, phi: the angles between n' and u as seen from P
;        psi:   the angle between  the plane of the sky and u
;
;
;
; HISTORY:
;       version 1, 24-May-2022 Giulio Del Zanna (GDZ) 
;       v.2,  5-June-2022  introduced the / sqrt(!pi) term within the integral
;       v.3 added thermal speed as array
;       v.4 Major re-write, fixing a bug in the angular part and using arrays.
;           thermal speed is now a scalar.  18 Jan 2025 GDZ
;       v.5 added regridding in 0.002 Angstroms an input spectrum. 18 Jan 2025 GDZ
;       v.6, 12 May 2025 GDZ , reorganised parts.
;
; VERSION     :  6
;
;-


function doppler_factor_theta,  lambda_0, RPE, n_theta=n_theta, verbose=verbose

  
  
; r is the distance in solar radii units of P from the source.
  r=RPE.r
  if r le 1. then begin 
     print,'% DOPPLER_FACTOR: Error in input: distance from limb is zero! '
     return,-1
  endif 
  
; this is the radius of a circular source on the solar disk at Sun centre, in
; solar radius units. By default assume the whole Sun:
  
  if tag_exist(RPE, 'radius') then radius=RPE.radius else radius=1.

  if tag_exist(RPE,'lb_helio_angle') and  tag_exist(RPE,'lb_values') then begin

; this option does not work if you have defined a radius.
     if radius ne 1 then begin 
        print,'% DOPPLER_FACTOR: Error in input: limb-brightening is NOT compatible with radius!'
        return,-1
     endif 
  endif    

  
; u is radial outflow in cm/s  
  if tag_exist(RPE, 'u') then u=RPE.u else u=0.

  if tag_exist(RPE, 'psi') then psi=RPE.psi else psi=0.


; scattering factors: this is the angular dependence which
; depends on the type of transition. If not defined, set the
; scattering factor a+b*(n x n')^2 =1 , i.e. a=1 and b=0

  if tag_exist(RPE, 'a') then a=RPE.a else a=1.
  if tag_exist(RPE, 'b') then b=RPE.b else b=0.

  
  if tag_exist(RPE, 'radiance_ergs') and tag_exist(RPE, 'fwhm_a') then begin

;  fwhm_a is FWHM in Angstrom
; radiance_ergs is the total radiance in ergs

;-- Disk intensity profile units are phot cm-2 sr-1 s-1 A-1 
     

     radiance_ergs=RPE.radiance_ergs & fwhm_a=RPE.fwhm_a
     
     if n_elements(radiance_ergs) eq 1 and n_elements(fwhm_a) eq 1 then begin

        nbin=30.
        bin_a= double(fwhm_a/nbin) ; bin in wavlength, Angstrom
        
; take a wavelength range +/- 3 FHWM
        n_lambda= fix(((lambda_0+3.*fwhm_a)-(lambda_0-3.*fwhm_a))/bin_a)+1
        wavelength_a=double(lambda_0-3*fwhm_a+bin_a*findgen(n_lambda))
        sigma=fwhm_a/2./ sqrt(2.*alog(2.))
        spectrum_photons= radiance_ergs*exp(-(wavelength_a-lambda_0)^2/2./(sigma)^2)/ sigma/ sqrt(2.*!pi)/1.9866e-8*lambda_0

; convert the bin to an array:
        bin_a= dblarr(n_lambda) & bin_a[*]= double(fwhm_a/nbin)
        
     endif  else begin
        print,'% DOPPLER_FACTOR: Error in input radiance_ergs, fwhm_a !'
        return,-1
     end
     
  endif else  if tag_exist(RPE, 'disk_lambda') and tag_exist(RPE, 'disk_spectrum') then begin

     disk_lambda=RPE.disk_lambda & disk_spectrum=RPE.disk_spectrum
     
     if n_elements(disk_lambda) ne  n_elements(disk_spectrum) then begin 
        print,'DOPPLER_FACTOR: Error, wavelengths and spectrum arrays have different sizes !'
        return,-1
     end
     
; make sure they are ordered:
     isort=sort(disk_lambda)
     disk_lambda=disk_lambda[isort]
     disk_spectrum=disk_spectrum[isort]

; take a 6 Angstrom region centred on the line of interest:
     index=where(disk_lambda ge lambda_0-3. and disk_lambda le lambda_0+3.,n_lambda)

     wavelength_a=disk_lambda[index]
     spectrum_photons=disk_spectrum[index]

; bin in wavelength (Angstroms), array 
     pp=wavelength_a[indgen(n_lambda-1)+1]-wavelength_a[indgen(n_lambda-1)] 
     bin_a= [pp, pp[n_lambda-2]]

     if min(bin_a) gt 0.002 then begin
        
        print,'% DOPPLER_FACTOR: WARNING: regridding the spectrum with a spline interpolation onto a linear grid of 0.002 Angstroms...'
        
        grid_wavelength_a=min(wavelength_a)+0.002*findgen( (max(wavelength_a)-min(wavelength_a))/0.002+1 )
        grid_spectrum=interpol(spectrum_photons, wavelength_a,grid_wavelength_a,/spline) 

        if keyword_set(verbose) then begin
           
           window,/free
           set_line_color
           plot, chars=2, wavelength_a, spectrum_photons, psym=10,/xst
           oplot, grid_wavelength_a,grid_spectrum, col=3
           al_legend,['Input spectrum','Re-gridded spectrum (interpol,/spline)'],$
                     /top,/left,chars=1.6, textcolor=[0,3]

           
        end
        
        bin_a=fltarr(n_elements(grid_wavelength_a)) & bin_a[*]=0.002

; now overwrite
        wavelength_a=temporary(grid_wavelength_a)
        spectrum_photons=temporary(grid_spectrum)
        
     end 
     
     
  endif else begin
     print,'% DOPPLER_FACTOR: Error, input not defined ! '
     return,-1
  end
  

  c=2.997925d+10                ; speed of light , cm

; this is the angle (radians) of the Sun subtended by the point in the corona C:  
  sa=asin( radius/r)   

; arc-sin in radians -->  sin(sa)= 1/r   recall that r is the distance in Rsun

  if keyword_set(verbose) then  $
     print, 'angle subtended by C (sr): '+trim(asin(radius/r))
  
  
; number of bins in theta: 
  if n_elements(n_theta) eq 0 then  n_theta=300.d                     
  
  dtheta=sa/(n_theta-1)

; create a linear grid of points:  
  array_theta=0. + findgen(n_theta)*dtheta

; calculate the distance SC as a function of theta:
  distance= (r*(cos(array_theta)>0.) - sqrt((radius^2 - (r*sin(array_theta))^2)>0.))>0.

; calculate the heliocentric angle theta':
  array_theta_p=  asin(distance/radius *sin(array_theta) )


  distance_p= (sqrt(r^2 +radius^2 -2*r*radius* cos(array_theta_p)))>0.


  cos_gamma= ((r* cos(array_theta_p)-1.)>0.)/distance_p
  
  domega=2*!pi*sin(array_theta_p)* cos_gamma / distance_p^2 

  
  if tag_exist(RPE,'lb_helio_angle') and  tag_exist(RPE,'lb_values') then begin

; check that we have numbers to interpolate:

     if min(rpe.lb_helio_angle) gt min(array_theta_p) then begin
        print,'% DOPPLER_FACTOR: Error in input: limb-brightening minimum heliocentric angle is: '+$
              trim(min(rpe.lb_helio_angle))+ " which is NOT compatible with minimum of theta'="+trim(min(array_theta_p))
        return,-1
     endif
     
     if max(rpe.lb_helio_angle) lt max(array_theta_p) then begin
        print,'% DOPPLER_FACTOR: Error in input: limb-brightening maximum heliocentric angle is: '+$
              trim(max(rpe.lb_helio_angle))+ " which is NOT compatible with maximum  of theta'="+trim(max(array_theta_p))
        return,-1         
     endif
     
; do a linear interpolation over the grid values of the heliocentric angle theta prime 
     
     limb_br_values= interpol(rpe.lb_values, rpe.lb_helio_angle, array_theta_p)

     
  endif else limb_br_values=1.+fltarr(n_theta) 

  
  
  dl=(wavelength_a-lambda_0)*1.d-8 ; delta lambda in cm - array
  
; calculate the averaged thermal speed and the first factor:
  kb=1.38062e-16
  mp=1.672661e-24

  convertname, RPE.ch_ion_name,iz,ion
; Ion mass:
  ion_mass= mp*get_atomic_weight(iz)
  
  if RPE.model eq 1 then begin
; most probable speed in cm:
     thermal_v=sqrt((2.*kb*RPE.Tion)/(ion_mass)) 

; Note: the  thermal FWHM in Angstroms is 
;                      sqrt(4.*alog(2))*wavelength (Angstroms)* thermal_v/c

     
; Replicate t1 into a 2D array with dimensions [number_theta, number_wavelengths]
     t1= replicate(1., n_theta) # (dl / (lambda_0*1e-8*thermal_v)*c)
     
     IF keyword_set(verbose) THEN $
        print, '% DOPPLER_FACTOR: isotropic thermal velocity (cm): '+string(thermal_v, format='(e9.2)')

  endif else if RPE.model eq 2 then begin

     thermal_v=sqrt( 2.*kb*RPE.Tpar/(ion_mass)* (cos(array_theta))^2 +$
                     2.*kb*RPE.Tperp/(ion_mass)* (sin(array_theta))^2  )

; now  thermal_v is an array:
     
     t1=(1./thermal_v)#(dl*c / (lambda_0*1e-8))
     
     IF keyword_set(verbose) THEN $
        print, '% DOPPLER_FACTOR: bi-Maxwellian average thermal velocity (cm): '+string(thermal_v, format='(e9.2)')


  end 

  
; Replicate t2 into a 2D array with same dimensions 
  t2= transpose(replicate(1., n_elements(dl)) # (u/thermal_v*cos(array_theta)))

; replicate the spectrum into a 2D array [number_theta, number_wavelengths]
  sp1= replicate(1., n_theta) # (spectrum_photons*bin_a)

  
; now we integrate over the wavelength to get a value for each theta
; note that F is constant if u=0
  
  F=  total(sp1*exp(- (t1-t2)^2.), 2)

; angular part:
; if alpha is the angle between the radial and the line of sight n (x)
; the scattering factor is
; (a + b*(cos(alpha)*cos(theta)-sin(alpha)*cos(phi)*sin(theta))^2 )
; however here we take psi, the angle between the plane of the sky
; and the radial outflow direction, alpha=!pi - psi

  
  sf= a + b*( (sin(psi) *cos(array_theta))^2. + 1./2* (cos(psi)*sin(array_theta))^2.)   

; now we integrate over the array of theta values:
  
  tf=int_tabulated(array_theta_p, F*limb_br_values*domega*sf/ thermal_v) 
  
  doppler_factor=c/(lambda_0*1e-8)/ sqrt(!pi)* tf 


  return, doppler_factor

end
