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
;      CH_LOAD_ION_RATES
;       
; PURPOSE:
;          a function that returns a structure with the 
;          main matrices  of the rates within one ion, built
;          from the CHIANTI v.8 format files.
;
; CATEGORY:
;      CHIANTI
;
; INPUTS:
; 
;       INPUT     can either be the structure read by ch_setup_ion or
;                 just the CHIANTI name of an  ion, e.g. 'fe_13'
;       
;      	TEMP	 temperature(s) [K]. Can be either a single scalar or
;      	         an array of numbers.
;
; OPTIONAL INPUTS:
;
;       n_lev   the maximum number of levels to inclue in the model
;               ion. Note that this option is disallowed if the ion
;               has autoionising levels. 
;
;	RADTEMP	Specify background radiation temperature (default: 6000 K)
;
;	RPHOT   Distance from the centre of the star in stellar radius units.
;               I.e., RPHOT=1 corresponds to the star's surface. (Default is
;               infinity, i.e., no photoexcitation.)
;
;       RADFUNC         The name of a user-defined function that will generate
;                       a radiation spectrum as a function of temperature. 
;                       This radiation field will replace the black-body that
;                       is assumed when using the RADTEMP keyword in the call
;                       to pop_solver.
;
;       SUM_MWL_COEFFS  An array of coefficients of the same length as 
;                       the array of temperatures. Electron and proton rate 
;                       coefficients will be calculated at each temperature 
;                       and then a weighted sum of the coefficients is 
;                       performed using SUM_MWL_COEFFS. This allows 
;                       non-Maxwellian energy distributions to be 
;                       incorporated into the level balance equations.
;                       This keyword is not compatible with the PRESSURE
;                       keyword.
;
;	PATH	If specified, the routine will look for the atomic data in 
;		the PATH directory, rather than in the CHIANTI database
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
; KEYWORDS:
;
;       NO_AUTO:   do not include autionising levels even if they are
;                 present in the ion directory.
;
;       NOIONREC: If set, then level-resolved ionization/recombination
;                rates are not read for the ion, even if they exist.
;
;       NO_RREC: If set, do not read the level-resolved radiative
;                recombination (RR) files.
;
; OUTPUT:
;       A structure with the following tags:
;           n_levels  No. of levels in model.
;           aa    A-values (2D array)
;           aax   Photoexcitation/stimulated emission (2D array).
;           qq    Electron rate coefficients (3D array: n_temperatures,n_levels,n_levels)
;           ppr   Proton rate coefficients (3D array).
;
;           temp  Array of temperatures (the input TEMP).
;           ion_data
;                 The structure with the atomic data read by
;                 CH_SETUP_ION.
;           sumtst Takes value 1 or 0 if sum_mwl_coeffs has been set
;                  or not.
;           sum_mwl_coeffs Array with same size as TEMP. Contains
;                  SUM_MWL_COEFFS if it has been set,
;                  otherwise contains all 1's.  
;
; PROGRAMMING NOTES:
;
; CALLS:
;      CH_SETUP_ION, CONVERTNAME, ZION2FILENAME, PROTON_DENS
;
;
; HISTORY:
;
;        Version 1,  written by Giulio Del Zanna (GDZ)  12 May 2018 
;        Version 2, fixed a major bug inroduced when switching the 
;                  indices of the qq and ppr arrays. GDZ 14 May 2018 
;        Version 3, fixed another major bug when populating the qq
;                   array. GDZ - 21 May 2018
;        Version 4, Peter Young, 9 July 2018, various changes.
;        Version 5, GDZ,  19 July 2018, 
;        various new modifications: changed name;  reinstated the
;        keyword  no_auto; reinstated a a warning if the reduced
;        levels are for ions with autoinizing levels; removed the
;        reading of the total ionization and recombination rates as
;        they only need to be included in CH_COMBINE_RATES (no need to
;        read extra rates not used later on).
;        Version 6, GDZ, 24 Aug 2018.
;                   reinstated that the input can either be the
;                   structure read by ch_setup_ion or just the CHIANTI
;                   name of an  ion, e.g. 'fe_13'
;
;        v.7,       GDZ, 4 Oct 2018 
;                  added keywords to be passed to ch_setup_ion. The
;                  main ones are those where a wavelength range can be
;                  specified.  If there are no lines, the routine
;                  exits. 
;        v.8,      GDZ, 18 Nov 2018 
;                  Reinstated the option of 'dielectronic' files
;        v.9,      14-Dec-2018 GDZ, added NOIONREC, NO_RREC keywords.
;        v.10, 5-Mar-2019, Peter Young
;          Added sumtst and sum_mwl_coeffs to the output structure as
;          they are needed by ch_load_2ion_rates.
;
;        v.11, 12-May-2023, Peter Young
;          Added /quiet for call to proton_dens; message about flipping
;          transitions only printed if /verbose set.
;
;        v.12, 12-Feb-2025, Peter Young
;          Added a comment about the "mult" paramter. No change to code.
;
;        v.13, GDZ,  25-May-2022 
;              added optional input RPE, to add resonant photo-excitation
;
;        v.14 , GDZ Major re-write, adding the scattering factors 
;           18 Jan 2025 
;
;        v.15 , GDZ fixed a few minor issues/bugs  18 Jan 2025 
;        v.16, GDZ, 16 May 2025 changed the input to the RPE
;
;        v.17, GDZ, 12 Aug 2025, corrected the values of a few constants..
;
; VERSION     :  17,  12 Aug 2025
;
;-

function ch_load_ion_rates, input, temp, n_lev=n_lev, $
                            radfunc=radfunc, abund_file=abund_file, ioneq_file=ioneq_file, $ 
                            PATH=PATH, NOPROT=NOPROT, RPHOT=RPHOT, RADTEMP=RADTEMP,$
                            sum_mwl_coeffs=sum_mwl_coeffs, no_auto=no_auto, verbose=verbose,$
                            wvlmin=wmin,wvlmax=wmax, index_wgfa=anylines,obs_only=obs_only,$
                            noionrec=noionrec, no_rrec=no_rrec, RPE=RPE_input, error=error

  
  IF n_params() LT 2 THEN BEGIN
     print,'Use:  IDL> output = ch_load_ion_rates( input, temp )'
     print,''
     print,'  input is the CHIANTI name of ion (e.g., "o_6" )'
     return,-1
  ENDIF 

 IF n_elements(radfunc) NE 0 and n_elements(RPE_input) ne 0 THEN begin 
     print,'% CH_LOAD_ION_RATES: **** WARNING: both RADFUNC and RPE options have been set ! '
     print,'% CH_LOAD_ION_RATES: **** WARNING: make sure there is consistency.  '
  end
 
  
  IF n_elements(abund_file) EQ 0 THEN abund_file=!abund_file

; If INPUT is not a structure then assume it is an ion name (e.g.,
; 'fe_13'). 
;
IF n_tags(input) EQ 0 THEN BEGIN
  gname=input

     convertname, gname,iz,ion
     zion2filename,iz,ion,filename,name=name 

if n_elements(verbose) eq 0 then verbose=0
if keyword_set(verbose) then quiet=0 else quiet=1

     ion_data=ch_setup_ion(gname,rphot=rphot,radtemp=radtemp,noprot=noprot, $
                        ioneq_file=ioneq_file,abund_file=abund_file,path=path, $
                           quiet=quiet,  wvlmin=wmin,wvlmax=wmax, index_wgfa=anylines,$
                           obs_only=obs_only,noionrec=noionrec ) ; GDZ

  endif else begin 
ion_data=input

gname=ion_data.gname 

end 


if  is_number(ion_data) then begin 
if verbose then print,'% CH_LOAD_ION_RATES:  no lines in the selected range !'
anylines=-1

return, -1 
endif 


ion2spectroscopic,gname,snote, dielectronic=dielectronic

; for V.9, some ions still have the older 'd' directories.
; check that we either have the old v.8 files or the new v.9 files:
if dielectronic and  tag_exist(ion_data,'autostr') then begin 
print,'% ch_load_ion_rates:  ERROR ! '
anylines=-1
return, -1 
endif 

;  ioneq_file=ion_data.ioneq_file

;
; abund_file is not in ION_DATA by default, so need to check.
;
  IF tag_exist(ion_data, 'abund_file') THEN BEGIN
     abund_file=ion_data.abund_file
  ENDIF ELSE BEGIN
     abund_file=!abund_file
  ENDELSE 


  t=double(temp)
;
; need the following to turn t into an array if it only has 1 element
;
  IF n_elements(t) EQ 1 THEN BEGIN
     t0=t
     t=dblarr(1)
     t[0]=t0
  ENDIF
  nt=N_ELEMENTS(temp)           ; no. of temperatures


;
; SUMTST is a switch indicating whether SUM_MWL_COEFFS has been
; specified. 
;
  IF n_elements(sum_mwl_coeffs) EQ 0 THEN BEGIN
     sum_mwl_coeffs=dblarr(nt)+1.
     sumtst=0
  ENDIF ELSE BEGIN
     sumtst=1
  ENDELSE 


; ------
; The ion_data structure has a set of "core" tags and a set of "optional"
; tags. The core tags are:
;
  gname = ion_data.gname
  jj= ion_data.jj
  ecm= ion_data.ecm
  ecmth= ion_data.ecmth
  wvl= ion_data.wvl
  a_value= ion_data.a_value
  splstr= ion_data.splstr

; GDZ****  a fundamental tag is the level index in the energy file !
  level_index= ion_data.elvlcstr.data.index
  
; ------
; The optional tags are given below.
;
  IF tag_exist(ion_data, 'radtemp') THEN  radtemp= ion_data.radtemp
  IF tag_exist(ion_data, 'dilute') THEN  dilute= ion_data.dilute
  IF tag_exist(ion_data, 'prot_struc') THEN  prot_struc= ion_data.prot_struc
;
; The following is retained for backwards compatibility, as pe_ratio
; is no longer computed by ch_setup_ion.
;
  IF tag_exist(ion_data, 'pe_ratio') THEN  pe_ratio= ion_data.pe_ratio else $
     pe_ratio=proton_dens(alog10(t),abund_file=abund_file, ioneq_file=ioneq_file,/quiet)
;
; If pe_ratio is in the ion_data structure, then make sure it has the
; right size:
;
  IF n_elements(pe_ratio) NE n_elements(t) THEN BEGIN
     print,'% WARNING, pe_ratio size does not match temp. It will be recomputed.'
     print,'             ',n_elements(pe_ratio),n_elements(t)
     pe_ratio=proton_dens(alog10(t), ioneq_file=ioneq_file, abund_file=abund_file,/quiet)
  ENDIF


  IF tag_exist(ion_data, 'ionrec') THEN BEGIN
     ionrec_struc= ion_data.ionrec
     status=ionrec_struc.status
  ENDIF ELSE BEGIN
     status=0
  ENDELSE

;
; Note that ECM contains the best energies, so no need to check
; for "missing" levels.
;
  ecmc=ecm
  n_elvl=n_elements(ecm)


  IF keyword_set(verbose) THEN print,'% CH_LOAD_ION_RATES: the '+gname+' energy file has '+trim(n_elvl)+' levels'

  ; Multiplicity:
  ; PRY, 12-Feb-2025: the variable "mult" has been wrongly named as it 
  ; is not the multiplicity, but the weight.
  mult=2.*jj+1.
;
; hck=1.98648d-16/1.38062d-16
  hck=1.986446d-16/1.380649d-16 ; GDZ corrected values
  
  ryd2cm=109737.31534d
;

; Now we check the sizes of the files, and remove levels that are not 
; populated.
; NOTE: for the ions with autoionizing levels, we need to retain all of
; them, as many of them are not directly populated.
;**********************************************************************

;
; In the following we set N_LEVELS, which sets the size of the atomic
; model. It's possible that different data files have different
; numbers of levels, so we have to check the various arrays to find
; the minimum model size.
;
; GDZ: 
;
; The optional ion_data N_LEV allows the user to set the number of
; levels. This should not be done for the ions with autoionizing
; levels, as the satellite lines would be removed.  If only some of
; the levels are asked to be removed, then the routine removes all of
; them. If the routine is asked to have a maximum level within the
; number of bound levels, then the autoinizing levels are removed.
;

  n_wgfa=max([ion_data.wgfastr.lvl1,ion_data.wgfastr.lvl2])
  usize=max([splstr.data.lvl1,splstr.data.lvl2])

  IF tag_exist(ion_data,'autostr') THEN asize=max(ion_data.autostr.lvl2) ELSE asize=usize
  ausize=max([asize,usize])
  n_levels=min([n_elvl,ausize,n_wgfa])

; this is the user-defined cut in the number of levels:
  IF n_elements(n_lev) NE 0 THEN begin 

     IF tag_exist(ion_data,'autostr') and n_lev lt n_levels THEN begin 

        print,'% CH_LOAD_ION_RATES: WARNING: the model ion should **not** be reduced due to the'+$
              ' presence of autoionising states '

; count how many autoionizing levels we have: 
        lvl_auto= get_uniq(ion_data.autostr.lvl2, count=n_auto)
        n_bound_levels=n_levels-n_auto

; we do not want to have only some of autoionising states. GDZ 
        if n_lev gt n_bound_levels  then begin 
           print,'% CH_LOAD_ION_RATES: WARNING: the model ion should **not** be reduced when  autoionising states are present '
           print,'% CH_LOAD_ION_RATES: WARNING: all the autoionising states have been removed - keeping bound levels. ' 

           n_levels= n_bound_levels

; remove the autoionizing levels 
           ion_data= rem_tag(ion_data, 'autostr')

        endif else n_levels=min([n_levels,n_lev])

     endif  else n_levels=min([n_levels,n_lev])

  endif                         ; user-defined 

;
  IF keyword_set(verbose) THEN $
     print,'% CH_LOAD_ION_RATES: the model for '+gname+' contains '+trim(n_levels)+' levels.'


;;------------------------------[]
; Now we start filling up the rate arrays. Note the ordering of the
; indices:
;
; e.g., aa(0,19) will be zero (no 0 -> 19 A value)
;       aa(19,0) will be non-zero
;
;       qq(*,0,19) electron excitation
;       qq(*,19,0) electron de-excitation
;;------------------------------[]

;
; Create the AA and AAT (i.e., AA-transpose) arrays from the A_VALUE ion_data.
;
  aat=a_value(0:n_levels-1,0:n_levels-1) ; transpose of aa
  aa=DOUBLE(TRANSPOSE(aat))

;
; Create 2D energy arrays from the 1D energy vector.
;
; DELTA_E is defined such that:
;    delta_e[i,j] > 0 implies level i is energetically higher than j
;                     -> therefore aa[i,j] is non-zero
;
  ident=make_array(n_levels,val=1.) ; use for making de arrays
  ecmn=ecmc(0:n_levels-1)           ; n_levels
  den=ecmn#ident & dem=ident#ecmn
  delta_e=den-dem

;
; PRY, 27-Apr-2016
; The following code checks to see if there are A-values that have
; been assigned to the inverse transition. That is, if a level i is
; energetically lower than a level j, yet AA[i,j] is non-zero. If such
; cases exist then we transfer the A-value into the transpose matrix (AAT)
;
  i=where(delta_e LT 0.,ni)
  j=where(aa[i] GT 0.,nj)
  IF nj NE 0 THEN BEGIN
     print,'%WARNING: ['+trim(gname)+'] - '+trim(nj)+$
           ' A-values have been assigned to inverse transitions. These have been flipped.'
     k=i[j]
                                ;
                                ; The following prints out the transitions that have been
                                ; flipped. 
                                ;
     IF keyword_set(verbose) THEN BEGIN 
        FOR ii=0,nj-1 DO BEGIN
           ij=get_ij(k[ii],n_levels)
           print,format='(2i5,e15.3)',ij[0]+1,ij[1]+1,aa[k[ii]]
        ENDFOR
     ENDIF 
                                ;
     aa_temp=make_array(n_levels,n_levels,/double)
     aa_temp[k]=aa[k]
     aa_temp_transpose=transpose(aa_temp)
     aa=aa + aa_temp_transpose - aa_temp
     aat=transpose(aa)
     aa_temp=0                  ; tidy up
     aa_temp_transpose=0        ;
  ENDIF 

 
  
;;-----------------------------------------------------------------
; The following loads up the photoexcitation (pexc) and stimulated 
; emission (stem) arrays. These are combined and put into the AAX
; array. 
;
  aax=dblarr(n_levels,n_levels)

  IF N_ELEMENTS(dilute) EQ 0 THEN dilute=0.

  if dilute ne 0 then begin 
; define the arryas:
     stem=dblarr(n_levels,n_levels) & pexc=dblarr(n_levels,n_levels) & ede=dblarr(n_levels,n_levels)

                                ;
     multn=mult[0:n_levels-1]   ; in case mult and ecm are bigger than 
                                ;
     mm=TRANSPOSE(multn#(1/multn)) ; needed for photoexcitation
                                ;
                                ; Note below that there's a factor of 1d8^3. This is because the units
                                ; of lambda are angstroms, but I need a cm^-3 in the units of the energy
                                ; density.
                                ;
     IF n_elements(radfunc) NE 0 THEN BEGIN
        en_rf=abs(den-dem)
        i=where(en_rf EQ 0.)           ; the index i prevents underflows
        IF i[0] NE -1 THEN en_rf[i]=1d50 ; set to arbitrarily large value
                                ;
        lmbda=1d8/en_rf
                                ;
        bits=str_sep(radfunc,',')
        CASE n_elements(bits) OF
           1: result=call_FUNCTION(radfunc,lmbda)
           2: BEGIN
              rfunc=bits[0]
              a1=double(bits[1])
              result=call_FUNCTION(rfunc,lmbda,a1)
           END
           ELSE: BEGIN
              rfunc=bits[0]
              a1=double(bits[1])
              a2=double(bits[2])
              result=call_FUNCTION(rfunc,lmbda,a1,a2)
           END
        ENDCASE
                                ;
        IF i[0] NE -1 THEN result[i]=0d0
                                ;
        result=result*lmbda^5/(1d8)^3/8d0/!pi/1.986d-8
     ENDIF ELSE BEGIN
        dd=ABS(den-dem)*hck/radtemp
                                ;
                                ; the following lines are necessary to prevent infinities and floating
                                ; underflow errors
                                ;
        dd=dd < 150.
        i=where(dd EQ 0.)
        j=where(dd LE 1d-15 AND dd NE 0.)
        k=where(dd GE 1d-15)
                                ;
        ede[k]=exp(dd[k]) - 1.
        IF j[0] NE -1 THEN ede[j]=dd[j]
        ede[i]=1d50             ; arbitrarily large value
                                ;
        result=1d0/ede
        result[i]=0d0           ; set i's to zero since A-values are zero (thus don't
                                ; contribute to pexc and stem
     ENDELSE
                                ;
     ind=WHERE( (aat NE 0.) AND (result NE 0.) )
     IF ind[0] NE -1 THEN pexc[ind]=aat[ind]*dilute*mm[ind]*result[ind]
                                ;
     ind=where( (aa NE 0.) AND (result NE 0.) )
     IF ind[0] NE -1 THEN stem[ind]=aa[ind]*dilute*result[ind]
                                ;
     aax=pexc+stem
     
  ENDIF
  
;;----------------------------------------------------------------------------
; NEW section to calculate RPE for a single transition using various methods:
;;----------------------------------------------------------------------------
; Giulio Del Zanna 

  
  
 if  n_elements(RPE_input) gt 0 then begin

; copy the structure that is going to be modified below.
    
RPE=RPE_input
    
;-- Constants --
c=2.997925e+10
kb=1.38062e-16
h=6.626205e-27
mp=1.672661e-24
; rs=7.e+10

    
    if RPE.model eq 1 or  RPE.model  eq 2 then begin
       
 if tag_exist(RPE, 'radiance_ergs') and tag_exist(RPE, 'fwhm_a') then begin

    if rpe.ch_ion_name eq gname then begin
       
; include RPE only for  the ion and line requested !
;---------------------------------------------------
    
; wavelength in vacuum at rest, Angstroms: 
       lambda_0=input.wvl[RPE.lower_level-1,RPE.upper_level-1]
; 1.d8 /abs(ecmc[RPE.upper_level-1]-emcm[RPE.lower_level-1])

; The A-value is  aa[RPE.upper_level-1, RPE.lower_level-1]
; also input.a_value[RPE.lower_level-1,RPE.upper_level-1]
    
; the wavelength in cm is 1/DE
       
;IF keyword_set(verbose) THEN $
;print, '% CH_LOAD_ION_RATES: thermal FWHM (Angstroms): ',sqrt(4.*alog(2))*lambda_0 * thermal_v/c
       
; first find the IDL indices of the levels:
; Note: the level ordering could be strange.  level_index has the full list - 

ind_lower=where( level_index eq RPE.lower_level,n_lower)
       if n_lower ne 1 then message, '% CH_LOAD_ION_RATES: ERROR in input lower level number ! '
; fix this into a scalar:
ind_lower=ind_lower[0]

       ind_upper=where( level_index eq RPE.upper_level,n_upper)
       if n_upper ne 1 then message, '% CH_LOAD_ION_RATES: ERROR in input upper level number ! '
       ind_upper=ind_upper[0]
       
; calculate Einstein absorption coefficient B_lu coefficient for the transition:
; B_lu= gu/gl*Aul * (lambda*1e-8)^3 /(2*h*c) 
       
    B_lu= mult[ind_upper]/mult[ind_lower]* $
          aa[RPE.upper_level-1, RPE.lower_level-1]* (lambda_0*1d-8)^3./(2.*h*c)

IF keyword_set(verbose) THEN $
print, '% CH_LOAD_ION_RATES: B_lu:'+string(B_lu, format='(e10.2)')
    
    if not tag_exist(RPE, 'a') then begin 
    
; Calculate the scattering factor for the transition
; See Hamilton (1947) and Del Zanna (2024)
    
    j_l=jj[ind_lower] & j_u=jj[ind_upper]
    
    Delta_j = j_u -j_l

    case Delta_j of 
       1: begin
          E_1= 1./ 10 *  (26. *j_l^2+ 37 * j_l+10) / ((j_l+1)*(2*j_l+1)) 
          E_2= 1./ 10 *  ((2*j_l+5)*(j_l+2)) /  ((j_l+1)*(2*j_l+1))             

       end 
       0:begin
          E_1= 1./ 10 *  (12.* j_l^2+12*j_l+1)  / (j_l*(j_l+1))
          E_2= 1./ 10 *  ((2.*j_l-1) *(2*j_l+3)) / (j_l*(j_l+1))
       end 
       -1:begin
          E_1= 1./ 10 *  (26. * j_l^2+ 15*j_l-1) / (j_l*(2*j_l+1))
          E_2= 1./ 10 *  ((2.*j_l-3)*(j_l-1))  /  (j_l*(2*j_l+1))
       end       
            ELSE: BEGIN 
               print, '% CH_LOAD_ION_RATES: ERROR? Delta J of the requested transition is  ',delta_j
               print, '% CH_LOAD_ION_RATES:  setting coefficients a=1, b=0'
               E_1=4./3 & E_2=0.
            END            
       endcase

; add the info to the structure:
    RPE=add_tag(RPE, 3./4*E_1, 'a')
    RPE=add_tag(RPE, 3./4*E_2, 'b')

    
    if keyword_set(verbose) then print,'% CH_LOAD_ION_RATES, added RPE scattering: a='+trim(3./4*E_1)+$
                                       ' b='+trim(3./4*E_2)

 endif 
    
; Calculate the Doppler factor:
    doppler_factor=doppler_factor_theta(lambda_0,  RPE, verbose=verbose)

    if doppler_factor eq -1 then begin 
       error_msg='% CH_LOAD_ION_RATES: ERROR, RPE  model not recognised ! '
   print,error_msg
   return,-1
end 

; Neglect the stimulated emission.
; **** Note: the 1/(4.*!pi) factor is *not* included in pop_solver.

; include the RPE rate into the matrix:          
      aax[ RPE.lower_level-1, RPE.upper_level-1] = (lambda_0*1.d-8)*B_lu *h /(4.*!pi)*doppler_factor

      IF keyword_set(verbose) THEN $
print, ' RPE rate: ', string((lambda_0*1.d-8)*B_lu *h /(4.*!pi)*doppler_factor, format='(e10.2)')
      
   endif 
      
 endif   else if tag_exist(RPE, 'disk_lambda') and tag_exist(RPE, 'disk_spectrum') then begin
      
; if not defined, only include lines with  A-values above 1e10
      if  tag_exist(RPE, 'min_avalue') then min_avalue=RPE.min_avalue else begin
       IF keyword_set(verbose) THEN $
          print, '% CH_LOAD_ION_RATES: calculating by default RPE for transitions with  A-values > 1e10'
         min_avalue=1e10
      end

      
; select the lines within the input disk spectrum (+/-  0.1 Angstroms).
; Only include observed lines, as presumably the other ones are very weak: 

; *** cut out all lines above the cut in the number of levels:       
      index_lines=where( ion_data.wgfastr.wvl ge min(RPE.disk_lambda)+0.1 and $
                        ion_data.wgfastr.wvl le max(RPE.disk_lambda)-0.1 and $
                         ion_data.wgfastr.aval gt min_avalue and $
                         ion_data.wgfastr.lvl1 lt n_levels and $
                         ion_data.wgfastr.lvl2 le n_levels ,nlines)

      if nlines gt 0 then begin
        IF keyword_set(verbose) THEN print,'% CH_LOAD_ION_RATES: adding RPE to '+trim(nlines)+' transitions:'
        IF keyword_set(verbose) THEN print, ion_data.wgfastr[index_lines].wvl

       
        
         for ii=0, nlines-1 do begin
            
            lambda_0=ion_data.wgfastr[index_lines[ii]].wvl

; first find the IDL indices of the levels:
; Note: the level ordering could be strange.

       ind_lower=where( level_index eq ion_data.wgfastr[index_lines[ii]].lvl1,n_lower)
       if n_lower ne 1 then message, '% CH_LOAD_ION_RATES: ERROR in input lower level number ! '

       ind_upper=where( level_index eq  ion_data.wgfastr[index_lines[ii]].lvl2,n_upper)
       if n_upper ne 1 then message, '% CH_LOAD_ION_RATES: ERROR in input upper level number ! '

; fix these into scalar:       
ind_lower=ind_lower[0]
       ind_upper=ind_upper[0]
              
; Einstein absorption coefficient:
            B_lu= mult[ind_upper]/mult[ind_lower]* $
                  ion_data.wgfastr[index_lines[ii]].aval* (lambda_0*1d-8)^3./(2.*h*c)


    if not tag_exist(RPE, 'a') then begin 

       no_input_a=1
       
; Calculate the scattering factor for the transition
; See Hamilton (1947) and Del Zanna (2024)
    
    j_l=jj[ind_lower] & j_u=jj[ind_upper]
    
    Delta_j = j_u -j_l

    case Delta_j of 
       1: begin
          E_1= 1./ 10 *  (26. *j_l^2+ 37 * j_l+10) / ((j_l+1)*(2*j_l+1)) 
          E_2= 1./ 10 *  ((2*j_l+5)*(j_l+2)) /  ((j_l+1)*(2*j_l+1))             

       end 
       0:begin
          E_1= 1./ 10 *  (12.* j_l^2+12*j_l+1)  / (j_l*(j_l+1))
          E_2= 1./ 10 *  ((2.*j_l-1) * (2*j_l+3)) / (j_l*(j_l+1))
       end 
       -1:begin
          E_1= 1./ 10 *  (26. * j_l^2+ 15*j_l-1) / (j_l*(2*j_l+1))
          E_2= 1./ 10 *  ((2.*j_l-3)*(j_l-1))  /  (j_l*(2*j_l+1))
       end       
             ELSE: BEGIN 
               print, '% CH_LOAD_ION_RATES: ERROR? Delta J of the requested transition is  ',delta_j
               print, '% CH_LOAD_ION_RATES:  setting coefficients a=1, b=0'
               E_1=4./3 & E_2=0.
            END            
       endcase

; add the info to the structure:
    RPE=add_tag(RPE, 3./4*E_1, 'a')
    RPE=add_tag(RPE, 3./4*E_2, 'b')

    if keyword_set(verbose) then print,'% CH_LOAD_ION_RATES, RPE scattering: a='+trim(3./4*E_1)+$
                                       ' b='+trim(3./4*E_2)

 endif else no_input_a=0

    
; Calculate the Doppler factor:
    doppler_factor=doppler_factor_theta(lambda_0, RPE, verbose=verbose)

    if doppler_factor eq -1 then begin 
       error_msg='% CH_LOAD_ION_RATES: ERROR, RPE  model not recognised ! '
   print,error_msg
   return,-1
end 

; Neglect the stimulated emission.

    
; include the PE rate into the matrix:          
          aax[ion_data.wgfastr[index_lines[ii]].lvl1-1, ion_data.wgfastr[index_lines[ii]].lvl2-1] =$
             (lambda_0*1.d-8)*B_lu *h /(4.*!pi)*doppler_factor

; we need to reset the scattering factors for the next line           
          if no_input_a then begin
             RPE=rem_tag(RPE,'a')
             RPE=rem_tag(RPE,'b')
          endif 
          
          
       endfor ; ii list of lines to include 
         
      endif 
      
   endif else message,'Error, input not defined ! '

   
endif  else begin               ; isotropic or Bi-maxwellian options

   error_msg='% CH_LOAD_ION_RATES: ERROR, RPE  model not recognised ! '
   print,error_msg
   return,-1

   
end

    
 endif ; end of PE section
 

;----------------------------------------
; create a ppr array for the proton rates
;---------------------------------------------------

  ppr=MAKE_ARRAY(nt,n_levels,n_levels,/double)
  IF n_tags(prot_struc) NE 0 THEN BEGIN
                                ;
     FOR i=0,n_elements(prot_struc)-1 DO BEGIN
        l1=prot_struc[i].lvl1-1
        l2=prot_struc[i].lvl2-1
        de=ABS(prot_struc[i].de)
        descale_all,t,prot_struc,i,prate
        IF ecmc(l1) LT ecmc(l2) THEN BEGIN
           ppr[*,l1,l2]=prate*pe_ratio*sum_mwl_coeffs
           ppr[*,l2,l1]=prate*pe_ratio*mult[l1]/mult[l2]* $
                        exp(de*13.61/8.617/10.^(-5)/t)*sum_mwl_coeffs
        ENDIF ELSE BEGIN
           ppr[*,l2,l1]=prate*pe_ratio*sum_mwl_coeffs
           ppr[*,l1,l2]=prate*pe_ratio*mult[l2]/mult[l1]* $
                        exp(de*13.61/8.617/10.^(-5)/t)*sum_mwl_coeffs
        ENDELSE
     ENDFOR
  ENDIF

;--------------------------------------
; Create a qq array for electron rates
;---------------------------------------

  qq=MAKE_ARRAY(nt, n_levels,n_levels,/double)
;
  l1=splstr.data.lvl1-1
  l2=splstr.data.lvl2-1


; NOTE: versions older than 9 (with the dielectronic files)
; required an ip_cm in the definition:


IF dielectronic THEN $
  ip_cm=ion_data.ip   ELSE $             ; ioniz. potential in cm^-1 
  ip_cm=0

  kte=(hck*abs(ecmc[l1]-(ecmc[l2]-ip_cm))) # (1d0/temp)
;******************************************


;
; xx and yy contain all factors in the expression for the 
; collisional excitation rate coefficient, 
; except for the upsilon. They can be generated using array operations - the 
; upsilons need a for loop.
;

  xx=dblarr(n_elements(splstr.data),nt)
  yy=dblarr(n_elements(splstr.data),nt)

  ind_pos=where(ecmc[l2] GT ecmc[l1])
  ind_neg=where(ecmc[l2] LT ecmc[l1])

; GDZ: fixed the correct scaling factor, instead of   8.63d-6
  
  IF ind_neg[0] NE -1 THEN BEGIN
     xx[ind_neg,*]=(8.62913438d-6/(mult[l1[ind_neg]])#(1./sqrt(t)))
     yy[ind_neg,*]=8.62913438d-6* exp(-kte[ind_neg,*]) * $
                   1./( mult[l2[ind_neg]] # sqrt(t) )
  ENDIF
  IF ind_pos[0] NE -1 THEN BEGIN
     yy[ind_pos,*]=(8.62913438d-6/mult[l2[ind_pos]]) # (1./sqrt(t))
     
     xx[ind_pos,*]=8.62913438d-6* exp(-kte[ind_pos,*]) * $ ; excitation 
                   1./(mult[l1[ind_pos]] # sqrt(t))
  ENDIF


;
; this is the for loop for the upsilons
;
  FOR i=0L,n_elements(splstr.data)-1 DO BEGIN
     
     IF (l1[i] LE n_levels-1) AND (l2[i] LE n_levels-1) THEN BEGIN
        
; descale the scaled Upsilon at the temperature t:
        
        descale_scups,t,splstr,i,ups ; GDZ
        
        
; NOTE: THE CORRECTIONS DUE TO RECOMBINATION ARE NOT INCLUDED HERE
        
        qq[*,l1[i],l2[i]]=xx[i,*]*ups*sum_mwl_coeffs  ; excitation                  
        qq[*,l2[i],l1[i]]=yy[i,*]*ups*sum_mwl_coeffs  ; de-excitation
        

        
     ENDIF 
  ENDFOR

  if keyword_set(verbose) and  n_elements(RPE) gt 0 then begin

     print,'% CH_LOAD_ION_RATES, temp:',arr2str(temp)
      if tag_exist(RPE, 'radiance_ergs')  then $
         print,'% CH_LOAD_ION_RATES, CE rates: '+ $
               string(qq[*, RPE.lower_level-1, RPE.upper_level-1] , format='(e10.2)')

  end
  
;
; Get the total radiative and dielectronic recombination rates.
; Note: this isn't used in the standard, single-ion atomic
; models, but is used in the multi-ion models.
;
;; tot_rr=recomb_rate(ion_data.gname,t,/rad)
;; tot_dr=recomb_rate(ion_data.gname,t,/diel)
;; tot_iz=ioniz_rate(ion_data.gname,t)

;
; For non-Maxwellians (sum_mwl_coeffs set), we have to sum the
; collision rates (electrons and protons) over temperature.
;
  IF sumtst EQ 1 THEN BEGIN
     qq=total(temporary(qq),1)
     ppr=total(temporary(ppr),1)
     ;; tot_rr=total(temporary(tot_rr),1)
     ;; tot_dr=total(temporary(tot_dr),1)
     ;; tot_iz=total(temporary(tot_iz),1)
  ENDIF 

  output={ n_levels: n_levels, $
           aa: aa, $
           aax: aax, $
           ppr:ppr, $
           qq:qq, $
           ;; tot_rr: tot_rr, $
           ;; tot_dr: tot_dr, $
           ;; tot_iz: tot_iz, $
           temp: temp, ion_data:ion_data, mult:mult, $
           sum_mwl_coeffs: sum_mwl_coeffs, $
           sumtst: sumtst}


;
; If sumtst=0 then need to make sure sum_mwl_coeffs is deleted before exiting.
;
  IF sumtst EQ 0 THEN junk=temporary(sum_mwl_coeffs)

  return,output

end 
