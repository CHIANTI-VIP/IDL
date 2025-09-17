;+
;
; PROJECT:  CHIANTI
;
;       This program was developed as part of CHIANTI-VIP. 
;       CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;       mantained by Giulio Del Zanna, to develop additional features and
;       provide them to the astrophysics community.
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the
;       University of Cambridge, Goddard Space Flight Center, and University of Michigan.
;
;
; NAME: CALC_PI_RATES
;
;
;
; PURPOSE:
; 
;  Calculate photo-ionization rates from PI cross-section files that have initial and
;  final level indexing, PI thresholds and cross-sections in cm^-2. The files 
;  have the names 'c_1_pi.lr', for example, for C I.
;
;   * Warning *  -  note cross-section from read_pi_cross has structure [neng, ntrans].
;                   The interpolated cs is changed here to the format
;                   [ntrans, neng].

; CATEGORY:
;
;       CHIANTI; ionization; recombination, ion fractions.
;
;
; CALLING SEQUENCE:
;
;	data=calc_pi_rate( filename, wvl, radiance)
;
;
; INPUTS:
;           cross_pi    -   a structure with the level-resolved PI cross-sections
;                           obtained e.g. with cross_pi = read_pi_cs(filename, ryd=ryd)
;
;           wvl         -   incident photon wavelength in A
;
;           radiance    -   radiance in  photons cm-2 s-1 sr-1 A-1
;
;
; OUTPUT:
;           rate_pi     -   PI rate coefficient in s-1
;
; KEYWORDS:
;
;           Ryd
;           input energies are in Rydbergs instead of eV. 
;
;           all
;           to return all the variables in the output
;
;           dilution
;
; MODIFICATION HISTORY:
;
;    v.1, 21-Aug 2025   Giulio Del Zanna (GDZ) and Roger Dufresne (RPD) 
;               DAMTP, University of Cambridge
;
; VERSION    v.1
;
;-

function calc_pi_rates, cross_pi, wvl, radiance, dilution=dilution, ryd=ryd, all=all


ntrans = n_elements(cross_pi.lvl1)
neng = n_elements(cross_pi.energies)
nwvl = n_elements(wvl)
level_i = cross_pi.lvl1
level_f = cross_pi.lvl2
final_energy = cross_pi.energies
incident_energy = dblarr(ntrans, neng)
interp_cross = dblarr(ntrans, nwvl)
rate_pi = dblarr(ntrans)

; convert spectral wavelengths to photon energy in eV for interpolating theoretical cross-sections
spectrum_e = 1.e8 / (wvl * 8065.545)

; loop to calculate rate coefficient for each transition
for ii = 0, ntrans-1 do begin
    
    ; calculate the incident photon energy of the theoretical data
    ;for ieng = 0, neng-1 do begin
        incident_energy[ii, *] = cross_pi.pe[ii] + cross_pi.energies[*]
    ;endfor
    
    ; establish which energies of the spectrum are: below the ionization threshold,
    ;       above the threshold but less than lowest theoretical energy, greater than the
    ;       highest theoretical energy, and within the theoretical values
    ind1=where(spectrum_e lt min(cross_pi.pe[ii]), n1)
    ind2=where(spectrum_e ge min(cross_pi.pe[ii]) and spectrum_e lt min(incident_energy[ii, *]), n2)
    ind3=where(spectrum_e gt max(incident_energy[ii, *]), n3)
    ind4=where(spectrum_e ge min(incident_energy[ii, *]) and spectrum_e le max(incident_energy[ii, *]), n4)
   
    ; output spectral cross-sections which are: zero below threshold, equal to the cross-section
    ;       at the lowest theoretical energy point if between threshold and the lowest energy,
    ;       zero above the highest theoretical energy,
    ;       and interpolated from the theoretical data if the energies lie within the theoretical
    if n1 gt 0 then interp_cross[ii, ind1] = 0.0
    if n2 gt 0 then interp_cross[ii, ind2] = cross_pi.cs[0, ii]
    ;if n3 gt 0 then interp_cross[ii, ind3] = cross_pi.cs[-1, ii]
    if n3 gt 0 then interp_cross[ii, ind3] = 0.0
    if n4 gt 1 then interp_cross[ii, ind4] = $
        10.^(interpol(alog10(reform(cross_pi.cs[*, ii])), reform(incident_energy[ii, *]), spectrum_e[ind4], /spline))

; just make sure there are no zeroes:
    good = where(interp_cross[ii, *] gt 0.0)

; calculate the PI rate coefficient from the interpolated cross-sections    
    rate_pi[ii] = 4. * !pi * int_tabulated(wvl[good], interp_cross[ii, good] * radiance[good], /sort)

;    rate_pi[ii] = 4. * !pi * int_tabulated(wvl, interp_cross[ii, *] * radiance)

endfor

if n_elements(dilution) eq 1 then rate_pi=rate_pi*dilution

if keyword_set(all) then $
return, {level_i:level_i, level_f:level_f, rate_pi:rate_pi, $
         spectrum_e:spectrum_e, interp_cross:interp_cross, $
                    incident_energy:incident_energy, final_energy:final_energy} else $
return, {level_i:level_i, level_f:level_f, rate_pi:rate_pi}



end
