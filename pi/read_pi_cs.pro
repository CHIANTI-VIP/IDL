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
; NAME: READ_PI_CS
;
;
; PURPOSE:
;       Read the Photoionisation cross sections stored in a file.
;      The file has in the first line ' Ejected electron energies (eV)',
;       then 
;       Initial CHIANTI level, final level (next ion stage),  photon ionisation energy (eV),
;       photoionisation cross-section (cm^2), and description of the initial and final states.
;
;
;  * Warning *  -  note cross-section has structure [neng, ntrans] .
;
; CALLING SEQUENCE:
;
;	st=read_pi_cs('c_3_pi.lr')
;
; INPUTS:
;        file name
; OUTPUT:
;           a structure with
;           energies:  Ejected electron energies (eV) 
;           lvl1: Initial CHIANTI level
;           lvl2: final  CHIANTI level in the next ion charge state
;           pe:   photon ionisation energy (eV)
;           cs:   photoionisation cross-section (cm^2)
;           extra: description of the initial and final states.
;           comments: references
;
; KEYWORDS:
;
;           Ryd
;           input energies are in Rydbergs instead of eV. 
;
;
; MODIFICATION HISTORY:
;
;    v.1, 21-Aug 2025   Giulio Del Zanna (GDZ) and Roger Dufresne (RPD) 
;               DAMTP, University of Cambridge
;
; VERSION    v.1
;
;-

  
function read_pi_cs, file, ryd=ryd


result=query_ascii(file,info)
nlines=info.lines
file_string=strarr(nlines)
openr,lin, file,/get_lun
readf,lin,file_string
free_lun,lin

; read energies - first line 
pp=str_sep(file_string[0],' ')
good=where(trim(pp) ne '')

energies= double(pp[good])
if keyword_set(ryd) then energies=energies*13.6056923

neng=n_elements(energies)

; find the start of the comments:  '-1' within the first 5 characters.
ind=where(trim(strmid(file_string,0,5)) eq '-1')
ind2=ind[1]
ind=ind[0]

; number of transitions:
ntrans=ind-1

comments=file_string[ind+1:ind2-1]

a=0 & b=0 & c=0.0 & r= dblarr(neng)

;lvl1= intarr(ntrans)           ; initial level number (ionising ion)
;lvl2= intarr(ntrans)           ;  final level number    
;pe= dblarr(ntrans)  ; ionization potential,  DE in eV 
; cs= dblarr( neng, ntrans) cross-section in cm^2

str={lvl1:0, lvl2:0, pe:0., cs:dblarr(neng) , extra:''}
output=replicate(str, ntrans)

reads, file_string[1:ntrans], output

if keyword_set(ryd) then output.pe=output.pe*13.6056923


return, {energies:energies, lvl1:output.lvl1, lvl2:output.lvl2, pe:output.pe, cs:output.cs,$
         extra:output.extra, comments:comments}


end 

