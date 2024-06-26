;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;
; NAME:
;	READ_DEM
;
; PURPOSE:
;
;	to read values the differential emission measure as a function 
;       of temperature
;
; CATEGORY:
;
;	science.
;
; CALLING SEQUENCE:
;
;       READ_DEM, File, T, Dem, Ref
;
;
; INPUTS:
;
;	File:  the name of the file containing the DEM values, usually in
;               !xuvtop/dem/*.dem	
;
;
; OUTPUTS:
;
;	T:  Log10 values of temperature (K)
;       Dem:  Log10 values of the differential emission measure
;       Ref:  the reference to the DEM values in the scientific literature
;
;
; OPTIONAL OUTPUTS:
;
;	Describe optional outputs here.  If the routine doesn't have any, 
;	just delete this section.
;
;
;
; EXAMPLE:
;
;             > read_dem,!xuvtop+'/dem/active_region.dem',t,dem,ref
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;       V.   3, 21-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;       Ver.4, 20-Nov-2002, Peter Young
;           DEM values of -1 caused errors, so this has been
;           corrected.
;       Ver.5, 11-Dec-2020, Peter Young
;           Made dem a double-precision array; reduced n_params
;           check. 
;
; VERSION     :   5, 11-Dec-2020
;
;-
pro read_dem,filename,t,dem,ref
;
;
;
if n_params() lt 3 then begin
   print,' type> read_dem,filename,temp,dem,ref'
   print,'         or '
   print,' type> read_dem,''',''',temp,dem,ref'
   return
endif
;
;
;
if strtrim(filename,2) eq '' then begin
;
    dir=concat_dir(!xuvtop,'dem')
    filename=dialog_pickfile(path=dir,filter='*.dem',title='Select DEM File')
    print,' selected:  ',filename
endif
;
;
openr,lu,filename,/get_lun
;
t=0.
dem=0d0
ref=''
;

string1=' '
chck1=0

WHILE eof(lu) NE 1 DO BEGIN
  readf,lu,string1
  IF strtrim(string1,2) EQ '-1' THEN BEGIN
    chck1=chck1+1
  ENDIF ELSE BEGIN
    CASE chck1 OF
      0: BEGIN
        reads,string1,a,b
        t=[t,a]
        dem=[dem,b]
      END
      1: ref=[ref,string1]
      2: 
    ENDCASE
  ENDELSE
ENDWHILE

t=t[1:*]
dem=dem[1:*]
ref=ref[1:*]

free_lun,lu

END
