
PRO READ_BEST_ENERGS, FILE, ENERGS, WEIGHTS, INDEX, $
                      ENCTH, ENRTH, enryd=enryd

;+
; NAME
;
;       READ_BEST_ENERGS
;
; PROJECTS
;
;       CHIANTI
;
; EXPLANATION
;
;	We often want to extract from the .elvlc file the best energies 
;	for the levels. This routine does this, returning an energy vector, 
;	with the ordering consistent with the level ordering in the .elvlc 
;	file.
;
;       It looks for the 3rd energy column which is present in some of 
;       the files, and uses this if needed. If the 3rd column exists,
;       but contains a zero, then this value is not used.
;
; INPUTS
;
;   FILE    Name of .elvlc file.
;
; OUTPUTS
;
;   ENERGS  The best energies in cm^-1 for the levels in the file.
;
;   WEIGHTS Statistical weights for the levels.
;
;   INDEX   This contains either a 0, 1 or 2 for each of the levels. The 
;           1 signifying the level energy is theoretical, and that it 
;           comes from the 3rd energy column. 2 signifies that the 
;           energy comes from the 2nd energy column.
;
;   ENCTH   Theoretical energies (from 2nd col.) in cm^-1
;
;   ENRTH   Theoretical energies (from 2nd col.) in Ryd
;
; OPTIONAL OUTPUTS
;
;   ENRYD   The best energy in Rydbergs.
;
; HISTORY
;
;   Ver.1, 16-Mar-01, Peter Young
;   Ver.2, 19-Nov-08, Peter Young
;      made small modifications.
;   Ver.3, 5-Sep-12, Peter Young
;      made sure all energies are double precision.
;-


;---------------------------------------------------------[]
; Copied the following from read_elvlc.pro
;
openr,lue,file,/get_lun           ; open the .elvlc file
;
l1a=intarr(1000)
confa=intarr(1000)
desiga=strarr(1000)
ssa=intarr(1000)
lla=intarr(1000)
spda=strarr(1000)
jja=fltarr(1000)
multa=intarr(1000)
ecma=dblarr(1000)
eryda=dblarr(1000)
ecmzha=dblarr(1000)
erydzha=dblarr(1000)
excolcm=dblarr(1000)
excolry=dblarr(1000)
;
;
string1=' '
l1=1  & desig=' ' & ll=1. & conf=1 & ss=1. & jj=1. & mult=1 & ecm=1.d
spd=' ' & eryd=1.d & ecmzh=1.d  & erydzh=1.d & excol=0d & excolr=0d
;
;

;;;
;;; I want to be able to read the 3rd energy column in the .elvlc file.
;;; To do this I check to see how long `string1' is. The standard length is 99 
;;; characters, so I'll make a condition that checks to see if the string is
;;; longer than this.
;;;
;;; Have changed this to 110, as Ken sometimes adds a label to the levels 
;;; (e.g., Mg IX)

WHILE (STRPOS(string1,'-1') LT 0) OR (STRPOS(string1,'-1') GT 10) DO BEGIN
  READF,lue,string1
  IF (STRPOS(string1,'-1') lt 0) or (STRPOS(string1,'-1') gt 10) THEN BEGIN
    IF STRLEN(string1) gt 110 THEN BEGIN
      READS,string1,l1,conf,desig,ss,ll,spd,jj,mult,ecm,eryd,ecmzh, $
        erydzh,excol, $
        format='(i3,i6,a15,2i3,a2,f5.1,i3,2(f15.3,f15.6),f15.3)'
    ENDIF ELSE BEGIN
      READS,string1,l1,conf,desig,ss,ll,spd,jj,mult,ecm,eryd,ecmzh, $
        erydzh,                                                     $
        format='(i3,i6,a15,2i3,a2,f5.1,i3,3(f15.3,f15.6))'
      excol=0d0
    ENDELSE
;  
    l=l1-1
    l1a(l)=l1
    confa(l)=conf
    desiga(l)=desig
    ssa(l)=ss
    lla(l)=ll
    spda(l)=spd
    jja(l)=jj
    multa(l)=float(mult)
    ecma(l)=ecm
    eryda(l)=eryd
    ecmzha(l)=ecmzh
    erydzha(l)=erydzh
    excolcm(l)=excol
  ENDIF
ENDWHILE

n=N_ELEMENTS(l1a)
energs=DBLARR(n)
enryd=DBLARR(n)
weights=float(multa)
index=intarr(n)
;
FOR i=1,n-1 DO BEGIN           ; note 1st level assumed to have zero energy
  IF ecma(i) NE 0 THEN BEGIN
    energs(i)=ecma(i)       ; best energy is observed energy
    enryd[i]=eryda[i]
  ENDIF ELSE BEGIN
    IF excolcm(i) NE 0 THEN BEGIN
      energs(i)=excolcm(i)  ; best energy is in 3rd column
      enryd[i]=excolcm[i]/109737.32d0
    ENDIF ELSE BEGIN
      energs(i)=ecmzha(i)   ; best energy is in 2nd column
      enryd[i]=erydzha[i]
      index[i]=index[i]+1
    ENDELSE
    index[i]=index[i]+1
  ENDELSE
  IF energs(i) EQ 0. THEN energs(i)=-1.
ENDFOR

ind=WHERE(energs EQ -1.)
energs=energs(0:ind(0)-1)
weights=weights(0:ind(0)-1)
index=index[0:ind[0]-1]
encth=ecmzha[0:ind[0]-1]
enrth=erydzha[0:ind[0]-1]

END
                    
