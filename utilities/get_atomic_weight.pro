;+
; PROJECT:
;
;        CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;        mantained by Giulio Del Zanna, to develop additional features and
;        provide them to the astrophysics community.
;        Contributions via github are welcomed. 
;
; NAME: get_atomic_weight 
;
; PURPOSE:
;
; Returns the average atomic weight or mass (in grams or kilograms) of an element.
;
; RESTRICTIONS:
;            Up to zinc (iz=30)
;
;  VERSION     : V 1,  24-May-2022 Giulio Del Zanna (GDZ) 
;
;
; http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=
;-


FUNCTION  get_atomic_weight, iz

weights = [1.00794, 4.002602, 6.941, 9.0122, 10.811, 12.0107, 14.0067, $
15.9994, 18.9984, 20.179, 22.98977, 24.3050, 26.9815, 28.086, $
30.9738, 32.065, 35.453, 39.948, 39.0983, 40.078, 44.956, $
47.90, 50.9414, 51.996, 54.9380, 55.845,$  ;' Fe !
 58.9332, 58.6934, 63.546, 65.37]


RETURN, weights[iz-1]

END
