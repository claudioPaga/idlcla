;file astro.pro - a collection of useful astronomical routines
;(c) Detlef Koschny, 1997 Sep 14 - GMST, LMST
;                    1997 Sep 15 - add normalize, XYZ2xxx...
;                    1997 Sep 16 - update normalize to any dimension
;                    1998 Jun 25 - add CROSS_PROD, SCAL_PROD
;                    1998 Jun 28 - add R_Earth, RA_hms (not yet working)
;                    1998 Jun 30 - add Cel_Dist - not yet thoroughly tested!!
;                    1998 Jul 02 - add trunc, round, needed for RAhms
;--------------------------------------------------------------------------

function normalize, vector
;+
; NAME:
;	normalize
;
; PURPOSE:
;       normalizes a vector
;       needed e.g. in XYZ2Cel
;
; CATEGORY:
;	Vector analysis
;
; CALLING SEQUENCE:
;	normalized_vector = normalize (vector)
;
; INPUTS:
;       vector:  A column vector with more than one element.
;
; OPTIONAL INPUT PARAMETERS:
;	None.
;
; OUTPUTS:
;       The normalized vector, i.e. its length is 1.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Hopefully none.
;
; RESTRICTIONS:
;       Returns an error if input is not a vector.
;	Gives a division by zero error if all coordinates are 0.
;
; MODIFICATION HISTORY:
;       1997 Sep 16, update to any dimension
;	1997 Sep 15, first version.
;       (c) Detlef Koschny
;-
  result = 0.
  info = size (vector)
  if info (0) NE 1 then print, 'NORMALIZE: not a one-dimensional vector!!!' $
  else if norm (vector) EQ 0 then print, 'NORMALIZE: vector length is 0!!!' $
  else begin
     l_of_vector = 0.
     for i = 0, info (1)-1 do l_of_vector = l_of_vector + vector (i)^2
     result = vector / sqrt (l_of_vector)
  endelse
return, result
end ;function normalize

;--------------------------------------------------------------------------

function cross_prod, vector1, vector2
;+
; NAME:
;	cross_prod
;
; PURPOSE:
;       take the cross product of the 3-dim. vectors, used for finding
;       the normal to these two vectors
;
; CATEGORY:
;	Vector analysis
;
; CALLING SEQUENCE:
;	normal_vector = cross_prod (vector1, vector2)
;
; INPUTS:
;       vector1:  A column vector with three elements.
;       vector2:  Another column vector with three elements.
;
; OPTIONAL INPUT PARAMETERS:
;	None.
;
; OUTPUTS:
;       The cross product of the two vectors.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Hopefully none.
;
; RESTRICTIONS:
;       Returns an error if input vectors do not have three elements.
;
; MODIFICATION HISTORY:
;       1998 Jun 25, first version.
;       (c) Detlef Koschny
;-
  result = [1.0,2.0,3.0]          ; dummy vector
  info1 = size (vector1)
  info2 = size (vector2)

  if info1 (0) NE 1 then $
     print, '% CROSS_PROD: vector 1 not a one-dimensional vector!!!' $
  else if info2 (0) NE 1 then $
     print, '% CROSS_PROD: vector 2 not a one-dimensional vector!!!' $
  else if info1 (1) NE 3 then $
     print, '% CROSS_PROD: vector 1 not a three-element vector!!!' $
  else if info2 (1) NE 3 then $
     print, '% CROSS_PROD: vector 2 not a three-element vector!!!' $
  else begin
     result (0) = vector1 (1) * vector2 (2) - vector1 (2) * vector2 (1)
     result (1) = - (vector1 (0) * vector2 (2) - vector1 (2) * vector2 (0))
     result (2) = vector1 (0) * vector2 (1) - vector1 (1) * vector2 (0)
  endelse
return, result
end ;function cross_prod

;--------------------------------------------------------------------------

function scal_prod, vector1, vector2
;+
; NAME:
;	scal_prod
;
; PURPOSE:
;       take the scalar product of the 3-dim. vectors, used for finding
;       the angle between two vectors
;
; CATEGORY:
;	Vector analysis
;
; CALLING SEQUENCE:
;	scalar_product = scal_prod (vector1, vector2)
;
; INPUTS:
;       vector1:  A column vector with any number of elements.
;       vector2:  Another column vector with the same number of elements.
;
; OPTIONAL INPUT PARAMETERS:
;	None.
;
; OUTPUTS:
;       The scalar product of the two vectors, defined as the sum of the
;       individual elements. Type is real.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Hopefully none.
;
; RESTRICTIONS:
;       Returns an error if input vectors do not have equal no. of elements.
;
; MODIFICATION HISTORY:
;       1998 Jun 25, first version.
;       (c) Detlef Koschny
;-
  result = 0.0          ; dummy value
  info1 = size (vector1)
  info2 = size (vector2)

  if info1 (0) NE 1 then $
     print, '% SCAL_PROD: vector 1 not a one-dimensional vector!!!' $
  else if info2 (0) NE 1 then $
     print, '% SCAL_PROD: vector 2 not a one-dimensional vector!!!' $
  else if info2 (1) NE info1 (1) then $
     print, '% SCAL_PROD: vector 1 and 2 do not have same no. of elements!!!' $
  else begin
     for i = 0, info1 (1)-1 do result = result + vector1 (i) * vector2 (i)
  endelse
return, result
end ;function scal_prod
;--------------------------------------------------------------------------

function vec_length, vector
;+
; NAME:
;	vec_length
;
; PURPOSE:
;       find the length of a vector
;
; CATEGORY:
;	Vector analysis
;
; CALLING SEQUENCE:
;	l_v = vec_length (vector)
;
; INPUTS:
;       vector:  A column vector with any number of elements.
;
; OPTIONAL INPUT PARAMETERS:
;	None.
;
; OUTPUTS:
;       The length of the vector, i.e. the square root of the sum of the
;       squares of all elements. Type is real.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Hopefully none.
;
; RESTRICTIONS:
;       none.
;
; MODIFICATION HISTORY:
;       1998 Jun 30, first version.
;       (c) Detlef Koschny
;-
  result = 0.0          ; dummy value
  info = size (vector)
  if info (0) NE 1 then print, '% VEC_LENGTH: not a one-dimensional vector!!!' $
  else if norm (vector) EQ 0 then print, '% VEC_LENGTH: vector length is 0!!!' $
  else begin
     l_of_vector = 0.
     for i = 0, info (1)-1 do l_of_vector = l_of_vector + vector (i)^2
     result = sqrt (l_of_vector)
  endelse
return, result
end ;function vec_length

;--------------------------------------------------------------------------

function GMST, JD0, UT
;+
; NAME:
;	GMST
;
; PURPOSE:
;       finds the Greenwich mean siderial time from the Julian date at 0h UT
;       and for the time of the day
;
; CATEGORY:
;	Astronomy
;
; CALLING SEQUENCE:
;	Result = GMST (JD0, UT)
;
; INPUTS:
;       JD0: Julian date at 0h UT (use JULDAY of standard userlib)
;       UT:  Universal time in decimal hours
;
; OPTIONAL INPUT PARAMETERS:
;	None.
;
; OUTPUTS:
;       GMST returns the Greenwich Mean Siderial Time for given UT, on
;       the meridian (i.e. longitude 0), in decimal hours
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Hopefully none.
;
; RESTRICTIONS:
;	Probably doesn't work before 1682 (calendar reformation).
;
; MODIFICATION HISTORY:
;	1997 Sep 14, first version. Adapted from Turbo Pascal
;       (c) Detlef Koschny
;-

  if (UT < 0) or (UT > 23.9999) then err_astro = 2
  help = 6.656306D + 0.0657098242D * (JD0 - 2445700.5D) + 1.002739093D * UT;

  ;range of GMST is 0..23.9999
  while help LT 0 do help = help + 24.
  while help GE 24 do help = help - 24.
  GMST = help

return, GMST
end ;function GMST

;--------------------------------------------------------------------------

function LMST, JD0, UT, lon

;+
; NAME:
;	LMST
;
; PURPOSE:
;       calculates local mean siderial time for given date, time, and longitude

; CATEGORY:
;	Astronomy
;
; CALLING SEQUENCE:
;	result = LMST (JD0, UT, lon)
;
; INPUTS:
;       JD0:     Julian date
;       UT:      Universal time in decimal hours (0 .. 23.9999)
;	lon:     longitude of point on Earth in rad, Eastern longitudes are
;                counted positive (-double (!dtor)*180 .. +double (!dtor)*180)
;
; OPTIONAL INPUT PARAMETERS:
;	None.
;
; OUTPUTS:
;       LMST returns the local sideral time, i.e. the right ascension of
;       the sky in the meridian of the current longitude.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Hopefully none.
;
; RESTRICTIONS:
;	Probably doesn't work before 1682 (calendar reformation).
;       Very important NOTE: I use: long > 0 is East!!!
;       This is unlike e.g. Montenbruck.
;
; MODIFICATION HISTORY:
;	1997 Sep 14, first version, adapted from Turbo Pascal.
;       (c) Detlef Koschny
;-

  help = GMST (JD0, UT) + double (!radeg) * lon/15.
  ;range of LMST is 0..23.9999
  while help LT 0 do help = help + 24.
  while help GE 24 do help = help - 24.
  LMST = help

return, LMST
end ;function LMST

;--------------------------------------------------------------------------

function R_Earth, lat
;+
; NAME:
;	R_Earth
;
; PURPOSE:
;	Calculate earth's radius as a function of latitude
;
; CATEGORY:
;	Astronomy
;
; CALLING SEQUENCE:
;	Result = R_Earth (lat)
;
; INPUTS:
;       lat:     latitude of point on Earth in rad, North is positive, South
;                is negative (-!dtor*90 .. +!dtor*90)
;
; OPTIONAL INPUT PARAMETERS:
;	None.
;
; OUTPUTS:
;       R_Earth returns the radius of the earth in meter for a given
;       latitude.
;
; COMMON BLOCKS:
;       none.
;
; SIDE EFFECTS:
;	Hopefully none.
;
; RESTRICTIONS:
;	none.
;
; MODIFICATION HISTORY:
;	1998 Jun 28, first version. (c) Detlef Koschny
;-

; define variables, set constants
r_Equ = 6378388.0D                        ;equatorial radius of Earth in m
r_Pol = 6356911.9D                        ;pol radius of Earth in m

; calculate r from geometry of ellipse, see paper doc no. I
result = sqrt ((r_Equ^2 - r_Pol^2) / ((r_Equ/r_Pol * TAN (lat))^2+1) + r_Pol^2)

return, result
end

;--------------------------------------------------------------------------

function Geo2XYZ, lon, lat, h, year, month, day, dectime
;+
; NAME:
;	Geo2XYZ
;
; PURPOSE:
;	Determine the x,y,z coordinates in a non-rotating reference frame
;	of a point on the Earth, given time and location.
;
; CATEGORY:
;	Astronomy
;
; CALLING SEQUENCE:
;	Result = Geo2XYZ (lon, lat, h, year, month, day, dectime)
;
; INPUTS:
;	lon:     longitude of point on Earth in rad, Eastern longitudes are
;                counted positive (-!dtor*180 .. +!dtor*180)
;       lat:     latitude of point on Earth in rad, North is positive, South
;                is negative (-!dtor*90 .. +!dtor*90)
;       h:       altitude above mean sea level in meter
;       year:    number of the desired year, four digits (1682..9999)
;       month:   number of the desired month (1 .. 12)
;       day:     number of the day of the month (1 .. 31)
;       dectime: time in decimal hours (0.0 .. 23.99)
;
; OPTIONAL INPUT PARAMETERS:
;	None.
;
; OUTPUTS:
;	Geo2XYZ returns a vector with three elements, giving x,y,z-coordinates
;	of the specified point in meters. The coordinate system is right handed,
;       +x is the direction to the vernal equinox, +z is North
;
; COMMON BLOCKS:
;       Needs GMST and LMST from file astro.pro.
;
; SIDE EFFECTS:
;	Hopefully none.
;
; RESTRICTIONS:
;	Probably doesn't work before 1682 (calendar reformation).
;
; MODIFICATION HISTORY:
;	1997 Sep 14, first version. (c) Detlef Koschny
;-

; define variables, set constants
GeoCoords = DBLARR (3)                 ;x, y, z coords will be stored here
r_Equ = 6378388.0D                        ;equatorial radius of Earth in m
r_Pol = 6356911.9D                        ;pol radius of Earth in m

; calculate r from geometry of ellipse, see paper doc no. I
r = sqrt ((r_Equ^2 - r_Pol^2) / ((r_Equ/r_Pol * TAN (lat))^2+1) + r_Pol^2) + h

; Now convert to x,y,z - paper doc no. II

; determine sid time
JD0 = JULDAY (month, day, year)
sidtime = LMST (JD0, dectime, lon)

; convert to an angle in radians
sidtime = double (!dtor) * sidtime*15

; convert spherical to orthogonal
GeoCoords (0) = r * cos (lat) * cos (sidtime)
GeoCoords (1) = r * cos (lat) * sin (sidtime)
GeoCoords (2) = r * sin (lat)

return, GeoCoords
end

;--------------------------------------------------------------------------

function Cel2XYZ, RA, Dec
;+
; NAME:
;	Cel2XYZ
;
; PURPOSE:
;	Determine the x,y,z coordinates in a universal reference frame
;	of a direction to a given RA and Declination.
;
; CATEGORY:
;	Astronomy
;
; CALLING SEQUENCE:
;	result = Cel2XYZ (RA, Dec)
;
; INPUTS:
;       RA:  right ascension in arcus (0 * !dtor .. 360 * !dtor)
;       Dec: declination in arcus (-90 * !dtor .. +90 * !dtor)
;
; OPTIONAL INPUT PARAMETERS:
;	None.
;
; OUTPUTS:
;       Cel2XYZ returns a vector with three elements, giving the
;       direction of the celestial coordinates entered. The length of
;       the vector is 1. x-direction is to vernal equinox, z-direction
;       is due North.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Hopefully none.
;
; RESTRICTIONS:
;	Don't know.
;       If Dec > 90 deg, a warning is printed on the screen.
;
; MODIFICATION HISTORY:
;	1997 Sep 14, first version. (c) Detlef Koschny
;-

CelCoords = DBLARR (3)

r = 1.0D

if double (!radeg) * Dec GT 90. then print, 'Cel2XYZ: Warning - Dec > 90 deg!'
CelCoords (0) = r * cos (Dec) * cos (RA)
CelCoords (1) = r * cos (Dec) * sin (RA)
CelCoords (2) = r * sin (Dec)

return, CelCoords
end  ;function Cel2XYZ

;--------------------------------------------------------------------------

function XYZ2Cel, OrthoCoords
;+
; NAME:
;	XYZ2Cel
;
; PURPOSE:
;       determines RA and Dec from orthogonal coordinate vector
;
; CATEGORY:
;	Astronomy, coordinate transformation
;
; CALLING SEQUENCE:
;	CelCoords = XYZ2Cel (OrthoCoords)
;
; INPUTS:
;       OrthoCoords:  A vector with three elements. Gives x, y, z values
;       of a viewing direction.
;
; OPTIONAL INPUT PARAMETERS:
;	Needs function 'normalize'.
;
; OUTPUTS:
;       The right ascension (CelCoords (0)) and declination (CelCoords (1))
;       of the direction of the vector in radians.
;       x-direction is to vernal equinox, z-direction is due North.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Hopefully none.
;
; RESTRICTIONS:
;	Gives a division by zero error in function normalize if
;       all coordinates are 0.
;
; MODIFICATION HISTORY:
;	1997 Sep 15, first version.
;       (c) Detlef Koschny
;-

; generate CelCoords as 2-element vector, double precision
CelCoords = DBLARR (2)

; norm OrthoCoords first
OrthoCoords = normalize (OrthoCoords)

; and now the conversion. The ATAN returns negative angles, we correct this
CelCoords (0) = ATAN (OrthoCoords (1), OrthoCoords (0))   ;RA
while (CelCoords (0) LT 0) do CelCoords (0) = CelCoords (0) + 2.*3.1415926535D
CelCoords (1) = ASIN (OrthoCoords (2))                    ;Dec

return, CelCoords
end    ;function XYZ2Cel

;--------------------------------------------------------------------------

function Cel_Dist, RA1, Dec1, RA2, Dec2
;+
; NAME:
;   Cel_Dist
;
; PURPOSE:
;   finds angular distance in radians between two points in the sky
;
; CATEGORY:
;   astro
;
; CALLING SEQUENCE:
;   distance = Cel_Dist (RA1, Dec1, RA2, Dec2)
;
; INPUTS:
;   RA1:  Right ascension of first point in radians
;   Dec1: Declination of first point in radians
;   RA2:  Right ascension of second point in radians
;   Dec2: Declination of second point in radians
;
; OPTIONAL INPUT PARAMETERS:
;   None.
;
; OUTPUTS:
;   Angular distance between point one and two in radians, type real.
;   The following procedure is used: The angle is the arccos ((A*B)/(|A||B|))
;   where A and B are the vectors in x,y,z of the viewing directions,
;   A*B the scalar product of these vectors, and |A| and |B| the length
;   of the vectors.
;
; COMMON BLOCKS:
;   None.
;
; SIDE EFFECTS:
;   Hopefully none.
;
; RESTRICTIONS:
;   Never know what the arcos does...
;
; MODIFICATION HISTORY:
;   1998 Jun 30, first version.
;   (c) Detlef Koschny
;-

  result = 0.0  ;define real number
  vec1 = fltarr (3)  ;define vectors with three elements
  vec2 = fltarr (3)

  vec1 = Cel2XYZ (RA1, Dec1)
  vec2 = Cel2XYZ (RA2, Dec2)
  result = acos (  (scal_prod (vec1, vec2)) $
                 / (vec_length (vec1) * vec_length (vec2)))
  return, result
end ;function Cel_Dist

;---------------------------------------------------------------------------

function RAh, RArad
;+
; NAME:
;	RAh
;
; PURPOSE:
;       extract hour value from an angle in radians
;       needed for the conversion of right ascension in h/m/s
;
; CATEGORY:
;	astronomy
;
; CALLING SEQUENCE:
;	RAh (RArad)
;
; INPUTS:
;       RArad: The right ascension in radians
;
; OPTIONAL INPUT PARAMETERS:
;	None.
;
; OUTPUTS:
;       The hours corresponding to the input angle, type is integer.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Hopefully none.
;
; RESTRICTIONS:
;
;
; MODIFICATION HISTORY:
;       1998 Jul 02, separate from RAhms.
;       1998 Jun 28, first version.
;       (c) Detlef Koschny
;-

  RA = RArad * double (!radeg) / 15 ;angle in decimal hours
  result = fix (RA)
return, result
end ;RAh

function RAm, RArad
  RA = RArad * double (!radeg) / 15 ;angle in decimal hours
  result = fix ((RA - fix (RA))*60)
  return, result
end ;RAm


function RAs, RArad
  RA = RArad * double (!radeg) / 15. ;angle in decimal hours
  result = float ((RA - RAh (RArad) - RAm (RArad)/60.)*3600)
return, result
end 
