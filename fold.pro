FUNCTION fold, time, period, time0
; Return the phase of a time or array of times relative to a given period and time
; All inputs must be in same units.
; time0 is optional.

IF N_ELEMENTS(time0) EQ 0 THEN time0=0

x = (double(time)-double(time0)) / double(period)

RETURN, x-long(x)

END
