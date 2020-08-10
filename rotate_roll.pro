pro rotate_roll,cra,cdec,ra,dec,theta,rap,decp

;cra & cdec are the coordinates of the central point to rotate around
;ra & dec are the coordinates of the position to be rotated
;theta is roll or something like (180-roll) in RADIANS
;rap & decp are the outputs

   cost=cos(theta)
   sint=sin(theta)
   vect2,cdec,cra,v
   vect2,dec,ra,v1

M=[[cost+(1-cost)*v[0]^2,(1-cost)*v[0]*v[1]-sint*v[2],(1-cost)*v[0]*v[2]+sint*v[1]],$
   [(1-cost)*v[1]*v[0]+sint*v[2],cost+(1-cost)*v[1]^2,(1-cost)*v[1]*v[2]-sint*v[0]],$
   [(1-cost)*v[2]*v[0]-sint*v[1],(1-cost)*v[2]*v[1]+sint*v[0],cost+(1-cost)*v[2]^2]]

   v1p=M#v1

   elaz2,v1p,decp,rap

   return
end
