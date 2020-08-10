pro shuffle_check,ppstfile=ppstfile

;Read PPST

if n_elements(ppstfile) eq 0 then ppstfile=pickfile(filter='PPST*')
  gapbegin=''
  gapend=''
  readcol,ppstfile,time_o,ppt_o,be_o,name_o,tid_o,seg_o,obsnum_o,ra_ppst_o,dec_ppst_o,roll_ppst_o,format='(a,a,a,a,a,a,a,f,f,f)',delimiter='|'

;Select the PPST start lines so I don't do the work twice

begin_lines=where(strtrim(be_o,2) eq 'Begin')
time=time_o(begin_lines)
be=be_o(begin_lines)
name=name_o(begin_lines)
tid=tid_o(begin_lines)
seg=seg_o(begin_lines)
obsnum=obsnum_o(begin_lines)
ra_ppst=ra_ppst_o(begin_lines)
dec_ppst=dec_ppst_o(begin_lines)
roll_ppst=roll_ppst_o(begin_lines)

roll_ppst_orig=roll_ppst

;Add the star tracker field - Spacecraft roll to the roll angles
;There is a 90degrees+17degres offset between spacecraft 
;coord and star tracker coord

roll_ppst=253.0-roll_ppst

;Read the star catalog
readcol,'/home/pagani/XDS/star_track/test_shuffle/fl_cat.txt',ra_cat,dec_cat,mag_cat,id_cat,format='(f,f,f,a)',skip=1,/silent

;DENSITY
radec6deg=intarr(n_elements(ra_ppst))
radecdens=intarr(n_elements(ra_ppst))
border=intarr(n_elements(ra_ppst))

;Check for stars in FoV 

openw,nomelogico,'log_st_density.txt',/get_lun
openw,warning4,'warning4.txt',/get_lun

for incppst=0l,(n_elements(ra_ppst)-1) do begin
    printares=[ra_ppst(incppst),dec_ppst(incppst),roll_ppst_orig(incppst),roll_ppst(incppst)]
;    print,'Target',name(incppst),'   Time:',time(incppst)
;    print,'Pointing: ',printares
    printf,nomelogico,'Target',name(incppst),'   Time:',time(incppst)
    printf,nomelogico,'Pointing:    RA         DEC        ROLL_PPST     MATRIX_ROT_ANGLE'
    printf,nomelogico,'     ',printares
    for inccat=0l,n_elements(ra_cat)-1 do begin
        res = sphdist(ra_cat(inccat), dec_cat(inccat), ra_ppst(incppst), dec_ppst(incppst),/degrees)  ;angles in DEGREES
        if res lt 6.0 then begin  ; if star closer than 6.0 degrees check if it's in the 8x8 camera field of view
            radec6deg(incppst)=radec6deg(incppst)+1
            
                                ;This rotate_roll rotate thing
                                ;clockwise.  It rotates the coordinate
                                ;of the neraby stars of the catalog stars
                                ;around the RA,DEC of the PPST poining
                                ;(the spacecraft pointing) by the
                                ;amound of the roll angle (modified so
                                ;that the roll angle definition and
                                ;this routing rotate by the same
                                ;amount)


            rotate_roll,ra_ppst(incppst),dec_ppst(incppst),ra_cat(inccat),dec_cat(inccat),roll_ppst(incppst)*!dtor,ra_cat_new,dec_cat_new
                                ;Here I do a simple check if there is
                                ;a situation like this: RA_PPST+359.8
                                ;and RA_CAT=0.5, then the two sources
                                ;are very close but H=359!!! not real!!
            if (abs(ra_cat_new-ra_ppst(incppst)) gt 340.) then begin
                if ra_cat_new gt ra_ppst(incppst) then ra_cat_new=ra_cat_new-360. else ra_cat_new=ra_cat_new+360
            endif
   ;         printf,nomelogico,'Catalog star RA,DEC after ROTATION:',ra_cat_new,dec_cat_new
            
            H=(ra_cat_new-ra_ppst(incppst))*cos(dec_ppst(incppst)*!dtor)   ;Star Tracker H,V stars coord
            V=dec_cat_new-dec_ppst(incppst)
            
            if abs(H) le 4.0 and abs(V) le 4.0 then begin
 ;               print,'H= ',H,'  V=',V,'   Mag=',mag_cat(inccat),'    ID:',id_cat(inccat)
                printf,nomelogico,'H= ',H,'   V=',V,'   Mag=',mag_cat(inccat),'    ID:',id_cat(inccat)
                radecdens(incppst)=radecdens(incppst)+1 ; That is if the star is in the field of view
                if  abs(H) gt 3.75 or abs(V) gt 3.75 then border(incppst)=border(incppst)+1
            endif
        endif

    endfor
    if radecdens(incppst) le 2 then printf,warning4,'WARNING!!!  Not enough stars in FoV for ',name(incppst),'   at ',time(incppst)
    if (radecdens(incppst)-border(incppst) le 2 and border(incppst) gt 0.0) then printf,warning4,'WARNING!!!  Stars at border causing problems!!! Not enough stars in FoV for ',name(incppst),'   at ',time(incppst)
    
  ;  print,'Number of stars within pointing 6 degrees / in star tracker Field of View: ',radec6deg(incppst),radecdens(incppst)
    printf,nomelogico,'Number of stars within pointing 6 degrees / in star tracker Field of View:  ',radec6deg(incppst),radecdens(incppst)
endfor

free_lun,nomelogico
free_lun,warning4


end
