pro prova
readcol,'/home/pagani/table.txt',a,b,format='(i,i)'
d=intarr(2,4)
d(0,*)=a
d(1,*)=b
image=0
mkhdr,hdr,image
sxaddpar,hdr,'Comment',''
mwrfits,d,'/home/pagani/prova.fits',hdr,/create
end

