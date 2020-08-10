pro plaw2_cla,lcname,pindex,offx

readcol,lcname,x,y,e,format='(d,d,d)'
w=1.0/y
x=x-offx
a0guess=y[1]/(x[1])^(pindex)
a=[a0guess,pindex]
fit=curvefit(x,y,w,a,paramerr,function_name='plaw2',yerr=yerr)
print,'a0: (should be = cr when t=1), tells how bright it was at trigger ',a[0]
print,'Slope: ',a[1]
print,'Parameter errors: ',paramerr

plot,x,y,psym=4
opx=(indgen(max(x)*10))/10.0
opy=a[0]*opx^(a[1])
oplot,opx,opy
end
