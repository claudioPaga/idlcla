pro plaw_cla,lcname,pindex,offx

readcol,lcname,x,y,e,format='(d,d,d)'
w=1.0/y
x=x-offx
a0guess=y[1]/(x[1])^(pindex)
a=[a0guess,pindex,0.1]
fit=mpcurvefit(x,y,w,a,paramerr,function_name='plaw',yerr=yerr)
print,'a0: (should be = cr when t=1), tells how bright it was at trigger ',a[0]
print,'Slope: ',a[1]
print,'a[2], when T big, cr=0, a[2] should be 0: ',a[2]
print,'Parameter errors: ',paramerr

plot,x,y,psym=4
opx=(indgen(max(x)*10))/10.0
opy=a[0]*opx^(a[1])+a[2]
oplot,opx,opy
end
