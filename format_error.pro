pro format_error,pname,p,perror
perrorlow=dblarr(n_elements(p)) 
perrorhigh=dblarr(n_elements(p)) 
for k=0,n_elements(p)-1 do perrorlow(k)=perror(0,k)
for k=0,n_elements(p)-1 do perrorhigh(k)=perror(1,k)
pname[2]='Br_T'
pst=string(p)
plowst=string(perrorlow)
phighst=string(perrorhigh)
param_arr=[transpose(pname),transpose(pst),transpose(plowst),transpose(phighst)]
print,param_arr                 ;Print best fit param 
norm_plateau_end=p[0]*p[2]^(p[1])   ;I calculate the normalization at t=-T_plateau_end
print,'Norm_plateau_end:',norm_plateau_end
print,'Post_plateau_PI,Plateau_end,Norm_plateau_end'
print,pst(3),pst(2),norm_plateau_end,format='(d4.2,d12.3,d12.3)'
checklist=findfile('lc_priorfit_noflares.txt')
if checklist eq '' then spawn,'cp lc_newout_noflares.txt lc_priorfit_noflares.txt'
end
