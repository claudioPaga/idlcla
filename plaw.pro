PRO plaw, X, A, F, pder
   F = A[0]*X^A[1]+A[2]
   IF N_PARAMS() GE 4 then pder=[[X^A[1]] ,[A[0]*X^A[1]*alog(X)], [replicate(1.0,n_elements(X))]]
END
