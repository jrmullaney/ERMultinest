FUNCTION comb_gam, x, plind, log_k, log_bmin, alpha

ytot = 0.*x

b = 10d^(log_bmin+0.3*findgen(20))
norm = b^(-1*(plind+1))
k = 10d^log_k

plot, x, ytot, /xlog, /ylog, xra=minmax(x), yra=[1e-10,1e4]

FOR i=0, n_elements(b)-1 DO BEGIN
   print, b[i]
   y = k*norm[i]*dgamma(x,alpha+1,b[i])
   
   ytot = ytot+y


   oplot, x, y
ENDFOR   

oplot, x, ytot

RETURN, ytot

END
