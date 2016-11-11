FUNCTION dgamma, x, a, b

l = (b^a)/gamma(a)
m = x^(a-1.)
n = exp(-1.*b*x)
p = l*m*n

o = WHERE(finite(l*m*n, /NAN), n_o)
IF n_o GT 0 THEN p[o] = 0.

;To copy R, return 0 if NAN (*x to ensure same format)
RETURN, p

END
