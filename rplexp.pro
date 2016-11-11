FUNCTION rplexp, seed, n, a, b, x0

;Generate x array:
log_x1 = 1./b
log_x0 = alog10(x0)
dlog_x = (log_x1-log_x0)/1d3
log_x = log_x0 + dlog_x*findgen(1d4)
x = 10d^log_x

;Calculate normalisation
igam = gamma(1.-a)*(1.-igamma(1.-a, b*x0))
norm = (b^(1.-a))/igam

;Generate random sample
y = norm*(x^(-a))*exp(-1*x/b)
log_r = randomf(seed, n, log_x, y)
r = 10d^log_r

RETURN, r

END
