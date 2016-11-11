FUNCTION rmf, n, mmin, mstar, phi1, a1, phi2, a2, seed=seed

log_m = mmin + 0.01*dindgen(600)
m = 10d^(mmin + 0.01*dindgen(600))
mstar = 10d^mstar

a = exp(-1.*m/mstar)
b = phi1*((m/mstar)^a1)
c = phi2*((m/mstar)^a2)

mf = m*a*(b + c)

log_rm = randomf(seed, n, log_m, mf)
rm = 10d^log_rm

RETURN, rm

END
