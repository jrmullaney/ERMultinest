FUNCTION ms_ssfr, m, z

log_m = alog10(m/1d9)
r = alog10(1.+z)
a = [1.5, 0.3, 2.5]
mstar = [0.5, 0.36]
b = 0. > log_m-mstar[1]-a[2]*r

log_sfr = log_m - mstar[0] + a[0]*r - a[1]*(b^2.)
sfr = 10d^log_sfr

RETURN, sfr/m

END
