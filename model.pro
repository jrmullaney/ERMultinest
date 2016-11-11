PRO model

;Get masses:
;Produces a random sample of galaxies
;with masses following the Ilbert+13 mass function:
seed = 1002
m = rmf(1d4, 9.76, 10.88, 1.68, -0.69, 0.77, -1.42, seed=seed)
;o = WHERE(alog10(m) GE 9.76, n_o)
;m = m[o]
mbh = m/500.

;Allocate AGN R_MS values according to
;log-normal R_MS distribution from Mullaney+15:
r_ms = rr_ms(n_elements(m), -0.369, 0.56, seed=seed)

;Calculate MS_sSFRs using Schreiber+15
ms_ssfr = ms_ssfr(m, 2.0) ;In units of yr^-1

;Calculate sSFR of AGN hosts by multiplying
;ms_ssfr by AGN R_MS values:
agn_ssfr = r_ms * ms_ssfr
sfr_agn = agn_ssfr*m

;Assign an alpha according to sSFR:
;alpha = 0.3*1d9*agn_ssfr - 0.65
alpha = 0.2*1d9*agn_ssfr - 0.65

;plothist, alpha, bin=0.01

;Assign Eddington ratio according to alpha:
;redd = dblarr(n_elements(m))
;FOR i=0L, n_elements(m)-1 DO $
;   redd[i] = rplexp(seed, 1, -1.*alpha[i], 1, 1d-4)
redd = rplexp(seed, n_elements(m), 0.2, 1, 1d-4)


;plothist, alog10(redd), bin=0.1, /ylog
;x = 10d^(-4.+0.01*findgen(500))
;oplot, alog10(x), x^(-0.65)/3.2

forprint, 1d9*agn_ssfr, redd, textout='SSFR_REDD.txt', /nocomm
stop

;y = comb_gam(x, -0.7, -3.2, -1.0, 2.0)
;oplot, alog10(x), 200*y, color=fsc_color('Red')
;forprint, x, y

;Calculate L_AGN
l_agn = redd*mbh*1.26e31*1d7
conv = (0.9/(0.1*9d16))*3600.*24.*365.25/!msun
mdot = l_agn*conv

sfr_agn = agn_ssfr*m

plot, l_agn, sfr_agn, /psym, /xlog, /ylog, yra=[1,1d4]

ll = 42+0.3*findgen(20)
ul = ll+0.3
plotsym, 0, 2, /fill
FOR i=0, n_elements(ll)-1 DO BEGIN
   o = WHERE(alog10(l_agn) GE ll[i] AND $
             alog10(l_agn) LT ul[i], n_o)
   IF n_o GT 0 THEN $
      plots, mean(l_agn[o]), mean(sfr_agn[o]), $
             psym=8, color=fsc_color('Red'), noclip=0
ENDFOR

stop

plot, redd, 1d9*agn_ssfr, /psym, /xlog, /ylog, xst=9;, $
      ;xra=[7d-4, 5], $
      ;yra=[0.7, 500], /yst
axis, xaxis=1, /xlog, xra=1d7*(10d^!x.crange)/conv, /xst

ll = -4.+0.3*findgen(20)
ul = ll+0.3
plotsym, 0, 2, /fill
FOR i=0, n_elements(ll)-1 DO BEGIN
   o = WHERE(alog10(mdot) GE ll[i] AND $
             alog10(mdot) LT ul[i], n_o)
   IF n_o GT 0 THEN $
      plots, mean(mdot[o]), mean(sfr_agn[o]), $
             psym=8, color=fsc_color('Red'), noclip=0
ENDFOR

ll = -2.+0.5*findgen(15)
ul = ll+0.5
FOR i=0, n_elements(ll)-1 DO BEGIN
   o = WHERE(alog10(sfr_agn) GE ll[i] AND $
             alog10(sfr_agn) LT ul[i], n_o)
   print, 10d^ll[i], 10d^ul[i], n_o
   IF n_o GT 0 THEN $
      plots, mean(mdot[o]), mean(sfr_agn[o]), $
             psym=8, color=fsc_color('Blue'), noclip=0
ENDFOR

oplot, 10d^!x.crange, (10d^!x.crange)*500, color=fsc_color('Red')


stop

window, 1
o = WHERE(sfr_agn GE 100)
plothist, alog10(redd), bin=0.01, /ylog
x = 10d^(-6.+0.01*findgen(700))
oplot, alog10(x), 0.02*x^(-0.65)



END
