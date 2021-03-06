FUNCTION randomf, seed, dims, x, y, double=double

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Returns a sample of random numbers the frequency of;;
;;which obeys the function described by x and y.     ;;
;;                                                   ;;
;;This uses the interpol function to interpolate the ;;
;;function given, so the given function should be    ;;
;;smoothly varying                                   ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ran = randomu(seed, dims, double=double)

;Cumulative must always be less than 1:
cum = total(y, /cum)/max(total(y, /cum))

RETURN, interpol(x, cum, ran)

END
