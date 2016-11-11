FUNCTION rr_ms, n, mu, sigma, seed=seed

log_r_ms = mu+sigma*randomn(seed, n)

RETURN, 10d^log_r_ms

END
