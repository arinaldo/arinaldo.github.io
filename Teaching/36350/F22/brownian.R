
    # Discrete covariance function k(s,t) = s*(1-t) for s < t
    n = 1000; t = 1:n/n
    Sig = t %o% (1-t)
    Sig = pmin(Sig, t(Sig))

    # Compute symmmetric square root (via eigendecomposition)
    eig = eigen(Sig)
    Sig.half = eig$vec %*% diag(sqrt(eig$val)) %*% t(eig$vec)

    # Generate discrete Brownian motion and plot it
    B = Sig.half %*% rnorm(n)
    plot(t, B, type="l")








