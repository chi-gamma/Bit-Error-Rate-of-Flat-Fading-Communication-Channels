import numpy as np
import matplotlib.pyplot as plt



def conic_svi_fit(k_vec, vol_vec, tenor, weights=None):
    """
    Fits SVI Raw parameters to the volatility smile for a given tenor using 
    the direct least-squares conic method (hyperbolic branch).

    Parameters
    ----------
    k_vec : np.ndarray
        Log-forward moneyness vector (log(strike / forward)).
    vol_vec : np.ndarray
        Black-Scholes implied volatilities corresponding to k_vec.
    tenor : float
        Option time to maturity (in years).
    weights : np.ndarray, optional
        Optional vector of weights for each strike. If None, all weights = 1.

    Returns
    -------
    svi_params : tuple
        Fitted SVI raw parameters (a, b, m, rho, sigma).
    fitted_vols : np.ndarray
        Volatility smile reconstructed from the fitted SVI parameters.
    """

    # -------------------------------
    # 0. Sanity checks & preprocessing
    # -------------------------------
    k_vec = np.asarray(k_vec)
    vol_vec = np.asarray(vol_vec)
    n_strikes = len(k_vec)

    if weights is not None:
        weights = np.asarray(weights).ravel()
        if len(weights) != n_strikes:
            raise ValueError("weights must have the same length as k_vec / vol_vec")
    else:
        weights = np.ones(n_strikes)  # uniform weighting

    # Convert vols -> total implied variance
    w_vec = tenor * vol_vec**2

    # -------------------------------
    # 1. Build design matrix D
    # -------------------------------
    # D = [k^2, w^2, k*w, k, w, 1]
    D = np.column_stack((k_vec**2, w_vec**2, k_vec*w_vec, k_vec, w_vec, np.ones(n_strikes)))

    # -------------------------------
    # 2. Compute scatter matrix S efficiently
    # -------------------------------
    # S = D.T @ W @ D, but W is diagonal
    # Equivalent to multiplying each row of D by weight first
    S = (D.T * weights) @ D

    # Partition scatter matrix into constrained (c) and unconstrained (u) parts
    Scc = S[0:2, 0:2]    # 2x2, quadratic coefficients
    Suc = S[2:6, 0:2]    # 4x2, cross terms
    Suu = S[2:6, 2:6]    # 4x4, unconstrained

    # -------------------------------
    # 3. Solve for constrained coefficients zc
    # -------------------------------
    # Compute the "reduced" matrix for eigenproblem
    M = Scc - Suc.T @ np.linalg.inv(Suu) @ Suc
    M = 0.5 * (M + M.T)  # enforce symmetry for numerical stability

    # Choose the hyperbolic branch: z1 = -sqrt(M[1,1]/M[0,0])
    d = M[1, 1] / M[0, 0]
    if d <= 0:
        raise ValueError("Non-positive discriminant, check input data")
    z1 = -np.sqrt(d)
    zc = np.array([z1, 1.0])

    # -------------------------------
    # 4. Recover unconstrained coefficients zu
    # -------------------------------
    zu = -np.linalg.inv(Suu) @ Suc @ zc
    zs = np.concatenate([zc, zu])  # full conic coefficient vector

    # -------------------------------
    # 5. Helper: reconstruct total variance w(k) from conic coefficients
    # -------------------------------
    def z_to_w(zs, x):
        z1, z2, z3, z4, z5, z6 = zs
        term1 = z3 * x + z5
        term2 = z1 * x**2 + z4 * x + z6
        discriminant = term1**2 - 4 * z2 * term2
        if np.any(discriminant < 0):
            raise ValueError("Negative discriminant in z_to_w, check conic coefficients")
        w = (-term1 + np.sqrt(discriminant)) / (2 * z2)
        return w

    # -------------------------------
    # 6. Helper: convert conic zs -> SVI raw parameters
    # -------------------------------
    def z_to_svi(zs):
        z1, z2, z3, z4, z5, z6 = zs

        # b must be positive
        disc_b = 0.25 * z3**2 - z1
        if disc_b <= 0:
            raise ValueError("Computed b^2 is non-positive; check input or branch selection")
        b = np.sqrt(disc_b)

        rho = -z3 / (2 * b)
        m = (z4 + b * rho * z5) / (2 * b**2)
        a = b * rho * m - 0.5 * z5

        sigma_num = 0.25 * z5**2 - b**2 * m**2 - z6
        if sigma_num < 0:
            raise ValueError("Negative value under sqrt for sigma; check input or data")
        sigma = np.sqrt(sigma_num) / b

        return (a, b, m, rho, sigma)

    # -------------------------------
    # 7. Compute fitted vols and SVI parameters
    # -------------------------------
    fitted_w = z_to_w(zs, k_vec)
    fitted_vols = np.sqrt(fitted_w / tenor)
    svi_params = z_to_svi(zs)

    return svi_params, fitted_vols








spot = 6941.81
forward = 7257.13229
discount = 0.948667
tenor = 1.6



strikes = np.array(
    [2082.543, 2429.6335, 2776.724, 3123.8145, 3470.905, 3817.9955, 4165.086, 4512.1765,
      4859.267, 5206.3575, 5553.448, 5900.5385, 6247.629, 6594.7195, 6941.81, 7288.9005,
      7635.991, 7983.0815, 8330.172, 8677.2625, 9024.353, 9371.4435, 9718.534, 10065.6245, 10412.715]
    )

vol_vec = np.array(
        [0.465101, 0.433829, 0.406067, 0.380986, 0.358017, 0.336751, 0.316882, 0.298179, 0.28046,
        0.263585, 0.247445, 0.231959, 0.217073, 0.202771, 0.189078, 0.176086, 0.163985, 0.153094,
        0.143868, 0.136788, 0.132136, 0.129803, 0.129362, 0.130299, 0.13217]
        )


k_vec = np.log(strikes / forward)

svi_params, fitted_vol_vec = conic_svi_fit(k_vec, vol_vec, tenor, weights=None)


plt.plot(k_vec, vol_vec, 'o', label='input vols')
plt.plot(k_vec, fitted_vol_vec, '-', label='fitted vols')
plt.xlabel('log-forward moneyness')
plt.ylabel('vol')
plt.legend()
plt.show()