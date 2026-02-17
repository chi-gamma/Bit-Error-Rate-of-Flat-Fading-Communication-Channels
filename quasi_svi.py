import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize


def svi_raw(parameters, k):
    """ SVI raw parametrization for total implied variance
    k: log forward moneyness for which total variance is evaluated
    
    :parameter a:       vertical translation of the smile
    :parameter b:       tightens the smile
    :parameter m:       translates the smile to the right
    :parameter rho:     counterclockwise rotation of the smile
    :parameter sigma:   controls the ATM curvature of the smile
    """
    a, b, m, rho, sigma = parameters
    return a + b * ( rho * (k - m) + np.sqrt((k - m)**2 + sigma**2) )


def solve_reduced_ls(X, y):
    """
    Solve least squares safely for possibly reduced matrix.
    """
    return np.linalg.lstsq(X, y, rcond=None)[0]


def svi_active_set(y_vec, v, tenor):
    """
    Active-set solver for:
        v ≈ a + d*y + c*sqrt(1+y^2)

    Returns:
        dict with a, d, c
    """

    g = np.sqrt(1.0 + y_vec**2)

    # Design matrix columns: [a, d, c]
    X_full = np.column_stack([
        np.ones_like(y_vec),
        y_vec,
        g
    ])

    a_max = np.max(v)
    c_max = 4.0 / tenor

    # Initial unconstrained solution
    beta = solve_reduced_ls(X_full, v)
    a, d, c = beta

    # Active set bookkeeping
    active = {
        "a": None,  # None means free, otherwise fixed value
        "d": None,
        "c": None
    }

    max_iter = 10
    for _ in range(max_iter):

        violated = False

        # ---- Check a ----
        if active["a"] is None:
            if a < 0.0:
                active["a"] = 0.0
                violated = True
            elif a > a_max:
                active["a"] = a_max
                violated = True

        # ---- Check c ----
        if active["c"] is None:
            if c < 0.0:
                active["c"] = 0.0
                violated = True
            elif c > c_max:
                active["c"] = c_max
                violated = True

        # ---- Check d (coupled constraint) ----
        # Only meaningful once c is known (fixed or current value)
        current_c = active["c"] if active["c"] is not None else c
        d_bound = min(current_c, c_max - current_c)

        if active["d"] is None:
            if d > d_bound:
                active["d"] = d_bound
                violated = True
            elif d < -d_bound:
                active["d"] = -d_bound
                violated = True

        # If no violations → done
        if not violated:
            break

        # Build reduced system
        free_vars = []
        fixed_vals = {}

        for i, name in enumerate(["a", "d", "c"]):
            if active[name] is None:
                free_vars.append(i)
            else:
                fixed_vals[i] = active[name]

        # Construct reduced regression
        if len(free_vars) > 0:
            X_free = X_full[:, free_vars]
            y_adj = v.copy()

            # subtract fixed contributions
            for idx, val in fixed_vals.items():
                y_adj -= X_full[:, idx] * val

            beta_free = solve_reduced_ls(X_free, y_adj)

            # Reconstruct full beta
            beta_new = np.zeros(3)

            j = 0
            for i in range(3):
                if i in free_vars:
                    beta_new[i] = beta_free[j]
                    j += 1
                else:
                    beta_new[i] = fixed_vals[i]

            a, d, c = beta_new
        else:
            # all variables fixed
            a = active["a"]
            d = active["d"]
            c = active["c"]
            break

    return a, d, c


def quasi_explicit_svi_fit(k_vec, vol_vec, tenor, cal_k_grid=None, prev_slice_params=None, init_guess=None):

    # -------------------------------
    # 0. Sanity checks & preprocessing
    # -------------------------------
    k_vec = np.asarray(k_vec)
    vol_vec = np.asarray(vol_vec)
    
    if prev_slice_params is not None:
        cal_k_grid = np.asarray(cal_k_grid) if cal_k_grid is not None else np.linspace(-1.5, 1.5, 50)
    
    # Convert vols -> total implied variance
    w_vec = tenor * vol_vec**2
    
    def sub_params_to_full(sub_params):
        m, sigma = sub_params
        y_vec = (k_vec - m) / sigma
        a_hat, d_hat, c_hat = svi_active_set(y_vec, w_vec, tenor)
        rho = d_hat / c_hat
        b = c_hat / sigma
        return (a_hat, b, m, rho, sigma)
    
    
    def obj_func(sub_params, cal_k_grid, prev_slice_params):
        m, sigma = sub_params
    
        if sigma <= 0:
            return 1e10
        
        y_vec = (k_vec - m) / sigma
        a_hat, d_hat, c_hat = svi_active_set(y_vec, w_vec, tenor)
        w_model = a_hat + d_hat*y_vec + c_hat*np.sqrt(1+y_vec**2)
        sse = sum((w_model - w_vec)**2)
        
        if prev_slice_params is not None:
            y_cal = (cal_k_grid - m) / sigma
            w_curr = a_hat + d_hat*y_cal + c_hat*np.sqrt(1 + y_cal**2)
            w_prev = svi_raw(prev_slice_params, cal_k_grid)
            crossedness = max(0.0, np.max(w_prev - w_curr))
            sse += 10e10 * crossedness**2
            
        return sse
            
    def calibrate(init_guess, cal_k_grid, prev_slice_params):  
        """ 
        finds optimal raw SVI parameters that fit martet prices
        
        """
        
        bnds = [(2 * min(k_vec), 2 * max(k_vec)), (1e-5, 1)]
        if init_guess is None:
            # sample initial guess
            N = 20
            sse = float('inf')
            sub_params = None
            for _ in range(N):
                x0 = [np.random.uniform(bnd[0], bnd[1]) for bnd in bnds]
                res = minimize(obj_func, x0, method='Nelder-Mead')
                if res.fun < sse:
                    sse = res.fun
                    sub_params = res.x
        else:
            res = minimize(obj_func, init_guess, method='Nelder-Mead')
            sub_params = res.x
            
        full_params  = sub_params_to_full(sub_params)
            
        return full_params

    
    svi_params = calibrate(init_guess, cal_k_grid, prev_slice_params)
    fitted_w = svi_raw(svi_params, k_vec)
    fitted_vols = np.sqrt(fitted_w / tenor)  

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

svi_params, fitted_vol_vec = quasi_explicit_svi_fit(k_vec, vol_vec, tenor)

k_other = np.linspace(-1.2, 1.0, 100)
other = np.sqrt(svi_raw(svi_params, k_other) / tenor)


plt.plot(k_other, other, '*', label='full')
plt.plot(k_vec, vol_vec, 'o', label='input vols')
plt.plot(k_vec, fitted_vol_vec, '-', label='fitted vols')
plt.xlabel('log-forward moneyness')
plt.ylabel('vol')
plt.title('SVI Quasi Calibration')
plt.legend()
plt.show()  
        

    

