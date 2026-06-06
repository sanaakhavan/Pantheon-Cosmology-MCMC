import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import emcee
import corner
from scipy.integrate import quad
from typing import Union, List, Tuple

# ==========================================
# GLOBAL VARIABLES (Configuration)
# ==========================================
PANTHEON_URL = "https://raw.githubusercontent.com/dscolnic/Pantheon/master/lcparam_full_long.txt"
M_FIXED = -19.35        # Absolute magnitude
C_LIGHT = 299792.458    # Speed of light in km/s

# MCMC Settings
N_WALKERS = 64
N_STEPS = 4000
DISCARD_STEPS = 1000
THIN_STEPS = 15
INITIAL_GUESS = [70.0, 0.3]  # [H0, Omega_m]

# ==========================================
# PHYSICS & MATH FUNCTIONS
# ==========================================


def luminosity_distance(z: Union[float, np.ndarray], H0: float, Om: float) -> Union[float, np.ndarray]:
    """
    Calculate the theoretical luminosity distance d_L for a Flat Lambda-CDM model.

    Args:
        z (float or np.ndarray): Redshift(s).
        H0 (float): Hubble constant in km/s/Mpc.
        Om (float): Matter density parameter (Omega_m).

    Returns:
        float or np.ndarray: Luminosity distance in Mpc.
    """
    OL = 1.0 - Om

    def integrand(x):
        return 1.0 / np.sqrt(Om * (1 + x)**3 + OL)

    if np.isscalar(z):
        integ, _ = quad(integrand, 0, z)
        return (C_LIGHT * (1 + z) / H0) * integ
    else:
        res = []
        for zi in z:
            i, _ = quad(integrand, 0, zi)
            res.append((C_LIGHT * (1 + zi) / H0) * i)
        return np.array(res)


def get_distance_modulus_model(params: List[float], z_array: np.ndarray) -> np.ndarray:
    """
    Calculate the theoretical distance modulus mu for given parameters.

    Args:
        params (list): Cosmological parameters [H0, Omega_m].
        z_array (np.ndarray): Array of redshifts.

    Returns:
        np.ndarray: Theoretical distance modulus array.
    """
    H0, Om = params
    dL = luminosity_distance(z_array, H0, Om)

    mu_model = np.zeros_like(dL)
    mask = dL > 0
    mu_model[mask] = 5 * np.log10(dL[mask]) + 25
    mu_model[~mask] = np.nan
    return mu_model

# ==========================================
# MCMC FUNCTIONS
# ==========================================


def log_prior(params: List[float]) -> float:
    """
    Define the prior probability for the parameters.
    H0 must be between 50 and 100, Omega_m between 0 and 1.
    """
    H0, Om = params
    if 50.0 < H0 < 100.0 and 0.0 < Om < 1.0:
        return 0.0
    return -np.inf


def log_likelihood(params: List[float], z: np.ndarray, mu_data: np.ndarray, err: np.ndarray) -> float:
    """
    Calculate the log-likelihood comparing the theoretical model to the data.
    """
    mu_model = get_distance_modulus_model(params, z)
    if not np.all(np.isfinite(mu_model)):
        return -np.inf
    sigma2 = err**2
    return -0.5 * np.sum((mu_data - mu_model)**2 / sigma2)


def log_probability(params: List[float], z: np.ndarray, mu_data: np.ndarray, err: np.ndarray) -> float:
    """
    Calculate the total log-probability (prior + likelihood).
    """
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, z, mu_data, err)

# ==========================================
# MAIN EXECUTION BLOCK
# ==========================================


def main():
    # 1. Load Data
    print(">>> Loading Pantheon Data...")
    try:
        df = pd.read_csv(PANTHEON_URL, sep=r'\s+')
        z_obs = df['zcmb'].values
        mb_obs = df['mb'].values
        dmb_obs = df['dmb'].values

        # Calculate observational distance modulus
        mu_obs = mb_obs - M_FIXED
        print(f"OK: Loaded {len(z_obs)} Supernovae.")
    except Exception as e:
        print("Error loading data:", e)
        return

    # 2. Setup MCMC
    ndim = len(INITIAL_GUESS)
    print(f">>> Running MCMC with {N_WALKERS} walkers for {N_STEPS} steps...")
    print("    (This might take 5-10 minutes due to heavy calculations, please wait...)")

    # Initial tight cluster of walkers around the initial guess
    pos = INITIAL_GUESS + 1e-3 * np.random.randn(N_WALKERS, ndim)

    sampler = emcee.EnsembleSampler(
        N_WALKERS, ndim, log_probability, args=(z_obs, mu_obs, dmb_obs)
    )

    # Run MCMC
    sampler.run_mcmc(pos, N_STEPS, progress=True)

    # 3. Process Results
    flat_samples = sampler.get_chain(
        discard=DISCARD_STEPS, thin=THIN_STEPS, flat=True)
    best_params = np.median(flat_samples, axis=0)

    print("\nRESULTS (High Precision):")
    labels = ["H0", "Omega_m"]
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(f"{labels[i]}: {mcmc[1]:.3f} (+{q[1]:.3f} / -{q[0]:.3f})")

    # 4. Plots
    # 4.1 Corner Plot
    fig = corner.corner(
        flat_samples, labels=labels, show_titles=True,
        color="darkblue",
        quantiles=[0.16, 0.5, 0.84],
        title_fmt='.3f',
        plot_density=True,
        plot_contours=True,
        levels=(1 - np.exp(-0.5), 1 - np.exp(-2.0)),
        smooth=1.0
    )
    fig.suptitle(
        f"Pantheon Constraints (Fixed M={M_FIXED})\nHigh Resolution Run", fontsize=14)
    plt.show()

    # 4.2 Hubble Diagram
    plt.figure(figsize=(10, 6))
    plt.errorbar(z_obs, mu_obs, yerr=dmb_obs, fmt=".k",
                 alpha=0.05, label="Pantheon Data ($m_B - M$)")

    z_fit = np.linspace(0.01, 2.3, 200)
    mu_fit = get_distance_modulus_model(best_params, z_fit)

    plt.plot(z_fit, mu_fit, "r-", lw=2.5,
             label=f"Best Fit ($\Omega_m$={best_params[1]:.3f})")

    plt.xlabel("Redshift (z)")
    plt.ylabel("Distance Modulus ($\mu$)")
    plt.title(
        f"Hubble Diagram\n$H_0$={best_params[0]:.2f}, $\Omega_m$={best_params[1]:.3f}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()


if __name__ == "__main__":
    main()
