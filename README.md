# Pantheon-Cosmology-MCMC
MCMC parameter estimation for cosmological models using the Pantheon Supernovae dataset.
# Cosmological Parameter Estimation using Pantheon Dataset

## Overview
This repository contains a Python-based cosmology project aimed at estimating the Hubble constant ($H_0$) and the matter density parameter ($\Omega_m$). The estimation is performed by fitting theoretical models to the **Pantheon Type Ia Supernovae dataset**.

## Methodology
The core of this project relies on **Markov Chain Monte Carlo (MCMC)** simulations to explore the parameter space and find the best-fit values. 
- **Dataset:** Pantheon Sample (Apparent magnitude vs. Redshift).
- **Libraries used:** `emcee` for the MCMC sampling, and `corner` for visualizing the posterior distributions.
- **Physical Model:** Flat $\Lambda$CDM cosmology.

## Project Structure
- `mcmc_64walkers.py`: The main Python script containing the standardized code (with type hinting and docstrings) to run the MCMC simulation using 64 walkers.
- `pantheon.txt`: The dataset containing redshift and distance modulus data.
- `hubble_diagram.png`: The plot showing the best-fit theoretical distance modulus curve against the Pantheon observational data.
- `corner_plot_32w.png` / `corner_plot_64w.png`: Corner plots visualizing the marginalized constraints on $H_0$ and $\Omega_m$ using 32 and 64 walkers respectively.

## Requirements
To run this code, you need Python 3.x and the following libraries:
```bash
pip install numpy matplotlib scipy emcee corner

## How to Run
Simply execute the main python script in your terminal:

bash
python mcmc_64walkers.py

## Results
The MCMC analysis successfully converges to estimate the cosmological parameters. The comparison between 32 and 64 walkers demonstrates the stability and accuracy of the sampling process. (Check the attached `.png` files for visual results).

