# ELS Monte Carlo Pricing Simulator

![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)
![Jupyter](https://img.shields.io/badge/Jupyter-Notebook-orange?logo=jupyter)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Status](https://img.shields.io/badge/Status-IN--PROGRESS-yellow?style=flat-square&logo=github)

This repository contains a Jupyter Notebook for pricing **Equity-Linked Securities (ELS)** using the **Monte Carlo simulation method**. The notebook (`ELS_MonteCarlo_1.ipynb`) walks through the entire process of modeling, simulating, and visualizing the fair value of an ELS product, both in theoretical settings and with real market data.

---

## Notebook Overview

### `ELS_MonteCarlo_1.ipynb`

| Section | Description |
|--------|-------------|
| **1. GBM Simulation** | Simulates a single stock price path using Geometric Brownian Motion (GBM). |
| **2. GBM Simulation (Multiple)** | Generates and visualizes multiple GBM paths to show stochastic variability. |
| **3. Final Price Distribution** | Simulates 10,000 paths and visualizes the distribution of final prices. |
| **4. Finite Difference Method** | Solves the Black-Scholes PDE using FDM on a non-uniform grid to compute European call option price. |
| **5. ELS Pricing via Monte Carlo** | Full ELS pricing model with early redemption and knock-in logic under simulated paths. |
| **6. ELS Fair Value Example** | Uses actual KOSPI data (2018–2021) to calculate daily ELS fair values. |

---

## Features

- **Monte Carlo Simulation** of asset price paths with early redemption and knock-in structure.
- **Black-Scholes PDE Solution** using finite difference method on a non-uniform grid.
- **Real Market Data Integration** for historical KOSPI index.
- **Daily Fair Value Tracking** of ELS across a 3-year time span.
- Clear **visualizations** for asset paths, price distributions, and ELS fair value evolution.

---

## Requirements

Install the required Python packages before running the notebook:

```bash
pip install numpy pandas matplotlib
```

IN PROGRESS
