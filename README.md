# ELS Pricing (Monte Carlo/Finite Difference Method)
 
![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)
![Jupyter](https://img.shields.io/badge/Jupyter-Notebook-orange?logo=jupyter)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Status](https://img.shields.io/badge/Status-IN--PROGRESS-yellow?style=flat-square&logo=github)

This repository contains Jupyter Notebooks for pricing **Equity-Linked Securities (ELS)** using **Monte Carlo simulation** and **Finite Difference Methods (FDM)**. The project models ELS pricing for varying numbers of underlying assets (1, 2, and 3) and is currently under development.

---

## Project Overview

The repository includes six Jupyter Notebooks, each focusing on a specific combination of asset count and pricing method, as outlined below:

|                | 1 Asset                            | 2 Assets                           | 3 Assets                           |
|----------------|------------------------------------|------------------------------------|------------------------------------|
| **Monte Carlo** | [A](./A_1Asset_MC.ipynb) | [B](./B_2Assets_MC.ipynb) | [C](./C_3Assets_MC.ipynb) |
| **Finite Difference** | [D](./D_1Asset_FDM.ipynb) | [E](./E_2Assets_FDM.ipynb) | [F](./F_3Assets_FDM.ipynb) |

---

## Notebooks

### Monte Carlo Simulation
1. **[A_1Asset_MC.ipynb](./A_1Asset_MC.ipynb)**  
   - Simulates stock price paths for a single asset using Geometric Brownian Motion (GBM).  
   - Includes multiple path simulations, final price distribution, and ELS pricing with early redemption and knock-in features.  
   - Visualizes asset paths and fair value evolution.

2. **[B_2Assets_MC.ipynb](./B_2Assets_MC.ipynb)**  
   - Extends Monte Carlo simulation to two correlated assets.  
   - Models ELS pricing with knock-in barriers and early redemption logic.  
   - Visualizes correlated price paths and fair value.

3. **[C_3Assets_MC.ipynb](./C_3Assets_MC.ipynb)**  
   - Simulates three correlated assets for ELS pricing.  
   - Implements complex knock-in and redemption structures.  
   - Includes visualizations for multi-asset paths and distributions.

### Finite Difference Method
4. **[D_1Asset_FDM.ipynb](./D_1Asset_FDM.ipynb)**  
   - Solves the Black-Scholes PDE for a single asset using FDM on a non-uniform grid.  
   - Computes European option prices and ELS fair values.  
   - Visualizes grid solutions and price evolution.

5. **[E_2Assets_FDM.ipynb](./E_2Assets_FDM.ipynb)**  
   - Applies FDM to price ELS with two correlated assets.  
   - Handles multi-dimensional PDEs for accurate pricing.  
   - Includes visualizations of price surfaces.

6. **[F_3Assets_FDM.ipynb](./F_3Assets_FDM.ipynb)**  
   - Extends FDM to three assets, solving high-dimensional PDEs.  
   - Computes ELS fair values with early redemption and knock-in features.  
   - Visualizes multi-dimensional results.

---

## Features

- **Monte Carlo Simulations**: Models asset price paths with GBM, supporting early redemption and knock-in structures for 1–3 assets.
- **Finite Difference Method**: Solves Black-Scholes PDEs on non-uniform grids for precise ELS pricing.
- **Real Market Data**: Integrates historical data (e.g., KOSPI 2018–2021) for realistic fair value calculations (in progress).
- **Visualizations**: Includes plots for price paths, final price distributions, and fair value evolution.
- **Modular Design**: Notebooks are structured for easy extension to additional assets or methods.

---

## Requirements

To run the notebooks, install the required Python packages:

```bash
pip install numpy pandas matplotlib scipy
```

---

## Status

This project is **in progress**. Current efforts focus on:
- Completing Monte Carlo simulations for three assets (`C_3Assets_MC.ipynb`).
- Refining FDM implementations for multi-asset cases (`E_2Assets_FDM.ipynb`, `F_3Assets_FDM.ipynb`).
- Integrating real market data across all notebooks.
- Enhancing visualizations and documentation.

Contributions and feedback are welcome!

---

## License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.

---
