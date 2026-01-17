# Antennas and Propagation Projects: Design, Simulation & Optimization

[![MATLAB](https://img.shields.io/badge/MATLAB-Simulation-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![NEC2](https://img.shields.io/badge/4NEC2-Antenna_Modeling-green.svg)](http://www.qsl.net/4nec2/)
[![University](https://img.shields.io/badge/AUTH-ECE-red.svg)](https://www.ece.auth.gr/en/)

This repository contains a comprehensive suite of antenna design and analysis projects developed for the **Antennas and Propagation** course at the **Aristotle University of Thessaloniki (AUTH)**. The project covers theoretical modeling, MATLAB numerical analysis, and professional antenna simulation using the NEC (Numerical Electromagnetics Code) engine.

---

## 🛰️ Project Overview

The project is divided into two main computational phases:

### Phase A: MATLAB Numerical Analysis & Array Optimization
Focused on the mathematical modeling of antenna arrays and the optimization of radiation characteristics.
* **2D Dipole Arrays:** Analysis of a $24 \times 12$ array of $\lambda/2$ dipoles. Calculation of radiation patterns and beam steering capabilities.
* **Parasitic Elements & Mutual Impedance:** Modeling a 3-element array with parasitic dipoles, calculating mutual impedance and complex input impedance.
* **Optimization:** Using **Genetic Algorithms (GA)** to minimize Side Lobe Levels (SLL) and optimizing reflection coefficients $|\Gamma|$ relative to element spacing and reflector distance.
* **Infinite Reflectors:** Application of the Method of Images to simulate the effect of an infinite reflector on antenna performance.

### Phase B: Computational Analysis with 4NEC2
Focused on the simulation of complex wire antennas using the Method of Moments (MoM).
* **Discone Antennas:** Analysis of broadband characteristics, impedance matching ($Z_{in}$), and radiation patterns across a wide frequency spectrum ($f_0$ to $4f_0$).
* **Folded Dipoles:** Study of the effect of conductor spacing on resonance and input impedance.
* **Traveling Wave Antennas:** Design of a horizontal wire antenna, evaluating the impact of termination resistance ($R_L$) and ground conditions (Perfect vs. Real Ground) on gain and front-to-back ratio.

---

## 🛠️ Tools & Technologies
* **MATLAB:** Custom scripts for array factor calculations, impedance matrices, and Genetic Algorithm optimization.
* **4NEC2:** Professional tool for MoM-based simulation and 3D radiation pattern visualization.
* **LaTeX:** Technical reporting and mathematical documentation.

---

## 📊 Key Findings
* Identified optimal matching regions for parasitic arrays ($|\Gamma| < 0.3$) for specific $d/\lambda$ and $h/\lambda$ ratios.
* Demonstrated how termination resistance $R_L$ eliminates standing waves in traveling wave antennas, achieving broadband performance.
* Optimized side lobe reduction through non-uniform element spacing and excitation.

---

## 📂 Repository Structure
```text
├── Phase_A_MATLAB/
│   ├── radiation_patterns.m    # 2D/3D radiation plots
│   ├── genetic_optimizer.m     # SLL minimization script
│   └── mutual_impedance.m      # Complex impedance calculations
├── Phase_B_NEC/
│   ├── discone_antenna.nec     # NEC input file for discone
│   ├── folded_dipole.nec       # NEC input file for folded dipole
│   └── traveling_wave.nec      # NEC input file for long-wire antenna
└── Reports/
    ├── Part_A_Analysis.pdf     # Detailed theoretical and MATLAB report
    └── Part_B_Simulation.pdf   # NEC results and performance evaluation
