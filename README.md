# Simulation and Modelling of Dynamic Systems

This repository contains a series of projects developed for the *Simulation and Modelling of Dynamic Systems* course, focusing on system identification, real-time parameter estimation, nonlinear modeling, and control.

## 🔧 Technologies
- MATLAB
- ODE Solvers
- System Identification
- Control Theory
- Numerical Methods

---

## 📁 Project Overview

The project is divided into three main parts, covering both linear and nonlinear dynamic systems.

---

### ⚙️ Part 1 – Parameter Estimation (Least Squares)
- Modeled a linearized pendulum system in state-space form.
- Simulated system dynamics using MATLAB ODE solvers.
- Applied Least Squares estimation to identify unknown parameters (mass, damping, length).
- Evaluated the impact of:
  - Measurement noise
  - Sampling period
  - Input signal characteristics
- Compared true and estimated system responses. :contentReference[oaicite:0]{index=0}

---

### 📡 Part 2 – Real-Time Parameter Estimation
- Implemented real-time adaptive estimation methods for a dynamic system:
  - Gradient-based estimation
  - Lyapunov-based adaptive control
- Designed parallel and composite (mixed) estimator structures.
- Analyzed:
  - Convergence behavior
  - Stability via Lyapunov functions
  - Sensitivity to noise and disturbances
- Demonstrated that composite estimators achieve higher accuracy, while parallel structures are more robust to noise. :contentReference[oaicite:1]{index=1}

---

### 🔄 Part 3 – Nonlinear System Modeling & Control
- Modeled unknown nonlinear dynamic systems using parametric representations.
- Performed model structure selection using basis functions.
- Applied cross-validation for model evaluation and selection.
- Designed a nonlinear control law for trajectory tracking.
- Ensured system stability using Lyapunov-based analysis.
- Evaluated system performance under external disturbances and uncertainties. :contentReference[oaicite:2]{index=2}

---

## 🚀 Getting Started
1. Open the MATLAB files (`.m`) provided in the repository.
2. Run simulations for each part separately.
3. Visualize system states, parameter estimates, and error dynamics.

---

## 📊 Key Concepts Covered
- Least Squares Estimation
- Adaptive Control & Real-Time Estimation
- Lyapunov Stability Analysis
- Nonlinear System Modeling
- Model Selection & Validation
- Robustness to Noise and Disturbances

---

## 📌 Highlights
- Developed real-time estimators with provable stability.
- Achieved accurate parameter convergence under different inputs.
- Analyzed trade-offs between robustness and accuracy.
- Designed controllers for nonlinear trajectory tracking.

---

## 👤 Author
Christos Alexopoulos
