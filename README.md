# Vehicle Side-Slip Angle Estimation using a Discrete-Time Interval Estimator

This repository provides a MATLAB implementation of a **Discrete-Time Interval Estimator** for robust estimation of a vehicle's **side-slip angle** and **yaw rate**, accounting for model uncertainties. The method is based on the work presented at the IFAC World Congress 2022.

> 📖 Reference:  
> H. Bessafa, “Discrete-Time Interval Estimator for Vehicle Side-Slip Angle Estimation,” *IFAC-PapersOnLine*, 2022.  
> [DOI: 10.1016/j.ifacol.2022.11.294](https://doi.org/10.1016/j.ifacol.2022.11.294)

---

## 📌 Features

- Implements a **robust interval observer** for uncertain vehicle dynamics.
- Includes two estimator configurations:
  - **Algorithm 1 (Algo1)**: Uses LMI optimization via **YALMIP** and **SDPT3**.
  - **Algorithm 2 (Algo2)**: Uses **pole placement** design (no solver needed).
- Handles **polytopic system representation** of parametric uncertainty.
- Includes 2D state estimation plots and bounded error visualization.

---

## 🛠️ Requirements

- MATLAB R2019b or later
- [YALMIP Toolbox](https://yalmip.github.io/download/)
- [SDPT3 Solver](https://github.com/sqlp/sdpt3) (only required for Algo1)

---

## ⚙️ Installation & Setup

1. **Clone the repository**:

```bash
git clone https://github.com/yourusername/vehicle-slip-angle-estimation.git
cd vehicle-slip-angle-estimation
