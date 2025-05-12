# Vehicle Side-Slip Angle Estimation using a Discrete-Time Interval Estimator

This repository contains a MATLAB implementation of a **Discrete-Time Interval Estimator** for estimating the **side-slip angle** and **yaw rate** of a ground vehicle. The algorithm is based on the work published in the 2022 IFAC World Congress by **Hichem Bessafa**.

> 📖 Reference:  
> H. Bessafa, “Discrete-Time Interval Estimator for Vehicle Side-Slip Angle Estimation,” *IFAC-PapersOnLine*, 2022.  
> [DOI: 10.1016/j.ifacol.2022.11.294](https://doi.org/10.1016/j.ifacol.2022.11.294)

---

## 📌 Features

- Implements a **robust interval observer** for uncertain vehicle dynamics
- Handles **parametric uncertainty** using a **polytopic system model**
- Uses **YALMIP** and **SDPT3** for solving LMIs
- Provides visual comparison of real and estimated vehicle states

---

## 🛠️ Requirements

- MATLAB R2019b or later
- [YALMIP Toolbox](https://yalmip.github.io/download/)
- [SDPT3 Solver](https://github.com/sqlp/sdpt3)

---

## ⚙️ Installation & Setup

1. **Clone the repository**:

```bash
git clone https://github.com/yourusername/vehicle-slip-angle-estimation.git
cd vehicle-slip-angle-estimation
