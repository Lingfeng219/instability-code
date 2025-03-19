# Baroclinic Eigen Mode Solver

This MATLAB function solves the **baroclinic linear normal mode instability problem** using the **matrix method (Smith et al., 2007)**. It includes submesoscale eddy mixing effects.

## 📌 Usage
```matlab
[B, A, om, phi] = baroclinic_eigen_mode(u, v, rho, k, l, lat, z, alphax, alphay, rho0, visc, qx, qy);
