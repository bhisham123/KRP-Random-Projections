# KRP-Random-Projection

This repository contains MATLAB codes to reproduce the results presented in the paper: 

**"Improved Analysis of Khatri-Rao Random Projections and Applications"**  by *Arvind K. Saibaba, Bhisham Dev Verma, and Grey Ballard.* (Submitted, 2025)

## 🔍 Folder Structure & Usage

- **`test_SysID.m`**  
  Main script to reproduce the experiments in **Section 4** (Compressing Block-structured matrices).  
  Depends on code and data in the **`sysID/`** folder.

- **`test_cauchy.m`**  
  Main script to reproduce the experiments in **Section 5** (Tucker compression with Khatri-Rao Products).  
  Depends on code in the **`tucker/`** folder.

- **`test_sensor_placement.m`**  
  Main script to reproduce the experiments in **Section 6** (Sensor Placement).  
  Depends on code in the **`sensor_placement/`** folder.
  (To run the code first download the simulation data from the following url: https://cgl.ethz.ch/Downloads/Data/ScientificData/tangaroa3d_nc.zip)

- **`Experimental_results/`**  
  Contains data and results corresponding to the figures in the paper.

---

## 📦 Dependencies

This code requires the following MATLAB toolboxes:

- **Tensor Toolbox** – Already included in the repository (tensor_toolbox).
- **Control System Toolbox** – Required for system identification experiments.
