# Dissertation Project: Radial Acceleration Relation (RAR) in Dwarf Galaxies

## Overview

This repository contains the code I wrote for my Final Year Project (dissertation) as part of my Physics BSc degree at the University of Surrey, supervised by Justin Read. The project focuses on investigating the Radial Acceleration Relation (RAR) within dwarf galaxies using data from the EDGE simulations.

## Data

The EDGE simulation data used within my project is not included in this repo as it is not mine to distribute. Information about the EDGE simulations can be found at [EDGE Simulation Website](https://edge-simulation.github.io/). Specifically, the 5 galaxy simulations used throughout my research are taken from the 2021 Matthew Orkney et al paper ([DOI: 10.1093/mnras/stab1066](https://doi.org/10.1093/mnras/stab1066)).

## Files

- `EDGEFYPMainCode.py`: Main script for compiling, storing, and visualizing data. Some calculations are completed here.
- `EDGELargeFileReducer.py`: Script for reducing large raw data files from EDGE simulations.
- `EDGEGCalcFINAL.py`: Script for completing gravitational acceleration calculations using refined data files.

## Methodology and Packages Used

- Functions were written for reading in raw/processed data and major calculations.
- **Packages Used:**
  - NumPy: For efficient array operations.
  - Matplotlib: For data visualization.

## Results

The results of my research show that dark matter is favored over MOND in explaining the dynamics seen within small galaxies. The latest observational data (Read et al, 2019) aligns much closer to the results from the EDGE simulations than the analytical MOND RAR.

### Baryonic vs Observed Gravitational Acceleration
![Baryonic vs Observed Gravitational Acceleration](https://github.com/tgholmes/FinalYearProject-RAR-EDGE/assets/148396727/8f35614c-41a2-496c-82e9-23a0d0c8c9ab)

### Residuals from MOND
![Residuals from MOND](https://github.com/tgholmes/FinalYearProject-RAR-EDGE/assets/148396727/1fce5559-d959-4f4e-94b9-c33e423a0635)

## Acknowledgments

I would like to thank my supervisor, Justin Read, for his continued support and guidance throughout this project. Additionally, I acknowledge the creators of the EDGE simulations and the authors of the referenced papers for making their data available.

## Contact

If you would like to contact me further about this project and the technical details/context, I would be more than happy to chat. [LinkedIn](www.linkedin.com/in/tom-holmes-278826214)
