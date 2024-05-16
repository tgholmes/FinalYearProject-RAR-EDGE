# Dissertation Project: Radial Acceleration Relation (RAR) in Dwarf Galaxies

## Overview

This repository contains the code I wrote for my Final Year Project (dissertation) as part of my Physics BSc degree at the University of Surrey, supervised by Justin Read. The project focuses on investigating the Radial Acceleration Relation (RAR) within dwarf galaxies using data from the EDGE simulations. The abstract from my dissertation report:

The radial acceleration relation (RAR) describes the relation between the baryonic and observed gravitational accelerations, 𝑔!"# and 𝑔$!% respectively, seen throughout a galaxy, at different radii. It has been proposed as a law that is followed by all galaxies and has been claimed to be a consequence of Milgrom’s Modified Newtonian Dynamics (MOND), an alternate theory removing the necessity for dark matter to explain galactic dynamics, such as rotation curves. In this project, I use the high-resolution Engineering-Dwarfs-at-Galaxy-Formations- Edge (EDGE) simulations of the smallest galaxies to predict the RAR on acceleration scales an order of magnitude lower than has been probed to date. Baryonic and observed gravitational accelerations were calculated for 5 simulations, constructed upon the cosmological parameters from the Planck Collaboration (2013), rooted in the Λ-Cold Dark Matter (𝛬CDM) model, and stellar masses between 105 𝑀⨀< 𝑀∗ < 107 𝑀⨀. I find that all 5 simulations do not lie along the RAR model suggested from MOND. Instead, they exhibit higher observed accelerations than predicted by the RAR at the same baryonic counterpart. ΛCDM can account for these higher observed values, whereas the MONDian interpretation cannot explain these different values that do not fall on the RAR. Comparing these findings with data for nearby dwarf spheroidal galaxies, I show that they are more well described by the ΛCDM EDGE simulations than the MOND regime.

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
