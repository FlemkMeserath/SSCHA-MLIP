# SSCHA-MLIP Library

This repository provides a library for the use of the **Stochastic Self-Consistent Harmonic Approximation (SSCHA)** method coupled with the active learning procedure used to train **Moment Tensor Potentials (MTPs)** via the **MLIP version 2** code. This method reduces the computational load of performing SSCHA calculations. Notably, we have extended support for self-consistent density functional theory (DFT) calculations performed using both the **Quantum ESPRESSO (QE)** package and the **Vienna ab initio simulation package (VASP)**. This development makes it possible to investigate anharmonic and quantum nuclear effects on a solid’s properties using a wide variety of exchange-correlation functionals, which can optionally be employed to train quantum-accurate MTPs enabling SSCHA calculations at a low computational cost. The reference work associated with this library can be found at https://doi.org/10.1038/s41524-025-01553-1.

## Code Structure

The **utils_sscha_mlip** directory contains essential utilities for managing the interface between SSCHA, MLIP2, and DFT calculations:

- **Cluster_Management.py:** Handles input-output operations with the computing cluster, as well as the creation and management of inputs and outputs for self-consistent calculations.
- **MTP_environment.py:** Provides the interface between **SSCHA** and **MLIP version 2**, facilitating the active learning procedure.

## Scripts and Examples

The **scripts** folder contains example scripts demonstrating how to use the code:

- **SSCHA.py:** Implements the main workflow for the automated SSCHA + MLIP protocol.
- **potential.mtp:** Example of a level 10 Moment Tensor Potential, one of the various levels available in **MLIP version 2**.
- **utils.py:** Includes input scripts for first-principles calculations and cluster handling.
- **PdCuH_0.dyn:** Dynamical matrices for the PdCuH<sub>2</sub> example 

To correctly execute the example, place **SSCHA.py**, **potential.mtp**, the **PdCuH_0.dyn** files, the **Pseudo** folder, and the **utils.py** file in the same directory. Make sure to edit the utils.py file according to the specific needs of the cluster.

## Code Requirements

This repository requires the **SSCHA** and **mlip-2** codes to be installed on the machine where the code is executed, as well as a version of either **Quantum Espresso** or **VASP** on the machine used for self-consistent calculations.

- **SSCHA:** https://github.com/SSCHAcode
- **MLIP version 2:** https://gitlab.com/ashapeev/mlip-2
- **Quantum Espresso:** https://github.com/QEF
- **VASP:** https://www.vasp.at/

The code requires the following libraries in addition to those required by the **SSCHA**:
numpy, sys, os, time, math, contextlib

## Installation
- Make sure to have a version of **Wheel** installed.
- Enter the **SSCHA-MLIP** folder where the **setup.py** file is located.
- Execute the following commands:

       python setup.py sdist bdist_wheel

       cd dist

       pip install utils_sscha_mlip-0.3-py3-none-any.whl
