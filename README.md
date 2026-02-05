# icaroverse

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Build](https://img.shields.io/github/actions/workflow/status/vascogennari/icaroverse/ci.yml?branch=main)](https://github.com/vascogennari/icaroverse/actions)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](#-contributing)

**icaroverse** is a flexible and extensible framework for automating **gravitational-wave population and cosmological inference** with [`icarogw`](https://github.com/simone-mastrogiovanni/icarogw).

It provides an end-to-end workflow to:  
- Initialize the **population model** and hierarchical likelihood  
- Run Bayesian samplers for population and cosmology  
- Postprocess results with built-in diagnostics and visualizations  

In addition, dedicated submodules support:  
- **Synthetic catalog generation** and sensitivity estimates for selection effects  
- **Parameter estimation** on gravitational-wave events with [`bilby`](https://lscsoft.docs.ligo.org/bilby/)  


## Installation

_Missing_


## Usage

To get started with the package, follow the steps below.

1. Generate LVK-like synthetic data
   - Update the PSD path in the configuration file (it should point to the package directory).
   - Run the configuration file with:
     ```bash
     icaroverse --generate-events config_file_population.ini
     ```
   * This generates a population of ~50 detected events in a few seconds.

2. Estimate the detectors’ sensitivity
   - Update the PSD path in the configuration file (it should point to the package directory).
   - Run the configuration file with:
     ```bash
     icaroverse --generate-injections config_file_injections.ini
     ```
   * This generates a set of ~$10^3$ detected events in ~2 minutes.

3. Run the hierarchical inference with `icarogw`
   - Update the population and injections paths in the configuration file (using the outputs from the previous steps).
   - Run the configuration file with:
     ```bash
     icaroverse --runner config_file_inference.ini
     ```
   * This step generates posterior samples for the population parameters in ~5 minutes, along with automatic diagnostic plots.

### Optional: Fully realistic simulations

To make the simulation fully realistic, you can provide `icarogw` with parameter estimation results for individual events.
   - Update the path in the configuration file to point to the output of Step 1.
   - Run the configuration file with:
     ```bash
     icaroverse --parameter-estimation config_file_PE.ini
     ```
   * This generates posterior samples for the observed events using `bilby`.
   - Feed these posterior samples into Step 3 and set `true-values = 0` in the configuration file.

You can change the population model and population parameters directly in the configuration file.
All available options are documented here.


## Documentation

Detailed documentation on package usage can be found here.


## Contribute

Contributions are welcome and greatly appreciated! If you’d like to improve the code, fix bugs, or add new features, please follow these steps:

1. Fork the repository and create a new branch from main.

2. Make your changes with clear, well-documented commits.

3. Ensure your code follows the existing style and passes any tests (if applicable).

4. Open a Pull Request describing:

   - what you changed
   - why the change is needed
   - any relevant context or issues

We are happy to review pull requests and provide feedback. Feel free to open an issue first if you’d like to discuss an idea before implementing it.


## Citation

If you use `icaroverse` in your research or publications, please cite the package as follows:

```bibtex
@software{icaroverse,
  author       = {Your Name and Collaborators},
  title        = {icaroverse: A toolkit for ...},
  year         = {2025},
  publisher    = {GitHub},
  url          = {https://github.com/vascogennari/icaroverse},
  version      = {v1.0.0},
}
```
