# icaroverse

**icaroverse** is a flexible and extensible framework for automating **gravitational-wave population and cosmological inference** with [`icarogw`](https://github.com/simone-mastrogiovanni/icarogw).

It provides an end-to-end workflow to:  
- Initialize the **population model** and hierarchical likelihood  
- Run Bayesian samplers for population and cosmology  
- Postprocess results with built-in diagnostics and visualizations  

In addition, dedicated submodules support:  
- **Synthetic catalog generation** and sensitivity estimates for selection effects  
- **Parameter estimation** on gravitational-wave events with [`bilby`](https://lscsoft.docs.ligo.org/bilby/)  


## Installation

To install the package, clone the git repository and then run:
```
pip install .
```

If you want to develop the package, you can do an editable installation:
```
pip install -e .
```

## Usage

### Getting started

The package comes with some example configuration files in folder `example_config_files/`.
To run some example, follow the steps below:

1. Generate LVK-like synthetic data
   - Run the configuration file with:
     ```bash
     iv_generate_events --config example_config_files/config_population_Powerlaw-Gaussian.ini
     ```
   * This generates a population of $\sim 50$ detected events in a few seconds.

2. Estimate the detectors’ sensitivity
   - Run the configuration file with:
     ```bash
     iv_generate_events --config-file example_config_files/config_injections_Powerlaw.ini
     ```
   * This generates a set of $\sim 10^3$ detected events in $\sim 2$ minutes.

3. Run the hierarchical inference with `icarogw`
   - Update the population and injections paths in the configuration file (using the outputs from the previous steps).
   - Run the configuration file with:
     ```bash
     iv_hierarchical_inference --config-file example_config_files/config_Powerlaw-Gaussian.ini
     ```
   * This step generates posterior samples for the population parameters in $\sim 5$ minutes, along with automatic diagnostic plots.

### Optional: Fully realistic simulations

To make the simulation fully realistic, you can provide `icarogw` with parameter estimation results for individual events.
   - Generate the parameter estimation config file with:
     ```bash
     iv_generate_pe_configs -d output_of_step_1
     ```
   - Update the path in the configuration file to point to the output of Step 1.
   - Run the configuration file with:
     ```bash
     iv_parameter_estimation --config config_file_PE.ini
     ```
   - Combine the results:
     ```bash
     iv_combine_pe_posteriors -d output_of_step_1 -o output_of_previous_step
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
