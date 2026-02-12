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

If you want to develop the package, you can use an editable installation:
```
pip install -e .
```

### Use lisabeta

If you want to use [`lisabeta`](https://gitlab.in2p3.fr/marsat/lisabeta) for SNR computation, use:

```bash
pip install '.[lisabeta_snr]'
```

Note that at this point, installation of `lisabeta` is a bit tricky.

## Usage

### Getting started

The package comes with some example configuration files in folder `example_config_files/`.
To run some example, follow the steps below:

1. Generate LVK-like synthetic data
   - Run the configuration file with:
     ```bash
     iv_generate_events --config-file example_config_files/config_population_LVK_BBH.ini
     ```
   * This generates a population of $\sim 50$ detected events in a few seconds.

2. Estimate the detectors’ sensitivity
   - Run the configuration file with:
     ```bash
     iv_generate_events --config-file example_config_files/config_injections_LVK_BBH.ini
     ```
   * This generates a set of $\sim 10^3$ detected events in $\sim 2$ minutes.

3. Run the hierarchical inference with `icarogw`
   - Update the population and injections paths in the configuration file (using the outputs from the previous steps).
   - Run the configuration file with:
     ```bash
     iv_hierarchical_inference --config-file example_config_files/config_inference_LVK_BBH.ini
     ```
   * This step generates posterior samples for the population parameters in $\sim 5$ minutes, along with automatic diagnostic plots.

### Optional: Fully realistic simulations

To make the simulation fully realistic, you can provide `icarogw` with parameter estimation results for individual events.
   - Generate the parameter estimation config files for each event of your simulated population with:
     ```bash
     iv_generate_pe_configs -d output_directory_of_step_1
     ```
   - Run parameter estimation for each eevnt with the previously generated config files (files names are given as examples here, they should in practice depend of the name of `output_directory_of_step_1`). Note that this step has not been callibrated for example use, so running PE on the events given in example may take hours.
     ```bash
     iv_parameter_estimation --config output_directory_of_step_1/parameter_estimation/config_file_PE_event_0000.ini
     iv_parameter_estimation --config output_directory_of_step_1/parameter_estimation/config_file_PE_event_0001.ini
     ...
     iv_parameter_estimation --config output_directory_of_step_1/parameter_estimation/config_file_PE_event_0050.ini
     ```
     A calibrated example for parameter estimation which runs in less than $\sim 10$ min is available by running:
     ```bash
     iv_parameter_estimation --config config_example_PE.ini
     ```
   - Combine the results from as many populations you generated in a single files:
     ```bash
     iv_combine_pe_posteriors -d output_directory_1_of_step_1 ... output_directory_N_of_step_1 -o where_to_store_the_combined_files
     ```
     Note that in most cases, you will only pass one `output_directory_1_of_step_1`-like argument to the -d option.
   * This generates posterior samples for the observed events using `bilby`.
   - Feed these posterior samples into Step 3 and set `true-values = 0` in the configuration file.

You can change the population model and population parameters directly in the configuration file.


## Documentation :construction:

The documentation is under construction. Please be patient!


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


## License and Citation

`icaroverse` is distributed under the EUROPEAN UNION PUBLIC LICENCE v. 1.2. See the [LICENSE file](./LICENSE).

If you use `icaroverse` in your research or publications, please cite the package as follows:

```bibtex
@software{icaroverse,
  author       = {Gennari, Vasco and Bertheas, Tom and Dubois, Mathieu},
  title        = {icaroverse},
  year         = {2025},
  publisher    = {GitHub},
  url          = {https://github.com/vascogennari/icaroverse},
  version      = {v1.0.0},
}
```
