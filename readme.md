# Information Theoretic Measures of Vicsek Flocking Model

This is one of the companion projects to my PhD thesis: Information Theoretic Measures of Transitions to Collective Behaviour. It contains the code relating to the Vicsek model.

Simulation requires C++11 and OpenMP, and analysis scripts are provided for MATLAB. Runs on Linux and Windows (and Mac, but not thoroughly tested).

# Compiling

The code works on both Linux and Windows. Simply run cmake to configre and compile:

```
mkdir build; cd $_
cmake ..
cmake --build .
```

# Running simulations

Helper scripts are provided in `/tools/running_model` which will run either a single job (over a range of noise values) or a batch job, running `N` repetitions.

A single job will run `2n` instances, where `n` is the number of noise steps, of the program sequentially. It will first run the lead-in steps at the highest noise value to reach equilibrium, using `./vicsek-entropy generate`, and then run the analysis time window using `./vicsek-entropy analyse`. Once this is done, the noise vale is dropped and the process repeats.

Use `./vicsek-entropy -h` for more information about all available components, or `./vicsek-entropy <component> -h` for help with a specific component.

# Visualisation

Visualisation can be performed in MATLAB using the scripts in `/tools/processing_results`. `domatlab.sh` is the main script, which will automate the process. Before processing however, run the `/tools/processing_results/extract.sh` script on the directory to extract key data from the `*.log` files (placed in corresponding `*.log.I` files).

To run manually, there is a three stage process:

1. Set the results base directory: `setenv DATADIR '/path/to/vicsek/results';`
2. Gather the data files: `[filenames, eta, runs] = get_files('path');`
3. Process and visualise: `display_data(filenames, eta, runs);`
 
