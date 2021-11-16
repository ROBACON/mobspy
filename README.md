
# Light Copasi Simulations using Python

This package allows for streamlined simulations using COPASI solvers, including merging many stochastic simulations and optimization to a target file. 

Organization:
1. Requirements
2. Directories
3. Basic Usage
4. Parameters
5. Plotting
6. Stochastic Runs
7. Data Format
8. Building a Model
9. Optimization
10. Experiments
11. Multi Runs
12. Issues and Future Work


## 1. Requirements

Copasi must be installed and added to the PATH environment variable, such that CopasiSE can be called via the command line from anywhere.

Python3 and the relevant python libraries must be installed. One can install all or most of the required python libraries via "pip3 install -r requirements.txt." Any remaining libraries must be installed manually using pip following the relevant error message.


## 2. Directories
- `params/`: parameter files
- `output/`: the default directory for data, plots, and optimization results
- `models/`: set of available models and a builder to translate them into Copasi
- `lib/`: utility scripts, suchs as statistics and parsing
- `experiments/`: more involved batch runs for comparing model behaviors
- `docs/`: documentation for parameter options, file descriptions, and further details
- `data/`: data file used for optimization, typically for empirical experiments
- `cps/`: files generated to run Copasi, only open to manually debug or save Copasi versions


## 3. Basic Usage

Running a single simulation:
	`python3 run.py params/PARAM_FILE`

With the default parameter file `params/AB.txt`, two output data files should appear in `output/` and a plot should pop up. 

Modify the chosen PARAM_FILE to rapidly pick a model and specify simulation parameters. For instance, one can change the algorithm used for simulation, initial concentrations of species, or kinetic rates of the model. 

See commands for plotting, experiments, optimization, and multi runs in their respective sections.


## 4. Parameters
There are three types of parameter files: basic, plotting, and optimization. Examples of existing parameter files can be found in `params/`, where the simplest basic, plotting and optimization files are `params/AB.txt`, `params/EPAplot.json`, and `params/EPAopt.json`, respectively. Further explanations of the format and parameters of all three can be found in their respective files in `doc/`, although they may not be up to date. 

In brief, the basic parameter file specifies a single run. It dictates the model, the parameter values, and the simulation algorithm. It must be tab-delimited for the parser to work. The plotting parameter file is optional and customizes plotting output. It uses a json format. Finally, the optimizer parameter file specifies which parameters are optimized and how to do so. It also uses a json format. Note that optimizer parameters reference a basic parameter file to specify the model to optimize and the values of parameters which are not optimized. 


## 5. Plotting
Plotting during a single run can be toggled with the "plot" option in the basic parameters. Without specifying seperate plotting parameters, default plots will generated. To further customize these plots, specify a file in "plot_json" parameter. This will be used to load plotting parameters, which can alter default plots or specify custom plots. More details can be found in `docs/plot_params.txt`.

Plotting is entirely seperated from simulation such that plots can be generated during or after simulation with ease. So long as the "save_data" parameter is ON, then plots can be generated from the resulting pickled data as follows: 
	`python3 plot.py PICKLED_FILE`

This will use the parameters of the original run and the plotting parameter file of the original run. Note that the plotting parameter file can be altered after the initial simulation, but before plotting the data. If a different plotting parameter file is desired, one can instead use:
	`python3 plot.py PICKLED_FILE PLOTTING_PARAM_FILE`

Note that during batch runs such as experiments or optimization, plotting and saving data for each individual run is already disabled through a function argument.


## 6. Stochastic Runs
Stochastic simulation is specified by changing the basic parameter simulation_method. The default option of "stochastic" uses the Direct Method of Gillespie. It is the most realistic, but also the slowest. Other options include 'TauLeap', which is described in the COPASI reference manual. Unclear how to use 'NextReaction', 'AdaptiveSSA' methods mentioned in COPASI manual. 

`WARNING`: using the wrong method is silently exchanged for deterministic by COPASI.

All output and plotting automatically handles stochastic and deterministic runs. Note that in the case of a deterministic run, the original COPASI tsd output and a pickle are saved. In the case of a stochastic run, only a pickled file with the merged data of many COPASI runs is saved.

Stochastic simulations are run multiple times to calculate sample mean and variance. Note that several more parameters are required for stochastic simulation. Confidence intervals about the sample mean and variance are calculated using a t-test with < 30 reptitions, and a normal distribution otherwise. Skewed confidence intervals can be used if the confidence intervals are suspected to be asymmetric.


## 7. Data Format
Output of a single run is structured such that stochastic and deterministic runs have a similar format. In this way, one can change the type of run without having to manually code their difference.

Stochastic runs contain more information on the average and variance of each species, as well as confidence intervals about the average and variance. As such each simulation data is a dictionary of species, each of which have {'avg':{'val':[...],'CI_min':[...],'CI_max':[...]},'var':{'val':[...],'CI_min':[...],'CI_max':[...]}}. The only "species" that is an exception is "Time":[], since time has no average or variance across different runs. See `lib/stat.py` where data for multiple stochastic runs are merged for details.

In contrast, deterministic runs contain only an average. The leads to the oddity species:{'avg':{'val':[...]}} for each species in the dictionary. One can avoid reformatting deterministic output, although it is generally not advised and may break other parts. To do so, set reformat=False in run.simulate_one().


## 8. Building a Model
One can follow the format of `models/AB.py`. A seperate .py file should be created that builds a dictionary of reactions, species, events, and parameters with the proper format for each. New models will be automatically detected if the model parameter matches the python file name.

Mappings are optional and used to reorganize groups of species for plotting and output data. Unmapped species will automatically be kept as their own group. Mappings can also be used to rename species between coding and plotting.


## 9. Optimization
Parameters can be tuned to match a given dataset. The loss function to minimize is based on a specified distance metric between each point of the dataset and the simulation. Error is normalized such that it is the average distance between a simulation point and an empirical datapoint. Currently the optimization uses basinhopping (implemented by Scipy) or particle swarm optimization (PSO, manually implemented).

Optimization relies on an optimizer parameter file, which also includes a reference to a basic parameter file. Any parameters defined in a basic parameter file can be tuned, and all non-tuned parameters are set to the values in the basic parameter file. See `docs/opt_params.txt` and `docs/basic_params.txt` for more info. 

To optimize a model, one uses:
	`python3 optimize.py params/OPT_PARAM_FILE`

An optimization run will place all of its output into one directory. Among other outputs, it will periodically pickle its progress. Doing so allows one to continue in the event of some failure. Continuing an optimization from a pickle is as follows:
	`python3 optimize.py PICKLED_FILE --pickle`
The final plots from an optimization run can be easily re-plotted using one of the output pickle files as follows:
	`python3 plot.py PICKLE_FILE --opt`


## 10. Experiments
WARNING: many experiments may be defunct due to large changes in stochastic runs and parameters files.

Experiments are custom scripts to compare particular parameter ranges or models. They are called directly:
	`python3 experiments/NAME.py params/PARAM.txt`

Note that their usage may vary by experiment since they are custom scripts. They are also often managed directly in the code in addition to the usual parameter files.


## 11. Multi Runs
For running multiple simulations to generate a single image, usually to compare to multiple time series of data:
	`python3 run.py params/PARAM_FILE --multi` 

This will scan the plot json referenced in the param file for any plots with the `multi:1` field (will do nothing if there is no plot json or no such plots). For such a plot, for each entry in the `dataset` field, it will run one simulation. Within `dataset`: 'simulation' is the name of the species (or mapping) to plot, 'empirical' is the column header of the data to plot, 'init_params' will be set to the value in 'init_vals' (overrides txt param file). If 'init_val' is 'T0' or 'Time0', it will use the first row of the 'empirical' column. See EPAplot.json for more examples.

Data from a multi run is also saved as a pickle, such that it can be easily re-plotted using:
	`python3 plot.py PICKLE_FILE --multi`


## 12. Issues and Future Directions
- errors from the copasi solver are not handled gracefully
- currently one has to make sure that a simulation has the same number of time points as the target data
- interpolation is not currently implemented
- weighting distance during optimization by confidence intervals or variance is not currently implemented
- if one uses the wrong simulation_method parameter, it is silently exchanged for deterministic by COPASI. It should be explicitly caught here instead
- optimization should really use training and validation sets
- optimization assumes floats, but this should be easy to change