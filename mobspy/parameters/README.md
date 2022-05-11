# Parameter Explanation

## Model paramters

- `volume` (float w/wo unit): sets the volume in the simulation (default is 1 litre)

- `repetitions` (int): sets the number of runs for stochastic
        
- `level` (int): sets the logging level for the simulator

## BasiCO parameters 

- `simulation_method`: sets the simulation method to use (deterministic, stochastic)

- `duration` (float w/wo unit): the duration in time units for how long to simulate)

- `start_time` (float): the output start time. If the model is not at that start time, a simulation
will be performed in one step, to reach it before starting to collect output.

- `step_size` (float): step_size of the output in simulation

- `seeds` (list): set the seeds that will be used 

- `a_tol` (float): the absolute tolerance to be used

- `r_tol` (float): the relative tolerance to be used

- `output_event` (bool): if true, output will be collected at the time a discrete event occurs.

## Output parameters

- `unit_x` (str): output unit for the x axis

- `unit_y` (str): output unit for the y axis

- `output_concentration` (bool): output in concentration or counts

- `output_dir` (str): directory to save the parameters 

- `output_event` (bool): if true, output will be collected at the time a discrete event occurs

- `output_file` (str): name of the output file

- `save_data` (bool): save the data in a JSON file or not

- `plot_data` (bool): perform standard plotting after simulation


