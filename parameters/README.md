# Parameter Explanation
For now, I have listed all the parameters basiCO allows us to use.
We must decide which ones we implement and give to the user

- update_model (bool) : sets whether the model should be updated, or reset to initial conditions

- use_initial_values (bool) : whether to use initial values

- update_model (bool) : sets whether the model should be updated, or reset to initial conditions

- simulation_method : sets the simulation method to use (deterministic, stochastic)

- duration (float): the duration in time units for how long to simulate)

- automatic (bool): whether to use automatic determined steps (True), or the specified interval / number of steps

- start_time (float): the output start time. If the model is not at that start time, a simulation
will be performed in one step, to reach it before starting to collect output.

- step_size (int): the output step size

- seeds(list): set the seeds that will be used 

- use_seed(bool): if true, the specified seed will be used

- a_tol(float) : the absolute tolerance to be used

- r_tol(float) : the relative tolerance to be used

- max_steps(int) : the maximum number of internal steps the integrator is allowed to used


