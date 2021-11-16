import sys, os, inspect
from pathlib import Path

local_dir = Path(inspect.getfile(inspect.currentframe()))
LIGHT_DIR = local_dir.absolute().parents[1]
sys.path.append(str(LIGHT_DIR))


IBIOSIM_PROPS = """reb2sac.diffusion.stoichiometry.amplification.value=1.0
ode.simulation.time.step={sim_stepsize}
ode.simulation.absolute.error={absolute_tolerance}
monte.carlo.simulation.runs=1
reb2sac.qssa.condition.1=0.1
reb2sac.abstraction.method=none
monte.carlo.simulation.start.index=1
reb2sac.operator.max.concentration.threshold=15
ode.simulation.print.interval={print_stepsize}
ode.simulation.out.dir=OUTPUT_DIR
simulation.initial.time=0.0
simulation.output.start.time=0.0
ode.simulation.time.limit={duration}
ode.simulation.relative.error={relative_tolerance}
simulation.printer.tracking.quantity=amount
ode.simulation.min.time.step=0
simulation.run.termination.decider=constraint
simulation.printer=tsd.printer
monte.carlo.simulation.random.seed=314159
reb2sac.rapid.equilibrium.condition.2=0.1
reb2sac.simulation.method=ODE
reb2sac.rapid.equilibrium.condition.1=0.1
reb2sac.generate.statistics=false"""


def ibiosim(params):
	#for key in ['duration','print_stepsize','sim_stepsize','absolute_tolerance','relative_tolerance']:
	#	IBIOSIM_PROPS = IBIOSIM_PROPS.replace(key.upper(),params[key])

	props = IBIOSIM_PROPS.replace('OUTPUT_DIR',os.path.join(LIGHT_DIR,params['output_dir']))
	with open(os.path.join(LIGHT_DIR,'cps/ibiosim_props.txt'), 'w') as f:
		f.write(props.format(**params))
