#!/usr/bin/env python3

import sys
import json
import re
import codecs
import argparse

import sys, inspect, os
from pathlib import Path


def cps_set_stochmethod(params, method_to="stochastic"):
	# rewrite and create model_new.cps
	ofile = codecs.open(params['cps_dir']+'model_new.cps', "w", "utf-8")
	with codecs.open(params['cps_dir']+'model.cps', "r", "utf-8") as fc:
		for line in fc:
			if f'<Method name="Stochastic"' in line:
				print(f'<Method name="Stochastic" type="{method_to}">', end="", file=ofile)
	ofile.close()

	# finally move new file to replace old file
	os.rename(params['cps_dir']+'model_new.cps', params['cps_dir']+'model.cps')



def light_version(params, species):

	# parse xml
	prog1 = re.compile(".*<ListOfTasks.*")
	xml_task = f"""<Task key="Task_999" name="Time-Course" type="timeCourse" scheduled="true" updateModel="false">
	<Report append="0" confirmOverwrite="0" reference="Report_999" target="{params['output_file']}"/>
		<Problem>
			<Parameter name="AutomaticStepSize" type="bool" value="0"/>
			<Parameter name="StepNumber" type="unsignedInteger" value="{int(params['duration']/params['sim_stepsize'])}"/>
			<Parameter name="StepSize" type="float" value="{params['sim_stepsize']}"/>
			<Parameter name="Duration" type="float" value="{params['duration']}"/>
			<Parameter name="TimeSeriesRequested" type="bool" value="1"/>
			<Parameter name="OutputStartTime" type="float" value="0"/>
			<Parameter name="Output Event" type="bool" value="0"/>
			<Parameter name="Start in Steady State" type="bool" value="0"/>
		</Problem>
	"""

	if params['simulation_method'] in ['LSODA','deterministic','Deterministic']:
		xml_task+= f"""
			<Method name="Deterministic (LSODA)" type="Deterministic(LSODA)">
				<Parameter name="Integrate Reduced Model" type="bool" value="0"/>
				<Parameter name="Relative Tolerance" type="unsignedFloat" value="{params['relative_tolerance']}"/>
				<Parameter name="Absolute Tolerance" type="unsignedFloat" value="{params['absolute_tolerance']}"/>
				<Parameter name="Max Internal Steps" type="unsignedInteger" value="{params['max_internal_steps']}"/>
				<Parameter name="Max Internal Step Size" type="unsignedFloat" value="{params['max_internal_step_size']}"/>
			</Method>
			</Task>
			"""
	elif params['simulation_method'] in ['Stochastic','stochastic','DirectMethod', 'NextReactionMethod','TauLeap', 'AdaptiveSSA']:
		if params['simulation_method'] in ['Stochastic','stochastic']:
			# method = 'DirectMethod'
			method = 'Stochastic'
		else:
			method = params['simulation_method']
		xml_task+= f"""
			<Method name="Stochastic (Direct method)" type="{method}">
				<Parameter name="Max Internal Steps" type="integer" value="{params['max_internal_steps']}"/>
				<Parameter name="Use Random Seed" type="bool" value="0"/>
				<Parameter name="Random Seed" type="unsignedInteger" value="1"/>
			</Method>
			</Task>
			"""
		# xml_task+= f"""
		# 	<Method name="Stochastic" type="{method}">
		# 		<Parameter name="Max Internal Steps" type="integer" value="{params['max_internal_steps']}"/>
		# 		<Parameter name="Subtype" type="unsignedInteger" value="2"/>
		# 		<Parameter name="Use Random Seed" type="bool" value="0"/>
		# 		<Parameter name="Random Seed" type="unsignedInteger" value="1"/>
		# 	</Method>
		# 	</Task>
		# 	"""
	else:
		print("\nERROR: unknown parameter for simulation_method!!\n", file=sys.stderr)
		assert(False)

	prog2 = re.compile(".*<ListOfReports.*")
	xml_plot_head = """
	<Report key="Report_999" name="Time-Course" precision="6" separator="&#9;" taskType="Time-Course">
	<Comment/>
		<Table printTitle="1">
		<Object cn="CN=Root,Model=NoName,Reference=Time"/>
	"""
	xml_plot_foot = """
		</Table>
	</Report>
	"""
	xml_lines = ""
	for s in species:
		xml_lines += (
			'<Object cn="CN=Root,Model=NoName,Vector=Compartments[c1],Vector=Metabolites[%s],Reference=Concentration"/>\n'
			% (s,)
		)

	# these file names matter and must match those in run.py
	ofile = codecs.open(params['cps_dir']+'model.cps', "w", "utf-8")
	with codecs.open(params['cps_dir']+'pre.cps', "r", "utf-8") as fc:
		for line in fc:
			if prog1.match(line):
				print(line + xml_task, end="", file=ofile)
			elif prog2.match(line):
				print(
					line + xml_plot_head + xml_lines + xml_plot_foot, end="", file=ofile
				)
			else:
				print(line, end="", file=ofile)
	ofile.close() 




def module_version(
		cps, react_json, settings_json, simulation, settings_as_py=False, write_cps_to=None
):
		# get species
		# with open(react_json, 'r') as infile:
		with codecs.open(react_json, "r", "utf-8") as infile:
				obj = json.load(infile)
		species = obj["species"].keys()

		# get settings
		if settings_as_py:
				settings = settings_json
		else:
				with codecs.open(settings_json, "r", "utf-8") as infile:
						settings = json.load(infile)["settings"]

		# parse xml
		prog1 = re.compile(".*<ListOfTasks.*")
		xml_task = f"""<Task key="Task_999" name="Time-Course" type="timeCourse" scheduled="true" updateModel="false">
	<Report append="0" confirmOverwrite="0" reference="Report_999" target="{simulation}"/>
				<Problem>
					<Parameter name="AutomaticStepSize" type="bool" value="0"/>
					<Parameter name="StepNumber" type="unsignedInteger" value="{int(settings['duration']/settings['stepsize'])}"/>
					<Parameter name="StepSize" type="float" value="{settings['stepsize']}"/>
					<Parameter name="Duration" type="float" value="{settings['duration']}"/>
					<Parameter name="TimeSeriesRequested" type="bool" value="1"/>
					<Parameter name="OutputStartTime" type="float" value="0"/>
					<Parameter name="Output Event" type="bool" value="0"/>
					<Parameter name="Start in Steady State" type="bool" value="0"/>
				</Problem>
				<Method name="Deterministic (LSODA)" type="Deterministic(LSODA)">
					<Parameter name="Integrate Reduced Model" type="bool" value="0"/>
					<Parameter name="Relative Tolerance" type="unsignedFloat" value="{settings['relative_tolerance']}"/>
					<Parameter name="Absolute Tolerance" type="unsignedFloat" value="{settings['absolute_tolerance']}"/>
					<Parameter name="Max Internal Steps" type="unsignedInteger" value="100000"/>
					<Parameter name="Max Internal Step Size" type="unsignedFloat" value="0"/>
				</Method>
	</Task>
	"""

		prog2 = re.compile(".*<ListOfReports.*")
		xml_plot_head = """
	<Report key="Report_999" name="Time-Course" precision="6" separator="&#9;" taskType="Time-Course">
		<Comment/>
			<Table printTitle="1">
				<Object cn="CN=Root,Model=NoName,Reference=Time"/>
	"""
		xml_plot_foot = """
			</Table>
	</Report>
	"""
		xml_lines = ""
		for s in species:
				xml_lines += (
						'<Object cn="CN=Root,Model=NoName,Vector=Compartments[c1],Vector=Metabolites[%s],Reference=Concentration"/>\n'
						% (s,)
				)

		ofile = codecs.open(write_cps_to, "w", "utf-8") if write_cps_to else sys.stdout
		with codecs.open(cps, "r", "utf-8") as fc:
				for line in fc:
						if prog1.match(line):
								print(line + xml_task, end="", file=ofile)
						elif prog2.match(line):
								print(
										line + xml_plot_head + xml_lines + xml_plot_foot, end="", file=ofile
								)
						else:
								print(line, end="", file=ofile)
		ofile.close()


if __name__ == "__main__":
		# from command line calls the module (heavyweight) version
		parser = argparse.ArgumentParser()
		parser.add_argument("cps", help="the cps file to rewrite")
		parser.add_argument("react_json", help="the model's reactions json file")
		parser.add_argument("settings_json", help="the settigs json file")
		parser.add_argument(
				"output", help="the name of the file to write simulation results to"
		)
		args = parser.parse_args()

		module_version(args.cps, args.react_json, args.settings_json, args.output, None)
