# this module could easily be repurposed outside of bacteria/simulation/light
# the only parts specific to 'light' is the callback_dict argument and update_fitness()

# note that it assumes that the target function is MINIMIZED

# a large improvement would be using all numpy matrix multiplication
import sys, inspect
from pathlib import Path
file_path = Path(inspect.getfile(inspect.currentframe())).resolve()

LIGHT_DIR = file_path.absolute().parents[1]
sys.path.append(str(LIGHT_DIR))

from copy import deepcopy
import numpy as np, os, pickle, math, random as rd
import optimize

# unimplemented: emp_uncert: will add to convert (using stats.py) to find avg & var, then write them to seperate files


def find_opt(params, orig_params, callback_dict, target_data, population=None):
	print("Starting PSO")
	# note that params refers to opt_params
	from types import SimpleNamespace
	bounds = callback_dict['bounds']

	optimum_dict = {'x':[0 for i in range(len(params["target_params"]))], 'fun':np.inf,'message':'Done.'}
	# just making it consistent with other scipy output
	# 'x' is array of param values, 'fun' is the error, 'message' is the exit message

	if population is not None: #i.e. continuing run from a pickle file
		print("Continuing PSO optimization from pickled file.")

	else:
		population = init_pop(params,orig_params,target_data,callback_dict,optimum_dict)
		# population.keys() = pos, vel, fitness, lbest, lbest_fitness, gbest, gbest_fitness, 
		# each are matrices, except for fitness, lbest_fitness and gbest_fitness which are vectors
		# note that if the topology is not a star, then gbest needs to be a matrix instead
		# and need a matrix for the topology


	for i in range(params["max_steps"]):
		print("PSO at iter",i)
		update_pos(params, population, callback_dict)
		dset_indxs = update_fitness(params,orig_params,population,target_data,callback_dict)
		
		update_best(params, population,callback_dict, dset_indxs)
		update_vel(params, population)

		end_opt = record_data(params,population,callback_dict,optimum_dict, orig_params)

		if end_opt:
			break

	if i==params['max_steps']-1:
		optimum_dict['message'] = 'Optimization exit due to max iterations.'

	min_index = np.argmin(population['gbest_fitness'])
	optimum_dict['fun'] = population['gbest_fitness'][min_index]
	assert(optimum_dict['fun']==min(population['gbest_fitness']))
	optimum_dict['x'] = deepcopy(population['gbest'][min_index])

	optimum = SimpleNamespace(**optimum_dict)
	print("PSO finished.")
	return optimum


##################################### HELPER FUNCTIONS ############################################

def init_pop(params,orig_params,target_data,callback_dict,optimum_dict):
	print("Init PSO")
	# mostly copied from optimize.default_step()

	n,k = params['population_size'], len(params['target_params'])
	bounds = callback_dict['bounds']
	callback_dict['% improved'], callback_dict['avg pos dist'], callback_dict['avg at boundary'] = [],[],[]
	if params['use_temperature']:
		callback_dict['% backwards steps'] = []

	callback_dict['population_fitness'] = [[] for i in range(n)]
	if params["use_dataset_weights"]:
		callback_dict['dataset_weights'] = [1 for i in range(len(params['dataset']))]
	if params['one_datapoint_at_a_time']:
		callback_dict['dataset_freq'] = []

	if params['log_uniform_steps'] in [1,10]: 
		pos = np.array([np.power(10,(np.ma.log10(bounds[:, 1]) - np.ma.log10(bounds[:, 0])) \
			* (np.random.random(k)) + np.ma.log10(bounds[:, 0])) for i in range(n)])
	elif params['log_uniform_steps']==2: 
		pos = np.array([np.power(2,(np.ma.log2(bounds[:, 1]) - np.ma.log2(bounds[:, 0])) \
			* (np.random.random(k)) + np.ma.log2(bounds[:, 0])) for i in range(n)])
	else:
		pos = np.array([(bounds[:, 1] - bounds[:, 0]) * np.random.random(k) + bounds[:, 0] for i in range(n)])
	
	lbest = deepcopy(pos)
	gbest = deepcopy(pos)
	lbest_fitness = [np.inf for i in range(n)]
	gbest_fitness = [np.inf for i in range(n)]
	fitness = [np.inf for i in range(n)]
	# TODO: someday add log uniform vel init?
	vel_range = bounds[:, 1] - bounds[:, 0]
	vel = np.array([2*vel_range*np.random.random(k) - vel_range for i in range(n)])*params['init_vel_percent']

	population = {'pos':pos, 'vel':vel, 'lbest':lbest, 'lbest_fitness':lbest_fitness, \
		'gbest':gbest, 'gbest_fitness':gbest_fitness, 'fitness':fitness}

	dset_indxs = update_fitness(params,orig_params,population,target_data,callback_dict)
	update_best(params,population,callback_dict, dset_indxs,init=True)

	if min(population['gbest_fitness']) < params['min_loss']:
		optimum_dict['message'] = 'WARNING: Optimization exit due to sufficient fitness FROM THE VERY START!'
		assert(False)

	optimize.console_log(params,'\nPSO init population: ' + str(population),verbose_lvl=2)

	lower_bound_matrix = np.array([bounds[:,0] for i in range(n)])
	upper_bound_matrix = np.array([bounds[:,1] for i in range(n)])
	callback_dict['lower_bound_matrix'] = lower_bound_matrix
	callback_dict['upper_bound_matrix'] = upper_bound_matrix

	return population


def update_pos(params, population, callback_dict):
	lower_bound_matrix, upper_bound_matrix = callback_dict['lower_bound_matrix'], callback_dict['upper_bound_matrix']
	#print('\nPSO pos before update_pos=',population['pos'], '\nvel=',population['vel'])
	population['pos'] = population['pos'] + population['vel']

	#print('\nPSO pos after +vel:',population['pos'])
	population['pos'] = np.maximum(population['pos'],lower_bound_matrix)

	population['pos'] = np.minimum(population['pos'],upper_bound_matrix)
	#print('\nPSO pos after bounding:',population['pos'])


def update_fitness(params,orig_params,population,target_data,callback_dict):

	n = params['population_size']
	pos, fitness = population['pos'], population['fitness']
	if not params['one_datapoint_at_a_time']:
		dset_indxs=[None for i in range(n)] 
	else:
		if params['use_dataset_weights']:
			weights=callback_dict['dataset_weights']
		else:
			weights=None
		dset_indxs = rd.choices([i for i in range(len(params['dataset']))], k=n, weights=weights )

	#print("PSO: using dataset indices:",dset_indxs,', from dset weights:',weights)
	for i in range(n):

		if params['toy_problem']:
			fitness[i] = toy_problem(pos[i], params, orig_params, target_data, callback_dict, dset_indxs[i])
		else:
			fitness[i] = optimize.sim_and_calc_loss(pos[i], params, orig_params, target_data, callback_dict, dset_indxs[i])
		#best_fitness=min(fitness[i],best_fitness)
		callback_dict['population_fitness'][i] += [fitness[i]]
	return dset_indxs

def update_best(params, population, callback_dict, dset_indxs, init=False):
	# assumes target function is MINIMIZED
	n = params['population_size']
	num_replaced = 0
	num_poss_back, num_back = 0,0 #for seeing how temp does

	for i in range(params['population_size']):
		replaced = False
		if not params['use_temperature']:
			if population['fitness'][i] < population['lbest_fitness'][i]:
				replaced=True
		else:
			rel_temp = params['temperature']*population['lbest_fitness'][i]
			# temp based on scipy's implmtn in basinhopping, which uses Metropolis criterion
			#print('PSO: prev_best=%s, now_fitness=%s, rd must be below %s' %(population['lbest_fitness'][i],population['fitness'][i],math.exp((population['lbest_fitness'][i] - population['fitness'][i])/params['temperature'])))
			if population['lbest_fitness'][i] >= population['fitness'][i]:
				replaced=True
			elif rd.random() < math.exp((population['lbest_fitness'][i] - population['fitness'][i])/rel_temp):
				replaced=True
				num_back+=1
				num_poss_back+=1
			else:
				num_poss_back+=1

		if replaced:
			num_replaced += 1
			population['lbest'][i] = deepcopy(population['pos'][i])
			population['lbest_fitness'][i] = population['fitness'][i]
			if params['use_dataset_weights'] and params["one_datapoint_at_a_time"]:
				#print("PSO reweighting dataset %s, indiv #%s before was %s" %(dset_indxs[i],i,callback_dict['dataset_weights'][dset_indxs[i]]))
				callback_dict['dataset_weights'][dset_indxs[i]]*= params['dataset_reweight'] 
				#print("after was",callback_dict['dataset_weights'][dset_indxs[i]])
	if not init:
		callback_dict['% improved'] += [num_replaced/n]
		if params['use_temperature']:
			callback_dict['% backwards steps'] += [num_back/max(num_poss_back,1)]

	for i in range(params['population_size']):
		if params['topology'] == 'star':
			#if min(population['lbest_fitness']) < population['gbest_fitness'][i]:
			
			# note that gbest is nondeterministic, and not strictly decreasing, if min(lbest) isn't
			indx = np.argmin(population['lbest_fitness'])
			population['gbest'] = np.array([deepcopy(population['pos'][indx]) for j in range(n)])
			population['gbest_fitness'] = [population['lbest_fitness'][indx] for j in range(n)]
	
		elif params['topology'] == 'ring':
			nghs_n_self = [population['lbest_fitness'][(i-1)%n],population['lbest_fitness'][i],population['lbest_fitness'][(i+1)%n]]
			indx = (i-1+np.argmin(nghs_n_self))%n #such that index is in terms of the full population array, not just ngh_n_self

			population['gbest'][i] = deepcopy(population['lbest'][indx])
			population['gbest_fitness'][i] = population['lbest_fitness'][indx]

		else:
			print("\nERROR: PSO unknown topology param:",params['topology'])
			assert(False)


def update_vel(params, population):
	# is R matrix is supposed to be multiplied element wise, hence * instead of np.multiply()
	lbest_rd, gbest_rd = np.random.random(population['lbest'].shape),np.random.random(population['gbest'].shape)
	phi = params['weight_best']
	if params['use_constriction_factor']:
		K = 2/abs(2-phi-math.pow(math.pow(phi,2)-4*phi,.5))
	else:
		K=1

	population['vel'] = K*(params['inertia']*population['vel'] + \
		params['weight_best']*(1-params['weight_gbest'])*lbest_rd*(population['lbest']-population['pos']) + \
		params['weight_best']*params['weight_gbest']*gbest_rd*(population['gbest']-population['pos']))

def toy_problem(pos, params, orig_params, target_data, callback_dict, dset_indx):

	if not params['one_datapoint_at_a_time']:
		# checked all of these:
		loss1 = abs(pos[0]*pos[1]-pos[2]*pos[3])
		#loss2 = abs(pos[0]+math.pow(pos[1]*pos[2],2)-math.log(pos[3]))
		#loss3 = abs(math.pow(2,pos[3])+pos[1]*pos[2]-pos[0])
		#loss4 = abs(math.pow(pos[1]-pos[3],2)/(pos[0]*pos[2]))
		#loss5 = max(abs(pos[1]*pos[3]-pos[2]*pos[0]),abs(pos[1]/pos[2]-pos[3]*pos[0]))
		#loss =  max(loss1,loss2,loss3,loss4,loss5)
		# typically this one alone is harder:
		loss = abs(max(pos[1]*pos[2]*pos[3]-pos[0],math.pow(pos[0],2)/pos[1]-pos[2]/pos[3]))
		
		#print('pos',pos,'->loss=',loss)
		
	else:
		if dset_indx == 0:
			loss = abs(math.pow(pos[0],2)-(pos[1]))
		else:
			loss = abs(math.pow(pos[1],2)-(pos[0]))

	
	return loss


def record_data(params,population,callback_dict,optimum_dict,orig_params):

	min_index = np.argmin(population['lbest_fitness'])
	callback_dict['loss_over_steps'] += [population['lbest_fitness'][min_index]]
	callback_dict['num_steps'] += 1
	callback_dict['min_loss'] = population['lbest_fitness'][min_index]
	callback_dict['best_tuned_param_values'] = population['lbest'][min_index]

	# more stats about the population distrib
	n,k = params['population_size'],len(params['target_params'])
	pos_dist = 0
	for i in range(n):
		for other in range(n):
		#other = (rd.choice([i for i in range(n)])+i)%n
			for j in range(k):
				pos_dist += abs(sum(population['pos'][i]-population['pos'][other]))/(params['target_params'][j]['bounds'][1])
			pos_dist /= k
			#pos_dist += np.linalg.norm(population['pos'][i]-population['pos'][other],ord=1)
	callback_dict['avg pos dist'] += [pos_dist/math.pow(n,2)]

	at_boundary = 0
	for i in range(n):
		for j in range(k):
			if population['pos'][i,j] in params['target_params'][j]["bounds"]:
				at_boundary += 1

	callback_dict['avg at boundary'] += [at_boundary/(n*k)]

	if callback_dict['num_steps'] % (1/params['console_freq']) == 0:
		txt = "Optimizer at step %s with an error of %s\ncurr best params are %s." %(callback_dict['num_steps'], callback_dict['min_loss'],callback_dict['best_tuned_param_values'])
		optimize.console_log(params,txt,verbose_lvl=1)

	if callback_dict['num_steps'] % (1/params['pickle_freq']) == 0:
		pickle_file = os.path.join(params['opt_output_dir'],'step' + str(callback_dict['num_steps']) + '.pickle')
		with open(pickle_file, 'wb') as f:
			data = {'callback_dict':callback_dict,'opt_params':params, 'orig_params':orig_params,'population':population}
			pickle.dump(data,f)

	if population['lbest_fitness'][min_index] < params['min_loss']:
		# could add other reasons to exit optimization early here
		optimum_dict['message'] = 'Optimization exit due to sufficient fitness.'
		return True
	else:
		return False