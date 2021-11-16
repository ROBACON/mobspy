#!/usr/bin/env python3
import sys, inspect
from pathlib import Path
import argparse

file_path = Path(inspect.getfile(inspect.currentframe())).resolve()
sim_dir = file_path.parents[1]

LIGHT_DIR = file_path.absolute().parents[0]
sys.path.append(str(LIGHT_DIR))
sys.path.append(str(sim_dir))

import scipy.optimize as opti
import run
from time import time, strftime, localtime
from utilities import parse, PSO
import codecs, json, os, shutil, pickle
import numpy as np, math, random as rd
from copy import deepcopy



def full_run(param_file):

    print("\nStarting optimization with param file: %s ...\n" %(param_file))
    with open(param_file,'r') as file:
        opt_params = json.load(file)
    orig_params, target_data, callback_dict, start_time, x0 = init_optimization(opt_params)
    console_log(opt_params,"Beginning optimization with %s..." %(opt_params['global_opt_method']),verbose_lvl=1)
    continue_optimization(opt_params,  orig_params, target_data, callback_dict, start_time, x0)

def pickle_run(pickled_data):
    print("\nStarting optimization with pickled file: %s ...\n" %(pickle_file))
                
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)
        for k in ['opt_params','orig_params','callback_dict']:
            if k not in data.keys():
                print("\nERROR pickled data is missing key:\n",k)
                assert(False)

    if data['opt_params']["global_opt_method"] == 'PSO':
        population = data['population']
    else:
        population = None

    start_time = time()
    opt_params = data['opt_params']
    x0 = [opt_params['target_params'][i]['initial_value'] for i in range(len(opt_params['target_params']))]
    target_data = parse.read_csv_target_file(opt_params['target_file'], time_units=opt_params['target_data_time_units'], transpose=False)

    console_log(data['opt_params'],"Continuing scipy basinhopping optimization...",verbose_lvl=1)
    continue_optimization(data['opt_params'],  data['orig_params'], target_data, data['callback_dict'], start_time, x0, PSOpopulation=population)


##################################### PRIMARY FUNCTIONS ###########################################

def init_optimization(opt_params):
    orig_params = parse.params(opt_params['original_params'])
    if not init_check(opt_params, orig_params):
        assert(False)
    
    target_data = parse.read_csv_target_file(opt_params['target_file'], time_units=opt_params['target_data_time_units'], transpose=False)
    #if opt_params["target_uncert_file"]:
    #   callback_dict['target_uncertainty'] = parse.read_csv_target_file(opt_params['target_uncert_file'], time_units=opt_params['target_data_time_units'], transpose=False)
    
    start_time = time()

    x0 = [opt_params['target_params'][i]['initial_value'] for i in range(len(opt_params['target_params']))]

    callback_dict = {'min_loss':np.inf, 'loss_over_steps': [], 'num_steps':0, 'actual_num_runs':0}
    bounds = np.array([opt_params['target_params'][i]['bounds'] for i in range(len(opt_params['target_params']))],dtype=np.float64)
    callback_dict['bounds'] = bounds
    callback_dict['best_tuned_param_values'] = {opt_params['target_params'][i]['name']:x0[i] for i in range(len(x0))}

    for i in range(len(bounds)):
        if bounds[i][0] > bounds[i][1]:
            console_log(opt_params,"ERROR: bounds must be formatted as [MIN, MAX]\n", verbose_lvl=-1)
            assert(False)
        elif x0[i] < bounds[i][0] or x0[i] > bounds[i][1]:
            console_log(opt_params,"ERROR: initial parameter value is not in bounds\n", verbose_lvl=-1)
            assert(False)

    if opt_params['decay']:
        cutoff = int(round(orig_params['duration']/orig_params['sim_stepsize']))+1
        decay_vec = np.ones(cutoff)
        for i in range(math.ceil(cutoff/opt_params['decay_step'])):
            decay_vec[i*opt_params['decay_step']:]*=opt_params['decay_factor']  
        callback_dict['decay_vec'] = decay_vec * cutoff/np.linalg.norm(decay_vec,ord=1)

    # initialize output directory and some of its contents
    # only after other initialization to avoid output trash on faulty runs
    opt_params['opt_output_dir'] =  os.path.join(LIGHT_DIR,opt_params['opt_output_dir'].replace('TIMESTAMP',orig_params['timestamp']))
    if not os.path.exists(opt_params['opt_output_dir']):
        os.makedirs(opt_params['opt_output_dir'])
    with open(os.path.join(opt_params['opt_output_dir'],'log.txt'),'w') as file:
        file.write('OPTIMIZATION LOG\n')
    with open(os.path.join(opt_params['opt_output_dir'],'orig_params.txt'),'w') as file:
        for param in orig_params:
            file.write(param+'\t'+str(orig_params[param]) + '\n')
    with open(os.path.join(opt_params['opt_output_dir'],'opt_params.json'),'w') as file:
        json.dump(opt_params, file)

    orig_params['solver_log'] = os.path.join(opt_params['opt_output_dir'],'solver_log.txt')
    with open(orig_params['solver_log'],'w') as file:
        file.write('SOLVER LOG\n')

    console_log(opt_params,"Starting optimization, writing to %s ... " %(opt_params['opt_output_dir']),verbose_lvl=1)
    console_log(opt_params,"Initializing and running initial parameters...",verbose_lvl=1)

    if opt_params['global_opt_method'] != 'PSO':
        # this is a bit hacky, but PSO inits a rd population instead of one spc starting point
        init_loss = sim_and_calc_loss(x0, opt_params, orig_params, target_data, callback_dict)
        callback_dict['min_loss'] = init_loss
        callback_dict['loss_over_steps'] += [init_loss]

    return orig_params, target_data, callback_dict, start_time, x0


def continue_optimization(opt_params, orig_params, target_data, callback_dict, start_time, x0, PSOpopulation=None):
    if opt_params["global_opt_method"] == 'PSO':
        optimum = PSO.find_opt(opt_params, orig_params, callback_dict, target_data,population=PSOpopulation)

    elif opt_params["global_opt_method"] == 'dual_annealing':
        optimum = opti.dual_annealing(
            func=sim_and_calc_loss,
            x0=x0,
            maxiter=opt_params['max_steps'],
            #T=opt_params['temperature'],
            #niter_success=opt_params['max_steps_without_improvement'],
            bounds=callback_dict['bounds'],
            callback=default_callback(opt_params, orig_params, callback_dict),
            local_search_options={
                "method": opt_params['local_opt_method'],
                "args": (  # passed to func each call
                    opt_params,
                    orig_params,
                    target_data,
                    callback_dict
                ),
                "bounds": callback_dict['bounds'],
                #"callback":default_local_callback(opt_params, callback_dict),
                "options": {
                    "maxiter": opt_params['max_local_steps']
                    # futher kwargs for L-BFGS-B: maxcor, eps, see docs 
                    # would putting bounds both here and above fix step problem?
                }
            },
        )   
    elif opt_params["global_opt_method"] == 'basinhopping':
        optimum = opti.basinhopping(
            func=sim_and_calc_loss,
            x0=x0,
            niter=opt_params['max_steps'],
            T=opt_params['temperature'],
            niter_success=opt_params['max_steps_without_improvement'],
            minimizer_kwargs={
                "method": opt_params['local_opt_method'],
                "args": (  # passed to func each call
                    opt_params,
                    orig_params,
                    target_data,
                    callback_dict
                ),
                "bounds": callback_dict['bounds'],
                #"callback":default_local_callback(opt_params, callback_dict),
                "options": {
                    "maxiter": opt_params['max_local_steps']
                    # futher kwargs for L-BFGS-B: maxcor, eps, see docs 
                    # would putting bounds both here and above fix step problem?
                }
            },
            take_step=default_step(opt_params, callback_dict),
            callback=default_callback(opt_params, orig_params, callback_dict),
            #stepsize=opt_params['stepsize'], #use for dynamic stepsize, not sure if works with custom step
            #interval= opt_params['interval'], #only include interval if dynamic step size used
        )
    else:
        console_log(opt,params,'\nERROR: Unknown global optimize parameter ' +  opt_params["global_opt_method"] + '\n',verbose_lvl=-1)
        assert(False)

    stop_time = time()

    txt = "\nOptimisation finished after %s hours, using %s steps and %s actual trials." %((stop_time-start_time)/3600, callback_dict['num_steps'], callback_dict['actual_num_runs']) 
    txt += f"\nMinimum distance found is {optimum.fun} with parameters:{optimum.x}"
    txt += "Optimization exit message:" + str(optimum.message) + '\n'
    console_log(opt_params, txt, verbose_lvl=1)

    with open(os.path.join(opt_params['opt_output_dir'],'optimizer_results.json'), 'w') as fp:
        optimized_params = {opt_params['target_params'][i]['name']:optimum.x[i] for i in range(len(optimum.x))}
        data = {'execution time (hours)':(stop_time-start_time)/360, 'distance':float(optimum.fun), \
        'message':optimum.message,'optimized_params':optimized_params,'steps':int(callback_dict['num_steps']),'trials':int(callback_dict['actual_num_runs'])}

        json.dump(data, fp)

    if opt_params['plot']:
        console_log(opt_params,"\nGenerating plots to summarize optimization...",verbose_lvl=1)
        plot.opt(opt_params, orig_params, callback_dict, target_data, optimum.x)

    console_log(opt_params,'\n...Done.\n',verbose_lvl=1)


##################################### HELPER FUNCTIONS ###########################################

def init_check(opt_params, orig_params):
    for param in opt_params['target_params']:
        if param['name'] not in orig_params.keys():
            console_log(opt_params,"\nERROR attempting to optimize a parameter that does not exist in the model:", param['name'],'\n',verbose_lvl=0)
            return False
    return True


def console_log(opt_params,txt,verbose_lvl=1):
    if opt_params['verbose'] >= verbose_lvl:
        print(txt)
        with open(os.path.join(opt_params['opt_output_dir'],'log.txt'),'a') as file:
            file.write('\n'+ txt)

################################### RUN MODEL AND COMPARE WITH TARGET DATA ########################################

def sim_and_calc_loss(x, opt_params, orig_params, target_data, callback_dict,dset_indx=None):

    callback_dict['actual_num_runs'] += 1
    sim_data = copasi_with_targets(x, opt_params, orig_params, target_data, callback_dict,dset_indx=dset_indx)
    if sim_data is None:
        return np.inf
    loss = calc_loss(opt_params, orig_params, sim_data, target_data, callback_dict)
    return loss


def copasi_with_targets(x, opt_params, orig_params, target_data, callback_dict, dset_indx=None):

    for i in range(len(x)):
        if opt_params['round'] == -1:
            val = x[i]
        else:
            val = round(x[i],opt_params['round'])       
        txt = 'setting ' + opt_params['target_params'][i]['name'] + ' <- ' + str(val)
        console_log(opt_params,txt,verbose_lvl=2)
        orig_params[opt_params['target_params'][i]['name']] = val

    if dset_indx is None:
        return run.sim_multi(opt_params, orig_params, target_data, opt=True)
    else:
        these_params = deepcopy(opt_params)
        these_params['dataset'] = [these_params['dataset'][dset_indx]]
        callback_dict['dataset_freq'] += [dset_indx]
        return run.sim_multi(these_params, orig_params, target_data, opt=True)



def calc_loss(opt_params, orig_params, sim_data, target_data, callback_dict):

    losses = []
    for i in range(len(sim_data)):
        if sim_data[i] is None or len(sim_data[i]['Time']) < 2: # copasi error
            # TODO: better detection of copasi err's, ex if partial complete
            # TODO: better handling of copasi err's, exit that step of something
            console_log(opt_params, "\nWARNING: optimize.py inferred copasi error, setting loss to inf\n", verbose_lvl=0)
            loss = np.inf #callback_dict['min_loss']*10000000
        
        else:
            loss = 0
            for j in range(len(opt_params["dataset"][0]['sim'])):
                emp_species, sim_species = opt_params["dataset"][i]['empirical'][j], opt_params["dataset"][i]['sim'][j]
                loss += calc_distance(opt_params, orig_params, target_data, sim_data[i], emp_species, sim_species, callback_dict)
            loss /= len(opt_params["dataset"][0]['sim'])    

        losses += [loss]

    if opt_params['sim_distance_norm'] == 'inf':
        total_loss = np.linalg.norm(losses, ord=np.inf)
    else:
        order = opt_params['sim_distance_norm']
        total_loss = np.linalg.norm(losses, ord=order)
        total_loss = total_loss/math.pow(len(losses),1/order)

    return total_loss

def calc_distance(opt_params, orig_params, target_data, sim_data, emp_species, sim_species, callback_dict):

    cutoff = int(round(orig_params['duration']/orig_params['sim_stepsize']))+1

    sim = np.array(sim_data[sim_species]['avg']['val'], dtype=np.float64)
    emp = np.array(target_data[emp_species], dtype=np.float64)[:cutoff]
    if opt_params['debug']:
        assert(len(sim)==len(emp))

    # sim is initialized to match the first datapoint in target_data, so time 0 for sim is the first time entry in the target data
    #target_time = np.array(target_data['Time'])[1:cutoff]
    #target_time = np.array(target_data['Time']) - target_data['Time'][0]
    #sim_time = np.array(sim_data['Time'])

    # linear interpolation from previous modules optimizer
    if False: #TODO fix interpolation
        linear_approx_fun = (
            lambda x: (x - int(x)) * sim[int(x)+1]
            + (1 + int(x) - x) * sim[int(x)])
        linear_approx = np.array(list(map(linear_approx_fun, target_time)), dtype="float64")

    #for i in range(target):
    #   linear_approx = (target[i] - int(target[i])) * sim[int(i)]


    #sim_approx = []
    #for i in range(len(target_time)):
    #   indx = np.where(sim_time==int(target_time[i]))
    #   #print('index=',indx[0][0],'vs len(sim_time)=',len(sim_time))
    #   sim_approx += [sim[indx[0][0]]]

    if opt_params['distance_logscaled']:
        X = np.array(np.ma.log10(sim) - np.ma.log10(emp)) 
        #masking (np.ma) to avoid taking log of 0
    else:
        X = np.array([sim-emp])

    if opt_params['decay']:
        X*=callback_dict['decay_vec']
    if opt_params['empirical_uncertainty']:
        X*=callback_dict['emp_uncertainty']

    if opt_params['distance_norm'] == 'inf':
        order = np.inf
    else:
        order = opt_params['distance_norm']
    dist = np.linalg.norm(X, ord=order)

    if opt_params['distance_norm'] != 'inf':
        dist = dist/math.pow(len(X),1/opt_params['distance_norm'])
        # TODO: check that this normz makes intuitive sense

    return dist

###################################### BASINHOPPING FUNCTIONS ###########################################


def default_callback(opt_params, orig_params, callback_dict):
    def call(tuned_param_values, loss, accept):
        if loss < callback_dict['min_loss']:
            callback_dict['min_loss'] = loss
            callback_dict['best_tuned_param_values'] =  {opt_params['target_params'][i]['name']:tuned_param_values[i] for i in range(len(tuned_param_values))}
        callback_dict['loss_over_steps'] += [callback_dict['min_loss']]

        callback_dict['num_steps'] += 1
        if callback_dict['num_steps'] % (1/opt_params['console_freq']) == 0:
            txt = "Optimizer at step %s with an error of %s\ncurr best params are %s." %(callback_dict['num_steps'], callback_dict['min_loss'],callback_dict['best_tuned_param_values'])
            console_log(opt_params,txt,verbose_lvl=1)
        if callback_dict['num_steps'] % (1/opt_params['pickle_freq']) == 0:
            pickle_file = os.path.join(opt_params['opt_output_dir'],'step' + str(callback_dict['num_steps']) + '.pickle')
            with open(pickle_file, 'wb') as f:
                data = {'callback_dict':callback_dict,'opt_params':opt_params, 'orig_params':orig_params}
                pickle.dump(data,f)

        if callback_dict['min_loss'] <= opt_params['min_loss'] or callback_dict['num_steps'] >= opt_params["max_steps"]:
            # note that cannot rely on basinhopping's exit for #steps, since may have restarted run from a pickle
            return True #basin-hopping will stop

    return call


def default_local_callback(opt_params, callback_dict):
    def local_call(state):
        aadfion88nldf = 0 #nothing here yet, but could count local steps, ect
        # if use also uncomment from the bassinhopping params
    return local_call

def default_step(opt_params, callback_dict):
    def step(target_param_values):
        bounds = callback_dict['bounds']
        # one step for the basin hopping macro solver
        if opt_params['log_uniform_steps'] in [1,10]: 
            new_param_vals = np.power(10,(np.ma.log10(bounds[:, 1]) - np.ma.log10(bounds[:, 0])) * (np.random.random(len(target_param_values))) + np.ma.log10(bounds[:, 0]))
        elif opt_params['log_uniform_steps']==2: 
            new_param_vals = np.power(2,(np.ma.log2(bounds[:, 1]) - np.ma.log2(bounds[:, 0])) * (np.random.random(len(target_param_values))) + np.ma.log2(bounds[:, 0]))
        
        else:
            new_param_vals = (bounds[:, 1] - bounds[:, 0]) * np.random.random(len(target_param_values)) + bounds[:, 0]
        
        if callback_dict['num_steps'] % (1/opt_params['console_freq']) == 0:
            txt='\nBasinhopping testing new param values: ' + str(new_param_vals)
            console_log(opt_params,txt,verbose_lvl=2)
        return new_param_vals
    return step



##################################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Optimizer for CRN parameters.')
    parser.add_argument('PARAM_OR_PICKLE_FILE',
                        nargs=1,
                        type=str,
                        help='the parameter file or the pickle file (when used with --pickle)')
    parser.add_argument('--pickle',
                        action=argparse.BooleanOptionalAction,
                        help='run multiple simulations')
    args = parser.parse_args()

    param_or_pickle_file = args.PARAM_OR_PICKLE_FILE[0]
    if not os.path.isfile(param_or_pickle_file):
        print(f"\nERROR: Cannot find param or pickle file {param_or_pickle_file}. Check its path.\n")
        exit(1)

    else:
        if args.pickle:
            pickle_run(param_or_pickle_file)
        else:
            full_run(param_or_pickle_file)



