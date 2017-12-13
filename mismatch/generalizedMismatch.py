#Generalized mismatch simulations
#Eric Schulz, 2017

import numpy as np
import pandas as pd
import os
import sys
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, Matern, ConstantKernel as C
from bayesian_optimization import (BayesianOptimizer, GaussianProcessModel, 
                                   MinimalRegretSearch, EntropySearch, UpperConfidenceBound,
                                   ExpectedImprovement, ProbabilityOfImprovement,
                                   GPUpperConfidenceBound)
from bayesian_optimization.utils.optimization import global_optimization #https://github.com/jmetzen/bayesian_optimization

#Cluster id
clusterId = sys.argv[1]
#create meshgrid of all 100 teacher student combos
teacherStudentCombos = np.array(np.meshgrid(np.arange(0, 1.1, 0.1), np.arange(0, 1.1, 0.1))).T.reshape(-1,2)

l0,l1 = teacherStudentCombos[clusterId] #assign teacher and student lambdas



def generate_function(kernel, boundaries, n_train_points, random_state):
    """Draws target function from GP."""
    # Create GP prior
    gp_prior = GaussianProcessRegressor(kernel=kernel)
    
    # Sample function from GP by first sampling X values and then
    # sampling the y values for this X values for a specific function
    # from the GP
    X_train = np.random.RandomState(random_state).uniform(
        boundaries[:, 0], boundaries[:, 1], (n_train_points, boundaries.shape[0]))
    y_train = gp_prior.sample_y(X_train, 1, random_state=random_state)

    # Fit GP to the (X, y) samples of the sampled function
    gp_target = GaussianProcessRegressor(kernel=kernel)
    gp_target.fit(X_train, y_train)
    
    # Use the mean of the fitted GP as target function
    def target_fct(x):
        return gp_target.predict(x[np.newaxis, :])[0, 0]
    
    return target_fct

def create_optimizer(kernel, acquisition_function):
    """Create Bayesian optimizer for given GP kernel and acquisition function."""
    model = GaussianProcessModel(kernel=kernel, alpha=1e-3)
    if acquisition_function == "mrs":
        acquisition_function = \
            MinimalRegretSearch(model=model, n_gp_samples=1000, n_candidates=25,
                                n_trial_points=250, n_samples_y=51, point=False)
    if acquisition_function == "mrs_point":
        acquisition_function = \
            MinimalRegretSearch(model=model, n_gp_samples=1000, n_candidates=25,
                                n_trial_points=250, n_samples_y=51, point=True)
    elif acquisition_function == "es":
        acquisition_function = \
            EntropySearch(model=model, n_gp_samples=1000, n_candidates=25,
                          n_trial_points=250, n_samples_y=51)
    elif acquisition_function == "ucb":
        acquisition_function = \
            UpperConfidenceBound(model=model, kappa=5.0)
    elif acquisition_function == "gp_ucb":
        acquisition_function = \
            GPUpperConfidenceBound(model=model, const=0.5) #Beta=0.5
    elif acquisition_function == "ei":
        acquisition_function = \
            ExpectedImprovement(model=model)
    elif acquisition_function == "pi":
        acquisition_function = \
            ProbabilityOfImprovement(model=model)
    bayes_opt = BayesianOptimizer(model=model, optimizer="direct+lbfgs",
                                  acquisition_function=acquisition_function)
    return bayes_opt
    
    
    # We do not impose any model mismatch between the GP from which the 
# function was sampled and the GP used internally in Bayesian optimization.
# Moreover, we assume that GP in BO already "knows" the true length scales

#create lambdas here.

kernel = RBF(length_scale=l0, length_scale_bounds="fixed")
kernel_bo = RBF(length_scale=l1, length_scale_bounds="fixed")

n_train_points = 250  # used only for sampling the target function
boundaries = np.array([[0.0, 1.0], [0.0, 1.0]])  # Search space

def perform_run(setting, n_trials):
    """Execute a single run of a specific setting."""
    # Make results reproducible by using run index as random seed
    np.random.seed(q)
    
    y_regret = np.empty((n_trials / 2))  # Compute simple regret every 5 steps
    X_dist = np.empty((n_trials / 2)) # Compute distance of recommendation from optimum every 5 steps
    X_query = np.empty((n_trials, boundaries.shape[0]))  # Remember query points
    
    # Generate target function and compute (approximately) its optimum (X_opt) and the 
    # maximal value y_opt
    target_fct = generate_function(kernel, boundaries, n_train_points, random_state=q)
    X_opt = global_optimization(target_fct,  boundaries=boundaries, 
                                optimizer="direct", maxf=1000)
    y_opt = target_fct(X_opt) 

    # Create Bayesian optimizer and perform run
    bayes_opt = create_optimizer(kernel_bo,acquisition_function=setting)
    for i in range(n_trials):
        # Delect query point and generate noisy observation
        query = bayes_opt.select_query_point(boundaries)
        X_query[i] = query
        result = target_fct(query) + 1e-3 * np.random.random()
        # Update Bayesian optimizer
        bayes_opt.update(query, result)
        
        if i % 2 == 1:
            # Every 5 time steps: determine recommendation of 
            # Bayesian optimizer and compute its regret and distance to optimum
            X_sel = global_optimization(lambda x: bayes_opt.model.gp.predict(x[np.newaxis, :]), 
                                        boundaries=boundaries, 
                                        optimizer="direct", maxf=1000)
            y_sel = target_fct(X_sel)

            y_regret[i / 2] = y_opt - y_sel
            X_dist[i / 2] = np.sqrt(((X_opt - X_sel)**2).sum())
    
    # Store results of individual runs in a log file (not stricly necessary)
    #f = open("log/%s_%s" % (setting, run), 'w')
    #cPickle.dump((setting, run, y_regret, X_dist, X_query), f, protocol=-1)
    #f.close()
    
    return setting, y_regret, X_dist, X_query

#Running the simulation

#Settings
settings = ["gp_ucb"]   # acquisition functions

#simulate for 20 steps    
setting, y_regret, X_dist, X_query = perform_run(acq, 20)
#create file for output
outfile=pd.DataFrame({'regret': y_regret, 'distance': X_dist, 'setting': setting})
#save to scratch
outfile.to_csv('/home/ucabchu/Scratch/out'+str(q)+'.csv')
