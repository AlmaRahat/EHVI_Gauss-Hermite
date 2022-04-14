#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions to run experiments for EHVI calculations.

:Author:
   Alma Rahat   <a.a.m.rahat@swansea.ac.uk>
:Date:
   14 April 2022
:Copyright:
   Copyright (c)  Alma Rahat, Swansea University, 2022
"""

# imports
import numpy as np
from scipy.stats import wishart
import pandas as pd
from gauss_hermite import mGauss_hermite
from ehvi_wrapper import generate_txt_input_pf, run_ehvi
from evoalgos.performance import FonsecaHyperVolume as FH
from scipy.stats import multivariate_normal as MVN

def calc_hyp(P, samples, hpv, current_hpv, weights=None):
    """
    A function to calculate the hypervolume improvement for a range of samples. 

    Args.
    P: Pareto front, and r x c numpy array, where c is the number of objectvies,
    and r is the number of solutions. [r for row, c for columns]
    samples: from Gauss-hermite, or MC. Again, r x c numpy array.
    hpv: FonsecaHyperVolume object from evoalgos. 
    current_hpv: current hypervolume value for P; intentionally done this way to
    reduce the time of this particular routine.

    Returns. 
    EHVI for samples. 
    """
    improvement = [hpv.assess_non_dom_front(
                    np.concatenate([P,i[np.newaxis,:]])) - 
                    current_hpv for i in samples]
    if weights is not None: # in case of MC
        return np.dot(improvement, weights)
    else: # in case of weighted samples, e.g. GH
        return np.average(improvement), np.std(improvement)

def run_exp(mean, cov, pf, ref_vec, n_mc, n_gh, prune=0.2, 
            code_base = "EHVI/EHVI_3D/", conf_filename = "test.txt"):
    """
    A function to run a single experiment that produces approximation for EHVI
    using exact formula, MC and GH, for a given Gaussian predictive density,
    the current estimated Pareto front, and the reference vector. 

    Args. 
    mean: mean of the Gaussian predictive distribution. 
    cov: covariance of the Gaussian predictive distribution. 
    pf: Pareto front. 
    ref_vec: reference vector. 
    n_mc: number of MC samples to use. 
    n_gh: the number of GH points per dimension. This should be a list. 
    prune: the percentage of points to prune; vary between 0 and 1.
    code_base: exact formula code directory. 
    conf_filename: what to call the file that is expected by the exact code. 

    Returns.
    A list of results and the associated labels. 
    """
    hpv = FH(ref_vec)
    current_hv = hpv.assess_non_dom_front(pf)
    samples_mc = MVN.rvs(mean=mean, cov=cov, size=n_mc)
    mc_avg, mc_std = calc_hyp(pf, samples_mc, hpv, current_hv)
    mc_err = mc_std/np.sqrt(n_mc)
    gh_ests = []
    for ngh in n_gh:
        samples_gauss, weights = mGauss_hermite(ngh, mean, cov, prune=prune)
        gh_est = calc_hyp(pf, samples_gauss, hpv, current_hv, weights)
        gh_ests.append(gh_est)
    generate_txt_input_pf(pf, ref_vec, mean, np.sqrt(cov.diagonal()), filename=conf_filename)
    measured = run_ehvi(code_base, conf_filename)
    res = [mean, np.sqrt(cov.diagonal()), ref_vec, np.array(gh_ests), [mc_avg, mc_std, mc_err, measured]]
    labels = ["mean"]*len(mean) + ["std"]*len(cov) + ["ref"]*len(ref_vec) + [str(i) for i in n_gh] + \
                ["mc_avg", "mc_std", "mc_err", "measured"]
    return np.concatenate(res), labels

# random generation of mean

def get_mean(lb, ub):
    """Get a random mean vector within the box defined by lower (lb) and 
    upper (ub) bounds. 
    """
    return np.random.random(len(lb))*(ub-lb)+lb
def get_cov(lb, ub):
    """Generate a corrrelated random covariance matrix between [0,..,0] and ub-lb. 
    """
    return wishart.rvs(len(lb), np.eye(len(lb)))
def get_cov_indep(lb, ub):
    """Generate a diagonal covariance matrix between [0,..,0] and ub-lb. 
    """
    ndim = len(lb)
    cov = np.zeros((ndim, ndim)) 
    np.fill_diagonal(cov, np.random.uniform(low=np.zeros(ndim), high=ub-lb)) 
    return cov

def repeated_runs(pf, ref_vec, lb, ub, n_runs, data_file, prune=0.2, nmc = 10000, 
                  ngh = np.arange(3, 16, 1), code_base="EHVI/EHVI_3D/"):
    print(code_base)
    for i in range(n_runs):
        print("Starting experiment:", i+1)
        ndim = pf.shape[1]
        # mean behind
        mean = get_mean(lb, ub)
        cov = get_cov_indep(lb, ub)
        # current hpv
        hpv = FH(ref_vec)
        current_hv = hpv.assess_non_dom_front(pf)
        val = run_exp(mean, cov, pf, ref_vec, nmc, ngh, code_base=code_base)
        if os.path.exists(data_file):
            with open(data_file, 'ab+') as abc:
                np.savetxt(abc, val[0][:,np.newaxis].transpose(), delimiter=",")
        else:
            file = open(data_file, 'w+', newline="\n")
            with file:    
                write = csv.writer(file)
                write.writerow(val[1])
            with open(data_file, 'ab+') as abc:
                np.savetxt(abc, val[0][:,np.newaxis].transpose(), delimiter=",")
        print("==========")