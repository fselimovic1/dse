import numpy as np
from routines.dynamic_state_estimation import chow_sauer_pai_4th_ekf, redy_ekf, redy_ukf, redy_enkf, kui4th_ekf, psat2nd_ekf, psat2nd_ukf, psat4th_ekf, psat4th_ukf

def performetrics(results, simdata):

    modelorder = results["states"].shape[0]
    tsteps = results["states"].shape[1]

    rmse = np.zeros((tsteps));
    for i in range(tsteps):
        for j in range(modelorder):
            rmse[i] += (results["states"][j, i] - simdata[results["snames"][j]][i + 1]) ** 2
        rmse[i] /= modelorder
        rmse[i] = np.sqrt(rmse[i])

    return rmse  
    

def rundse(settings, ppc, simdata):

    if settings["dsemethod"] == "redy_ekf":
        results = redy_ekf.redy_ekf(settings, ppc, simdata);
    elif settings["dsemethod"] == "redy_ukf":
        results = redy_ukf.redy_ukf(settings, ppc, simdata);
    elif settings["dsemethod"] == "redy_enkf":
        results = redy_enkf.redy_enkf(settings, ppc, simdata);
    elif settings["dsemethod"] == "psat2nd_ekf":
        results = psat2nd_ekf.psat2nd_ekf(settings, ppc, simdata);
    elif settings["dsemethod"] == "psat4th_ekf":
        results = psat4th_ekf.psat4th_ekf(settings, ppc, simdata);
    elif settings["dsemethod"] == "kui4th_ekf":
        results = kui4th_ekf.kui4th_ekf(settings, ppc, simdata);
    elif settings["dsemethod"] == "psat2nd_ukf":
        results = psat2nd_ukf.psat2nd_ukf(settings, ppc, simdata);
    elif settings["dsemethod"] == "psat4th_ukf":     
        results = psat4th_ukf.psat4th_ukf(settings, ppc, simdata);
    elif settings["dsemethod"] == "chow_sauer_pai_4th_ekf":
        results = chow_sauer_pai_4th_ekf.chow_sauer_pai_4th_ekf(settings, ppc, simdata);
    
    
    if settings["metrics"]:
        results["rmse"] = performetrics(results, simdata)
    
    return results;