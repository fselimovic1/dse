import os, sys
sys.path.append(os.getcwd())
import add_paths
import matplotlib.pyplot as plt
from acquire_data import acquire_data
from routines.time_domain import td
from routines.dynamic_state_estimation import dse
from routines.time_domain.constZ2nd import constZ2nd
from process_data import process_data
import numpy as np


EPS = 1e-6
ms = 1e-3

powsys = "ACIEEE3_SG2_Dynamics.xml";

from power_systems import case39
ppc = process_data(case39.case39())

# SETTINGS
settings = { "time": 5 };
settings["tdmethod"] = "constZ2nd"
settings["nimethod"] = "trapezoidal"
settings["tstep"] = 1 * ms;
settings["v2plot"] = "w"
settings["metrics"] = True

# EVENT settings
# 1 Load on: {"etype": "loadOn", "power": 5}
# 2 Bus fault: {"etype": "bbfault", "noBus": 4}
# 3 Line fault: {"etype": "lfault", "noLine": 5 }
settings["event"] = {"etype": "lfault", "noLine": 23 };
settings["etime_s"] = 1
settings["etime_e"] = 1.05

dsemethods= [ "redy_ekf", "redy_ukf" ];
settings["estep"] = 10 * ms;
settings["nSG"] = 1;

if settings["tdmethod"] == "modelSolver":
    # GENERATE MEASUREMENTS
    # Generate time - domain simulation configuration file (XML)
    # To be DONE
    # To be DONE
    # To be DONE


    # Generate measurements
    outputFile = "tdSolution.txt" # File with simulation results (transient stability)

    # Run transient stability simulation
    os.chdir("modelSolver");
    msCommand = "modelSolver DAE real real\\" + powsys + " ..\\" + outputFile + " 0 " + str(settings["tstep"]) + " " + str(settings["time"]);
    os.system(msCommand);
    #os.system("python plotSol.py ..\\" + outputFile)
    os.chdir("..")

    # Aquire PMU measurement and exact values of variables from transient stability simulation data file
    simdata = acquire_data(outputFile)

elif settings["tdmethod"] == "constZ2nd":
    simdata = constZ2nd(settings, ppc);

# Dynamic State Estimation (DSE)
for i in range(len(dsemethods)):
    settings["dsemethod"] = dsemethods[i]
    ppc = process_data(case39.case39())    
    results = dse.rundse(settings, ppc, simdata)

    plt.figure(1)
    plt.plot(simdata["t[s]"][1:], np.log10(results["rmse"]), label = "RMSE_" + dsemethods[i])
    plt.figure(2)
    plt.plot(simdata["t[s]"][1:], results["states"][0, :], label = "Est_" + dsemethods[i])


# True state values
plt.figure(2)
plt.plot(simdata["t[s]"][1:], simdata["w1"][1:], label = "True")
         
plt.figure(1)
plt.xlabel("time [s]")
plt.ylabel(r'$RMSE$')

plt.figure(2)
plt.xlabel("time [s]")
plt.ylabel(r'$w$')

# Adding legends to the plots
plt.figure(1)
plt.legend()
plt.figure(2)
plt.legend()

# Show all figures
plt.show()