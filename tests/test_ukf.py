import os, sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.getcwd())
from power_systems import case3
from acquire_data import acquire_data
from routines import dse


EPS = 1e-6
ms = 1e-3

powsys = "ACIEEE3_SG4_Dynamics.xml";
ppc = case3.case3()

methods = [ "forwardEuler", "RK2" ];


# Settings
settings = { "time": 5 };
settings["etime"] = 0.5
settings["tstep"] = 20 * ms;
settings["method"] = "psat4th_ekf"
settings["nimethod"] = "forwardEuler"
settings["nSG"] = 1;
settings["metrics"] = True


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

# Dynamic State Estimation (DSE)
for i in range(len(methods)):
    settings["nimethod"] = methods[i]
    results = dse.rundse(settings, ppc, simdata)
    plt.plot(simdata["t[s]"][1:], results["rmse"], label = "RMSE_" + str(methods[i]))
    if i == 0:
        print("RESULT COMPARISON:")
    print("RMSE_" + str(methods[i]) + " = " + str(np.sum(results["rmse"])))   
    
plt.xlabel("time [s]")
plt.ylabel(r'$RMSE$')
plt.legend()
plt.show()