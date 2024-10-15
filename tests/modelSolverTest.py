import os, sys
sys.path.append(os.getcwd())
import add_paths
import matplotlib.pyplot as plt
from process_data import process_data
from acquire_data import acquire_data
from routines.time_domain.xml4DAE_CMNNE import xml4DAE
import numpy as np


from power_systems import case39
ppc = process_data(case39.case39())

EPS = 1e-6
ms = 1e-3

# Test variable 
nimethods = [ "RK2", "Trapezoidal" ]; 
# Dynamic Simulation Settings
settings = { "time": 5 };
settings["tstep"] = 20 * ms;
settings["filename"] = "IEEE3";
# simultanous or partitioned
settings["intType"] = "simultaneous";
settings["EPS"] = 1e-4;

# EVENT settings
# 1 Load on: {"etype": "loadOn", "power": 5}
# 2 Bus fault: { "etype": "bbfault", "noBus": 4 }
# 3 Line fault: {"etype": "lfault", "noLine": 25 }
settings["event"] = { "etype": "loadOn", "power": 20 }
settings["etime_s"] = 1
settings["etime_e"] = 6


# Generate measurements
outputFile = "tdSolution.txt" # File with simulation results (transient stability)

# Run transient stability simulation
os.chdir("modelSolver");
plt.figure();
os.chdir("..")
for i in range(len(nimethods)):
    # Generate time - domain simulation configuration file (XML)
    print(os.getcwd())
    settings["nimethod"] = nimethods[i]; 
    xml4DAE(settings, ppc)
    os.chdir("modelSolver");
    msCommand = "modelSolver DAE real real\\" + settings["filename"] + ".xml" + " ..\\" + outputFile + " 0 " + str(settings["tstep"]) + " " + str(settings["time"]);
    os.system(msCommand);
    os.chdir("..")
    simdata = acquire_data(outputFile)
    # Plot result
    plt.plot(simdata["t[s]"][1:], simdata["w1"][1:], label = "w_" + nimethods[i]);


plt.legend()
plt.show()



