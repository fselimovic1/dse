import os
import matplotlib.pyplot as plt
from routines.time_domain.xml4DAE import xml4DAE
from acquire_data import acquire_data
from process_data import process_data
import pandas as pd
import numpy as np
import json

from power_systems import case39
ppc = process_data(case39.case39())

saveSimData = True;
plotSimData = False;

EPS = 1e-6
ms = 1e-3

# Specify name of the prepared power system or type 'TBD' if you want to generate new power system XML file
# ThreeBusDynamics_SG4
powsys = "TBD";
nb = 9;
ng = 3;
# Insert filename in the case of new XML file
filename = "IEEE9_LoadOnTest";

# Simuluation SETTINGS
settings = { "time": 3 };
settings["tstep"] = 1 * ms;
settings["nimethod"] = "Trapezoidal";
settings["v2plot"] = "w";
settings["EPS"] = 1e-4;

# Generate XML settings
settings["filename"] = filename;
# simultaneous or partitioned
settings["intType"] = "partitioned";
settings["debug"] = False;
# EVENT settings
# 1 Load on: {"etype": "loadOn", "power": 5}
# 2 Bus fault: { "etype": "bbfault", "noBus": 4 }
# 3 Line fault: {"etype": "lfault", "noLine": 25 } -> to be updated
# 4 Line removal: {"etype": "lrem", "noLine": 25 }
settings["event"] = { "etype": "loadOn", "power": 2 }
settings["etime_s"] = 1.4
settings["etime_e"] = 4


if powsys == "TBD":
   # Generate time - domain simulation configuration file (XML)
    xml4DAE(settings, ppc)
else:
    settings["filename"] = powsys;

# Generate measurements
outputFile = "tdSolution.txt" # File with simulation results (transient stability)

# Run transient stability simulation
os.chdir("modelSolver");
msCommand = "modelSolver DAE real real\\" + settings["filename"] + ".xml" + " ..\\" + outputFile + " 0 " + str(settings["tstep"]) + " " + str(settings["time"]);
os.system(msCommand);
# os.system("python plotSol.py ..\\" + outputFile)
os.chdir("..")

# Save simulated data as CSV file
simdata = acquire_data(outputFile, False)

# Plot results
if plotSimData:
    plt.figure(1)
    for i in range(ng):
        vr = settings["v2plot"] + str(i  + 1); 
        plt.plot(simdata["t[s]"], simdata[vr], label = vr)
    #plt.plot(simdata["t[s]"][1:], np.array(simdata["theta1"][1:]) - np.array(simdata["theta2"][1:]), label = "w1")
    plt.legend()
    plt.show()

# Save simulation data for the purpose od DSE
if saveSimData:
    with open("data/" + settings["filename"] + "_simdata" + '.json', 'w') as file:
        # Complement with general data
        simdata["time"] = settings["time"];
        simdata["tstep"] = settings["tstep"];
        simdata["nb"] = ppc["bus"].shape[0]
        simdata["ng"] = ppc["gen"].shape[0]
        # Step 3: Write the dictionary to the file
        json.dump(simdata, file)
