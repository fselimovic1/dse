import os, sys
sys.path.append(os.getcwd())
import add_paths
import matplotlib.pyplot as plt
from process_data import process_data
from routines.time_domain import td
from ybus import ybus

from power_systems import case39
ppc = process_data(case39.case39())

EPS = 1e-6
ms = 1e-3

nimethods = [ "trapezoidal" ];

# Settings
settings = { "time": 50 };
settings["tdmethod"] = "constZ2nd";
settings["nimethod"] = "trapezoidal";
settings["tstep"] = 5 * ms;
settings["estep"] = 100 * ms;
settings["v2plot"] = "w";
settings["metrics"] = True;
# EVENT settings
# 1 Load on: {"etype": "loadOn", "power": 5}
# 2 Bus fault: {"etype": "bbfault", "noBus": 4}
# 3 Line fault: {"etype": "lfault", "noLine": 5 }
settings["event"] = {"etype": "loadOn", "power": 1 };
settings["etime_s"] = 1
settings["etime_e"] = 2

plt.figure(1)
# Run time domain simulations 
for i in range(len(nimethods)):
    settings["nimethod"] = nimethods[i];
    simdata = td.td(settings, ppc);
    # Plot results
    for j in range(1):
        vr = settings["v2plot"] + str(1 + 1); 
        plt.plot(simdata["t[s]"], simdata[vr], label = vr + "_" + nimethods[i])  

plt.legend()
plt.show()