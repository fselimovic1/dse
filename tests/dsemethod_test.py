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

from power_systems import case39
ppc = process_data(case39.case39())

# SETTINGS [SIMULATION]
settings = { "time": 5 };
settings["tdmethod"] = "constZ2nd"
settings["nimethod"] = "trapezoidal"
settings["tstep"] = 1 * ms;
settings["v2plot"] = "w"
settings["metrics"] = True

# EVENT settings
# 1 Load on: {"etype": "loadOn", "power": 5}
# 2 Bus fault: {"etype": "bbfault", "noBus": 4}
# 3 Line fault: {"etype": "lfault", "noLine": 23 }
settings["event"] = {"etype": "lfault", "noLine": 23 };
settings["etime_s"] = 1
settings["etime_e"] = 1.05

# SETTINGS [DSE]
dsemethods= [ "redy_ekf", "redy_ukf", "redy_enkf" ];
settings["estep"] = 10 * ms;
settings["ntrails"] = 20;
settings["metrics"] = True

# Simulate dynamic data
simdata = constZ2nd(settings, ppc);

# Dynamic State Estimation (DSE)
nmethods = len(dsemethods);
plt.figure();
for i in range(nmethods):
    rmse = np.zeros(int(len(simdata["t[s]"]) - 1));
    settings["dsemethod"] = dsemethods[i];
    for j in range(settings["ntrails"]): 
        results = dse.rundse(settings, ppc, simdata)
        rmse = rmse+ results["rmse"];
    rmse /= settings["ntrails"];
    # Plot result
    plt.plot(simdata["t[s]"][1:], rmse, label = "RMSE_" + dsemethods[i]);

plt.xlabel("Time [s]");
plt.ylabel("RMSE");
plt.legend()
plt.show();

