import add_paths
import matplotlib.pyplot as plt
from process_data import process_data
from routines.time_domain import td
from ybus import ybus
import json

from power_systems import case3
ppc = process_data(case3.case3())

saveSimData = True;
plotSimData = False;
filename = "IEEE3_TEST"

EPS = 1e-6
ms = 1e-3

# Settings
settings = { "time": 5 };
settings["tdmethod"] = "constZ2nd"
settings["nimethod"] = "trapezoidal"
settings["tstep"] = 10 * ms;
settings["estep"] = 10 * ms;
settings["v2plot"] = "qg"
settings["metrics"] = False
# EVENT settings
# 1 Load on: {"etype": "loadOn", "power": 5}
# 2 Bus fault: { "etype": "bbfault", "noBus": 4 }
# 3 Line fault: {"etype": "lfault", "noLine": 25 }
settings["event"] = {"etype": "loadOn", "power": 5}
settings["etime_s"] = 1
settings["etime_e"] = 6

# Run time domain simulations 
simdata = td.td(settings, ppc);

# Plot results
plt.figure(1)
for i in range(ppc["gen"].shape[0]):
    vr = settings["v2plot"] + str(i + 1); 
    plt.plot(simdata["t[s]"], simdata[vr], label = vr)
plt.xlabel("time [s]")
plt.show()


"""
dL = r"$\delta_{30}$";
tL = r"$\theta_{30}$"
plt.plot(simdata["t[s]"], simdata["delta1"], label = dL)
plt.plot(simdata["t[s]"], simdata["theta30"], label = tL)
plt.xlabel("time [s]")
plt.ylabel("Angle [rad]")
plt.legend()
plt.savefig("DynAlg", dpi=400)
"""

# Save simulation data for the purpose od DSE
if saveSimData:
    with open("data/" + filename + "_simdata" + '.json', 'w') as file:
        # Complement with general data
        simdata["time"] = settings["time"];
        simdata["tstep"] = settings["tstep"];
        simdata["nb"] = ppc["bus"].shape[0]
        simdata["ng"] = ppc["gen"].shape[0]
        # Step 3: Write the dictionary to the file
        json.dump(simdata, file)
