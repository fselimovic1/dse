import os
import add_paths
import json
import matplotlib.pyplot as plt
from process_data import process_data
from acquire_data import acquire_data
from routines.dynamic_state_estimation import dse
from routines.time_domain.constZ2nd import constZ2nd
from power_systems import case3


filename = "ThreeBusDynamics_SG4";

EPS = 1e-6
ms = 1e-3

# Import simulation data
with open("data/" + filename + "_simdata" + '.json', 'r') as file:
    simdata = json.load(file)

# Read general settings
settings = { "time": simdata["time"] };
settings["tstep"] = simdata["tstep"];

# DSE settings
saveEstData = True
settings["dsemethod"] = "chow_sauer_pai_4th_ekf";
settings["nimethod"] = "forwardEuler";
settings["estep"] = 20 * ms;
settings["nSG"] = 1;

# Plot settings
toPlotResults = False;
settings["metrics"] = False;
settings["v2plot"] = "e1q";

# Dynamic State Estimation (DSE)
ppc = process_data(case3.case3())
results = dse.rundse(settings, ppc, simdata)
var = settings["v2plot"] + str(settings["nSG"]);

# Save estimation data
if saveEstData:  
    with open("data/" + filename + "_estdata" + '.json', 'w') as file:
        # Complement with general data
        results["nSG"] = settings["nSG"];
        # Step 3: Write the dictionary to the file
        json.dump(results, file);

# Plot estimation data
if toPlotResults:        
    # Estimation results
    var = settings["v2plot"] + str(settings["nSG"]);
    plt.plot(simdata["t[s]"][1:], simdata[var][1:], label = "True")
    plt.plot(simdata["t[s]"][1:], results["states"][3, :], label = "Est")

    # plt.plot(simdata["t[s]"][1:], simdata["w2"][1:], label = "w2")
    # plt.plot(simdata["t[s]"][1:], simdata[var][1:] + np.random.normal(size = simdata[var].size - 1, scale = 0.005), label = "Meas")
    plt.xlabel("time [s]")
    plt.ylabel(r'$w$')
    plt.legend()
    plt.show()

    #plt.plot(simdata["t[s]"][1:], simdata["pg1"][1:], label = "pg1")
    #plt.plot(simdata["t[s]"][1:], simdata["V1"][1:], label = "V1")
    #plt.plot(simdata["t[s]"][1:], simdata["vf1"][1:], label = "vf1")

    # plt.plot(simdata["t[s]"][1:], simdata["pg2"][1:], label = "pg")
    # plt.plot(simdata["t[s]"][1:], h, label = "pg_est")
    # plt.legend()
    # plt.show()
