import os, sys
sys.path.append(os.getcwd())
import add_paths
import matplotlib.pyplot as plt
import json

EPS = 1e-6
ms = 1e-3

filename = "IEEE9_LoadOnTest";
info = [ "EKF", "UKF" ];
nmethods = len(info)
variable2plot = "w6";


# Import simulation data
with open("data/" + filename + "_simdata" + '.json', 'r') as file:
    simdata = json.load(file)

# Import estimation data
estdata = {};
for i in range(len(info)):
    with open("data/" + filename + "_estdata" + info[i] + '.json', 'r') as file:
        estdata[info[i]] = json.load(file);

# Compare estimation results
plt.plot(simdata["t[s]"], simdata[variable2plot], label = "True")
for i in range(nmethods): 
    plt.plot(estdata[info[i]]["t[s]"], estdata[info[i]][variable2plot], label = variable2plot + info[i])   
plt.xlabel("time [s]")
plt.ylabel(rf'${variable2plot}$')
plt.legend()
plt.show()

# Plot errors of the estimators
plt.figure()
for i in range(nmethods): 
    plt.plot(estdata[info[i]]["t[s]"], estdata[info[i]]["errors"], label = variable2plot + info[i]) 

plt.legend()
plt.show()
