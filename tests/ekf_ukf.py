import os, sys
sys.path.append(os.getcwd())
import add_paths
import matplotlib.pyplot as plt
import json

saveFigures = True

EPS = 1e-6
ms = 1e-3

filename = "IEEE9_LoadOnTest";
info = [ "EKF", "UKF" ];
nmethods = len(info)
colors = [ "orange", "slategrey" ]
variable2plot = "w1";


# Import simulation data
with open("data/" + filename + "_simdata" + '.json', 'r') as file:
    simdata = json.load(file)

# Import estimation data
estdata = {};
for i in range(len(info)):
    with open("data/" + filename + "_estdata" + info[i] + '.json', 'r') as file:
        estdata[info[i]] = json.load(file);

# Compare estimation results
plt.plot(simdata["t[s]"], simdata[variable2plot], label = "True", color="darkgreen", linewidth=2)
for i in range(nmethods): 
    plt.plot(estdata[info[i]]["t[s]"], estdata[info[i]][variable2plot], label = info[i], color=colors[i])   
plt.xlabel("time [s]")
plt.ylabel(rf'${variable2plot}$')
plt.legend()

if saveFigures:
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(7, 5)
    fig.savefig(fname="EKFvsUKF_w.png", dpi=800, bbox_inches='tight')
else: 
    plt.show()

# Plot errors of the estimators
plt.figure()
for i in range(nmethods): 
    plt.plot(estdata[info[i]]["t[s]"], estdata[info[i]]["errors"], label = info[i], color=colors[i]) 
plt.ylabel(r'$\xi$')
plt.xlabel("time [s]")
plt.legend()

if saveFigures:
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(7, 5)
    fig.savefig(fname="EKFvsUKF_errors.png", dpi=800, bbox_inches='tight')
else: 
    plt.show()

