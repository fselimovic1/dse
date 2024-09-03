import os, sys
import matplotlib.pyplot as plt
sys.path.append(os.getcwd())
import add_paths
from process_data import process_data
from routines.time_domain import td

from power_systems import case9
ppc = process_data(case9.case9())

EPS = 1e-6
ms = 1e-3
saveFigure = True;
add_name = "_FE&RK4_";

# nimethods = [ "forwardEuler", "forwardEuler", "RK4" ]
#labels = [ "Forward Euler " + r"$\Delta t$" + "= 0.001 s", "Forward Euler " + r"$\Delta t$" + "= 0.01 s", "Runge-Kutta 4 " + r"$\Delta t$" + "= 0.01 s" ];
nimethods = [ "trapezoidal", "trapezoidal", "trapezoidal" ]
labels = [ "Trapezoidal " + r"$\Delta t$" + "= 0.001 s",  "Trapezoidal " + r"$\Delta t$" + "= 0.01 s", "Trapezoidal " + r"$\Delta t$" + "= 0.05 s" ];
# labels = [ "RK4 " + r"$\Delta t$" + "= 0.001 s", "Trapezoidal " + r"$\Delta t$" + "= 0.001 s" ];

tsteps = [ 1e-3, 1e-2, 5e-2 ];
var2plot = "w2";
xlabel = "time [s]";
ylabel = "Rotor speed "+ r"$w_2$";

# Settings
settings = { "time": 30 };
settings["tdmethod"] = "constZ2nd"

settings["tstep"] = 1 * ms;
settings["estep"] = 1 * ms;
settings["metrics"] = False
# EVENT settings
# 1 Load on: {"etype": "loadOn", "power": 5}
# 2 Bus fault: { "etype": "bbfault", "noBus": 4 }
# 3 Line fault: {"etype": "lfault", "noLine": 25 }
settings["event"] = { "etype": "bbfault", "noBus": 7 }
settings["etime_s"] = 1
settings["etime_e"] = 1.2

# Run time domain simulations 
for i in range(len(nimethods)):
    settings["nimethod"] = nimethods[i];
    settings["tstep"] = tsteps[i];
    settings["estep"] = tsteps[i];
    simdata = td.td(settings, ppc);
    plt.plot(simdata["t[s]"][1:], simdata[var2plot][1:], label = labels[i])

# Plot the results
plt.ylim([0.99, 1.035])
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend()

# Save the current figure
if saveFigure:
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(7, 5)
    fig.savefig(fname="nimethods" + add_name + "test.png", dpi=800, bbox_inches='tight')

plt.show()
