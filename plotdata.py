import matplotlib.pyplot as plt
import json

filename = "ThreeBusDynamics_SG4";
variable2plot = "w";

plot_estdata = True
printforall = True

# Import simulation data
with open("data/" + filename + "_simdata" + '.json', 'r') as file:
    simdata = json.load(file)

# Import estimation data
if plot_estdata:
    with open("data/" + filename + "_estdata" + '.json', 'r') as file:
        estdata = json.load(file);

# if decentralizde DSE
for i in range(simdata["ng"]): 
    if (plot_estdata or not printforall) and i + 1 != estdata["nSG"]:
         continue
    var = variable2plot + str(i + 1);
    if plot_estdata: 
        plt.plot(simdata["t[s]"][1:estdata["max_idx"] + 1:estdata["dT"]], simdata[var][1:estdata["max_idx"] + 1:estdata["dT"]], label = "True")
        plt.plot(simdata["t[s]"][1:estdata["max_idx"] + 1:estdata["dT"]], estdata[var], label = "Est")
    else:
        plt.plot(simdata["t[s]"], simdata[var], label = var)

plt.xlabel("time [s]")
plt.ylabel(rf'${variable2plot}$')
plt.legend()
plt.show()



