import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# Define unknown
t = sp.symbols('t');

def eval_fun(time, fun):

    n = len(time);
    efun = np.empty((n, 1));
    for i in range(n):
        point = {t: time[i]}
        efun[i] = fun.subs(point).evalf();
    return efun;

# Define constants
t0 = 0.05;
dT = 0.2;
edT = 1e-3;
rF = 0.25;

# Define x (state variable) as a function of t (time)
x =  -10 * t ** 2 + 2 * t + 2;
# Evaluate at t = t0;
point = { t: t0 };
x0 = x.subs(point).evalf();
t1 = t0 + dT;
point = { t: t1 };
x1 = x.subs(point).evalf();

# Plot function
time = np.arange(t0 - rF, t0 + rF, edT)
plt.plot(time, eval_fun(time, x), color = "slategrey", linewidth = 3);

# Plot trapezoidal
ts = [ t0, t0, t1, t1 ];
xs = [ 0, x0, x1, 0 ];
plt.fill(ts, xs, color = "wheat")

plt.ylim(0, None)

# Hiding numerical values on x and y axes
plt.xticks([])  # Hide x-axis numerical values
plt.yticks([])  # Hide y-axis numerical values
#plt.legend()

# Save the current figure
fig = plt.gcf()  # Get the current figure
fig.set_size_inches(8, 6)
fig.savefig(fname="illustrate_implicit.png", dpi=600)

# Display the plot
#plt.grid(True)
plt.show()