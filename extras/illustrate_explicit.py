import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

t0 = 0.3;
dT = 0.15;
edT = 1e-3;
rF = 0.3;
rD = 0.2;
t = sp.symbols('t');

def eval_fun(time, fun):

    n = len(time);
    efun = np.empty((n, 1));
    for i in range(n):
        point = {t: time[i]}
        efun[i] = fun.subs(point).evalf();
    return efun;


# Define x (state variable) as a function of t (time)
x =  0.6 * (t) ** 2;
dx = sp.Derivative(x);
# Evaluate at t = t0;
point = { t: t0 };
x0 = x.subs(point).evalf();
t1 = t0 + dT;
point = { t: t1 };
x1 = x.subs(point).evalf();
point = { t: t0 };
k = dx.subs(point).evalf()
y = k * t + (x0 - k * t0);
point = { t: t0 };
y0 = y.subs(point).evalf()
point = { t: t1 };
y1 = y.subs(point).evalf()

# Plot function
time = np.arange(t0 - rF, t0 + rF, edT)
plt.plot(time, eval_fun(time, x), color= 'slategrey', linewidth = 3); # tangent_line(x), label='Tangent at x_0 = 1', color='red', linestyle=':')

# Plot derivative
time = np.arange(t0 - rD, t0 + rD, edT)
plt.plot(time, eval_fun(time, y), color = 'goldenrod', linewidth = 1); # tangent_line(x), label='Tangent at x_0 = 1', color='red', linestyle=':')


# Plot dotted borders
# x0
time = np.arange(0, t0, edT)
n = len(time);
plt.plot(time, x0 * np.ones((n, 1)), color = "black", linestyle=':');
# y1
time = np.arange(0, t1, edT)
n = len(time);
plt.plot(time, y1 * np.ones((n, 1)), color = "black", linestyle=':');

# t0
arrayX0 = np.arange(0, x0, edT);
n = len(arrayX0);
time = t0 * np.ones((n, 1));
plt.plot(time, arrayX0, color = "black", linestyle=':');
# y1
arrayX1 = np.arange(0, y1, edT);
n = len(arrayX1);
time = t1 * np.ones((n, 1));
plt.plot(time, arrayX1, color = "black", linestyle=':');

# Add labels and title
# Hiding numerical values on x and y axes
plt.xticks([])  # Hide x-axis numerical values
plt.yticks([])  # Hide y-axis numerical values
#plt.legend()

plt.xlim(0, None)
plt.ylim(0, None)

# Save the current figure
fig = plt.gcf()  # Get the current figure
fig.set_size_inches(8, 6)
fig.savefig(fname="illustrate_explicit.png", dpi=600)

# Display the plot
plt.grid(True)
plt.show()

