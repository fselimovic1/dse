import numpy as np
import sympy as sp
import math

def psat2nd_ekf(settings, ppc, simdata):
    wn = math.pi * 2 * ppc["fn"]
    ws = 1;
    noSG = settings["nSG"];
    modelorder = 2;
    nm = 2;
    t = simdata["t[s]"]
    tstep = settings["tstep"]
    tsteps = len(t) - 1 

    # Synchronous Generator parameters
    M = ppc["sg"][noSG - 1, 17];
    D = ppc["sg"][noSG - 1, 18];
    ra = ppc["sg"][noSG - 1, 6];
    x1d = ppc["sg"][noSG - 1, 8];

    vars = ["delta" + str(noSG), "w" + str(noSG), "e1d" + str(noSG), "e1q" + str(noSG)]


    # Constants
    K = 1 / (ra**2 + x1d ** 2);
    c1 = ra * K;
    c2 = x1d * K;
    c3 = x1d * K;

    # KALMAN FILTER PARAMETERS - intitialization
    delta, w = sp.symbols('delta, w')
    states = np.empty((modelorder, tsteps))
    states_plus = np.array([0, 1.1]);
    F = np.zeros((modelorder, modelorder));

    # covariance matrices
    sigmaw = 0.04
    sigmav = 0.1
    Q = sigmaw **2 * np.eye(modelorder)
    R = sigmav**2 * np.eye(nm);
    P = 100 * np.eye(modelorder);

    for i in range(1, tsteps + 1):
        # PMU measurements
        Pe = simdata["pg" + str(noSG)][i] * (1 + np.random.normal(scale=0.02));
        Qe = simdata["qg" + str(noSG)][i] * (1 + np.random.normal(scale=0.02));
        #wm =  simdata["w" + str(noSG)][i] * (1 + np.random.normal(scale=0.001));
        z = np.vstack((Pe, Qe));
        V = simdata["V" + str(noSG)][i] #* (1 + np.random.normal(scale=0.005));
        theta = simdata["theta" + str(noSG)][i] #* (1 + np.random.normal(scale=0.005));

        # Input paramters:
        Pm = simdata["pm" + str(noSG)][i] #* (1 + np.random.normal(scale=0.005));
        Vf = simdata["vf" + str(noSG)][i] #* (1 + np.random.normal(scale=0.005));


        # Define simbolic functions
        id = -V * (c1 * sp.sin(delta - theta) + c3 * sp.cos(delta - theta)) + c3 * Vf;
        iq = V * (c2 * sp.sin(delta - theta) - c1 * sp.cos(delta - theta)) + c1 * Vf;
        
        #CORRECTION STEP
        if i != 1:
            point = {delta: states_minus[0].item(), w:states_minus[1].item()};
            # Calculate expected measurement value
            hp = V * (id * sp.sin(delta - theta) + iq * sp.cos(delta- theta))
            hq = V * (id * sp.cos(delta - theta) - iq * sp.sin(delta- theta))
            # hw = w;
            sph = sp.Matrix([hp, hq]);
            h = np.array(sph.subs(point))
            spH = sph.jacobian([delta, w]);
            H = np.array((spH.subs(point)).evalf()).astype(float)

            # Calculate Kalman gain
            K = P @ np.transpose(H) @ np.linalg.inv(H @ P @ np.transpose(H) + R)
            # Calculate estimated state values
            states_plus = states_minus +  (K @ (z - h)).reshape(-1, 1)
            # Calculated estimation covariance matrix
            P = (np.eye(modelorder) - K @ H) @ P @ np.transpose(np.eye(modelorder) - K @ H) + K @ R @ np.transpose(K)

        # PREDICTION STEP
        # Calculate predicted states
        point = {delta: states_plus[0].item(), w:states_plus[1].item()};
        if settings["nimethod"] == "forwardEuler":
            fdelta = delta + tstep * wn * (w - ws)
            fw = w + tstep/M * (Pm - V * (id * sp.sin(delta - theta) + iq * sp.cos(delta- theta)) - ra * (id**2 + iq **2) - D * (w - ws))
        
        elif settings["nimethod"] == "RK2":
            # K1
            k1d = wn * (w - ws) * tstep;
            k1w = tstep/M * (Pm - V * (id * sp.sin(delta - theta) + iq * sp.cos(delta- theta)) - ra * (id**2 + iq **2) - D * (w - ws));
            # K2
            id = -V * (c1 * sp.sin(delta + k1d - theta) + c3 * sp.cos(delta + k1d - theta)) + c3 * Vf;
            iq = V * (c2 * sp.sin(delta + k1d - theta) - c1 * sp.cos(delta + k1d - theta)) + c1 * Vf;
            k2d = wn * (w + k1w - ws) * tstep
            k2w = tstep/M * (Pm - V * (id * sp.sin(delta + k1d - theta) + iq * sp.cos(delta + k1d - theta)) - ra * (id**2 + iq **2) - D * (w + k1w - ws));

            fdelta = delta + 1/2 * (k1d + k2d )
            fw = w +  1/2 * (k1w + k2w)

        elif settings["nimethod"] == "RK4":
            # K1
            k1d = wn * (w - ws) * tstep;
            k1w = tstep/M * (Pm - V * (id * sp.sin(delta - theta) + iq * sp.cos(delta- theta)) - ra * (id**2 + iq **2) - D * (w - ws));
            # K2
            id = -V * (c1 * sp.sin(delta + k1d/2 - theta) + c3 * sp.cos(delta + k1d/2 - theta)) + c3 * Vf;
            iq = V * (c2 * sp.sin(delta + k1d/2 - theta) - c1 * sp.cos(delta + k1d/2 - theta)) + c1 * Vf;
            k2d = wn * (w + k1w/2 - ws) * tstep
            k2w = tstep/M * (Pm - V * (id * sp.sin(delta + k1d/2 - theta) + iq * sp.cos(delta + k1d/2 - theta)) - ra * (id**2 + iq **2) - D * (w + k1w/2 - ws));
            # K3
            id = -V * (c1 * sp.sin(delta + k2d/2 - theta) + c3 * sp.cos(delta + k2d/2 - theta)) + c3 * Vf;
            iq = V * (c2 * sp.sin(delta + k2d/2 - theta) - c1 * sp.cos(delta + k2d/2 - theta)) + c1 * Vf;
            k3d = wn * (w + k2w/2 - ws) * tstep
            k3w = tstep/M * (Pm - V * (id * sp.sin(delta + k2d/2 - theta) + iq * sp.cos(delta + k2d/2 - theta)) - ra * (id**2 + iq **2) - D * (w + k2w/2 - ws));
            # K4
            id = -V * (c1 * sp.sin(delta + k3d - theta) + c3 * sp.cos(delta + k3d - theta)) + c3 * Vf;
            iq = V * (c2 * sp.sin(delta + k3d - theta) - c1 * sp.cos(delta + k3d - theta)) + c1 * Vf;
            k4d = wn * (w + k3w - ws) * tstep
            k4w = tstep/M * (Pm - V * (id * sp.sin(delta + k3d - theta) + iq * sp.cos(delta + k3d - theta)) - ra * (id**2 + iq **2) - D * (w + k3w - ws));        

            fdelta = delta + 1/6 * (k1d + 2 * k2d + 2 * k3d + k4d)
            fw = w +  1/6 * (k1w + 2 * k2w + 2 * k3w + k4w)

        # **** IMPLEMENTATION NOT FINISHED ****
        elif settings["nimethod"] == "trapezoidal":
            delta_f = wn * (w - ws)
            w_f = 1/M * (Pm - V * (id * sp.sin(delta - theta) + iq * sp.cos(delta- theta)) - ra * (id**2 + iq **2) - D * (w - ws))
            f = sp.Matrix([delta_f, w_f]);
            fx = np.array((f.subs(point)).evalf())
            dFx = f.jacobian([delta, w]);
            Fx = np.array((dFx.subs(point)).evalf()).astype(float)
            Ac = np.eye(modelorder) - tstep/2 * Fx
            dx = -np.linalg.inv(Ac) @ (-tstep * fx)
            fdelta = delta + dx[0]
            fw = w + dx[1]
            F = Fx


        spf = sp.Matrix([fdelta, fw]);
        states_minus = np.array((spf.subs(point)).evalf()).astype(float)
        if settings["nimethod"] =="trapezoidal":
            F = Fx
        else:
            spF = spf.jacobian([delta, w]);
            F = np.array((spF.subs(point)).evalf())

        P = (F @ P @ np.transpose(F) + Q).astype(float);

        states[:, i - 1] = states_plus.reshape(modelorder);

    results = {"states": states}
    results["snames"] = vars

    return results;