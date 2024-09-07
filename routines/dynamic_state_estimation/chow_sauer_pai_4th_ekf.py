import numpy as np
import sympy as sp
import math
EPS = 1e-4

def chow_sauer_pai_4th_ekf(settings, ppc, simdata):
    wn = math.pi * 2 * ppc["fn"]
    ws = 1;
    noSG = settings["nSG"];
    modelorder = 4;
    nm = 2;
    dT = int(settings["estep"]/settings["tstep"]);
    tsteps_sim = len(simdata["t[s]"]);
    tsteps_est = math.floor((tsteps_sim - 1)/dT);
    tstep = settings["estep"];

    # Synchronous Generator parameters
    M = ppc["sg"][noSG - 1, 17];
    D = ppc["sg"][noSG - 1, 18];
    ra = ppc["sg"][noSG - 1, 6];
    xd = ppc["sg"][noSG - 1, 7];
    x1d = ppc["sg"][noSG - 1, 8];
    T1d = ppc["sg"][noSG - 1, 10];
    xq = ppc["sg"][noSG - 1, 12];
    x1q = ppc["sg"][noSG - 1, 13];
    T1q = ppc["sg"][noSG - 1, 15];

    vars = ["delta" + str(noSG), "w" + str(noSG), "e1d" + str(noSG), "e1q" + str(noSG)]

    # Constants
    K = 1 / (ra**2 + x1d * x1q);
    c1 = ra * K;
    c2 = x1d * K;
    c3 = x1q * K;

    # KALMAN FILTER PARAMETERS - intitialization
    delta, w, e1d, e1q = sp.symbols('delta, w, e1d, e1q')
    states = np.empty((modelorder, tsteps_est))
    states_plus = np.array([0, 1, 0, 1]);
    # states_plus = np.array([0.163405815186294, 1, 0, 1.13941982995182]);
    F = np.zeros((modelorder, modelorder));

    # covariance matrices
    sigmaw = 1e-2;
    sigmav = 1e-1;
    Q = sigmaw **2 * np.eye(modelorder)
    R = sigmav**2 * np.eye(nm);
    P = 100 * np.eye(modelorder);

    l = 0
    for i in range(1, tsteps_sim, dT):
        # PMU measurements
        Pe = simdata["pg" + str(noSG)][i] * (1 + np.random.normal(scale=0.005));
        Qe = simdata["qg" + str(noSG)][i] * (1 + np.random.normal(scale=0.005));
        wm = simdata["w" + str(noSG)][i] * (1 + np.random.normal(scale=0.005));
        if nm == 3:
            z = np.vstack((Pe, Qe, wm));
        else:
            z = np.vstack((Pe, Qe));
        V = simdata["V" + str(noSG)][i] #* (1 + np.random.normal(scale=0.005));
        theta = simdata["theta" + str(noSG)][i] #* (1 + np.random.normal(scale=0.005));

        # Input paramters:
        Pm = simdata["Pm" + str(noSG)][i] #* (1 + np.random.normal(scale=0.005));
        Vf = simdata["Vf" + str(noSG)][i] #* (1 + np.random.normal(scale=0.005));


        # Define simbolic functions
        id = c1 * (e1d - V * sp.sin(delta - theta)) + c3 * (e1q - V * sp.cos(delta - theta));
        iq = -c2 * (e1d - V * sp.sin(delta - theta)) + c1 * (e1q - V * sp.cos(delta - theta))
        # id = -V * (c1 * sp.sin(delta - theta) + c3 * sp.cos(delta - theta)) + c1 * e1d + c3 * e1q;
        # iq = V * (c2 * sp.sin(delta - theta) - c1 * sp.cos(delta - theta)) - c2 * e1d + c1 * e1q;
        
        #CORRECTION STEP
        if i != 1:
            point = {delta: states_minus[0].item(), w:states_minus[1].item(), e1d:states_minus[2].item(), e1q:states_minus[3].item()}; 
            # Calculate expected measurement value
            hp = V * (iq * sp.cos(theta-delta) - id * sp.sin(theta-delta))
            hq = V * (id * sp.cos(theta-delta) + iq * sp.sin(theta-delta))
            if nm == 2:
                sph = sp.Matrix([hp, hq]);
            else:
                hw = w;
                sph = sp.Matrix([hp, hq, hw]);
            h = np.array(sph.subs(point))
            spH = sph.jacobian([delta, w, e1d, e1q]);
            H = np.array((spH.subs(point)).evalf()).astype(float)

            # Calculate Kalman gain
            K = P @ np.transpose(H) @ np.linalg.inv(H @ P @ np.transpose(H) + R)
            # Calculate estimated state values
            states_plus = states_minus +  (K @ (z - h)).reshape(-1, 1);
            # Calculated estimation covariance matrix
            P = (np.eye(modelorder) - K @ H) @ P @ np.transpose(np.eye(modelorder) - K @ H) + K @ R @ np.transpose(K)

        # PREDICTION STEP
        # Calculate predicted states
        point = {delta: states_plus[0].item(), w:states_plus[1].item(), e1d:states_plus[2].item(), e1q:states_plus[3].item()};       
        if settings["nimethod"] == "forwardEuler":
            fdelta = delta + tstep * wn * (w - ws)
            # fw = w + tstep/M * (Pm - V * (id * sp.sin(delta - theta) + iq * sp.cos(delta- theta)) - ra * (id**2 + iq **2) - D * (w - ws))
            fw = w + tstep/M * (Pm - e1d * id - e1q * iq - (x1q - x1d) * id * iq - D * (w - ws))
            fe1d = e1d + tstep / T1q * (-e1d + (xq - x1q) * iq);
            fe1q = e1q + tstep / T1d * (-e1q + (xd - x1d) * id + Vf);
        elif settings["nimethod"] == "RK2":
            # K1
            k1d = wn * (w - ws) * tstep;
            k1w = tstep/M * (Pm - V * (id * sp.sin(delta - theta) + iq * sp.cos(delta- theta)) - ra * (id**2 + iq **2) - D * (w - ws));
            k1e1d = tstep / T1q * (-e1d + (xq - x1q) * iq);
            k1e1q = tstep / T1d * (-e1q + (xd - x1d) * id + Vf);
            # K2
            id = -V * (c1 * sp.sin(delta + k1d - theta) + c3 * sp.cos(delta + k1d - theta)) + c1 * (e1d + k1e1d) + c3 * (e1q + k1e1q);
            iq = V * (c2 * sp.sin(delta + k1d - theta) - c1 * sp.cos(delta + k1d - theta)) - c2 * (e1d + k1e1d) + c1 * (e1q + k1e1q)
            k2d = wn * (w + k1w - ws) * tstep;
            k2w = tstep/M * (Pm - V * (id * sp.sin(delta + k1d - theta) + iq * sp.cos(delta + k1d - theta)) - ra * (id**2 + iq **2) - D * (w + k1w - ws));
            k2e1d = tstep / T1q * (-(e1d + k1e1d) + (xq - x1q) * iq);
            k2e1q = tstep / T1d * (-(e1q + k1e1q) + (xd - x1d) * id + Vf);

            # NI functions 
            fdelta = delta + 1/2 * (k1d + k2d)
            fw = w + 1/2 * (k1w + k2w)
            fe1d = e1d + 1/2 * (k1e1d + k2e1d)
            fe1q = e1q + 1/2 * (k1e1q + k2e1q)

        spf = sp.Matrix([fdelta, fw, fe1d, fe1q]);
        states_minus = np.array((spf.subs(point)).evalf()).astype(float) + 0
        spF = spf.jacobian([delta, w, e1d, e1q]);
        F = np.array((spF.subs(point)).evalf())
        P = (F @ P @ np.transpose(F) + Q).astype(float);
        states[:, l] = states_plus.reshape(modelorder);
        l = l + 1
        max_idx = i;
    
    # Save results
    results = {"vars" : vars};
    for i in range(len(vars)):
        results[vars[i]] = states[i, :].tolist();
    results["max_idx"] = int(max_idx);
    results["dT"] = dT;

    return results;