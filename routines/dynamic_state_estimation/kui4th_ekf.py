import numpy as np
import math

def kui4th_ekf(settings, ppc, simdata):
    wn = math.pi * 2 * ppc["fn"]
    ws = 1;
    noSG = 2;
    modelorder = 4;
    nm = 1;
    t = simdata["t[s]"]
    tstep = settings["tstep"]
    tsteps = len(t) - 1 

    # Synchronous Generator parameters
    M = ppc["sg"][noSG - 1, 17];
    D = ppc["sg"][noSG - 1, 18];
    xd = ppc["sg"][noSG - 1, 7];
    x1d = ppc["sg"][noSG - 1, 8];
    T1d = ppc["sg"][noSG - 1, 10];
    xq = ppc["sg"][noSG - 1, 12];
    x1q = ppc["sg"][noSG - 1, 13];
    T1q = ppc["sg"][noSG - 1, 15];


    # KALMAN FILTER PARAMETERS - intitialization
    states = np.empty((modelorder, tsteps))
    states_plus = np.array([0.6, 1, 0, 0]);
    F = np.zeros((modelorder, modelorder));
    h = np.zeros(tsteps);
    hx = simdata["pg" + str(noSG)][1]
    # covariance matrices
    sigmaw = 0.08
    sigmav = 0.2
    Q = sigmaw **2 * np.eye(modelorder);
    R = sigmav**2 * np.eye(nm);
    P = 10 * np.eye(modelorder);

    # Time loop
    for i in range(1, tsteps + 1):
        # PMU measurements
        Pe = simdata["pg" + str(noSG)][i] + np.random.normal(scale=0.01);
        V = simdata["V" + str(noSG)][i] + np.random.normal(scale=0.01);

        # Input paramters:
        Pm = simdata["pg" + str(noSG)][i] + np.random.normal(scale=0.01);
        Vf = simdata["vf" + str(noSG)][i] + np.random.normal(scale=0.01);

        #CORRECTION STEP
        if i != 1:
            # Calculate expected measurement value
            hx = V/x1d * states_minus[2] * np.sin(states_minus[0]) + V**2 * 1/2 * (1/xq - 1/x1d) * np.sin(2 * states_minus[0]);
            # Calulate measurement Jacobian
            H = np.zeros((1, 4)) 
            H[0, 0] = V * states_minus[2] * np.cos(states_minus[0]) / x1d + V**2 * (1/xq - 1/x1d) * np.cos(2 * states_minus[0])
            H[0, 2]= V * np.sin(states_minus[0]) / x1d;

            # Calculate Kalman gain
            K = P @ np.transpose(H) @ np.linalg.inv(H @ P @ np.transpose(H) + R)
            # Calculate estimated state values
            states_plus = states_minus +  (K @ (Pe - hx)).reshape(-1, 1)
            # Calculated estimation covariance matrix
            P = (np.eye(modelorder) - K @ H) @ P @ np.transpose(np.eye(modelorder) - K @ H) + K @ R @ np.transpose(K)

        # PREDICTION STEP
        # Calculate predicted states
        xp_delta = states_plus[0] + tstep * wn * (states_plus[1] - ws);
        xp_w = states_plus[1] + tstep/M * (Pm - (V/x1d * states_plus[2] * np.sin(states_plus[0]) + V**2 * 1/2 * (1/xq - 1/x1d) * np.sin(2 * states_plus[0])) - D * (states_plus[1] - ws))
        xp_eq = states_plus[2] + tstep/T1d * (Vf - states_plus[2] - (xd - x1d) * ((states_plus[2] - V * np.cos(states_plus[0]))/x1d))
        xp_ed = states_plus[3] + tstep/T1q * (-states_plus[3] + (xq - x1q) * (V * np.sin(states_plus[0])/xq));
        states_minus = np.vstack((xp_delta, xp_w, xp_eq, xp_ed));
        # Calculate dF/dx
        # ddelta/dx
        F[0, 0] = 1
        F[0, 1] = tstep * wn
        # dw/dx 
        F[1, 0] = -tstep/M * (V * states_plus[2] * np.cos(states_plus[0]) / x1d + V**2 * (1/xq - 1/x1d) * np.cos(2 * states_plus[0]))
        F[1, 1] = -tstep * D /M + 1
        F[1, 2] = -tstep/M * V * np.sin(states_plus[0]) / x1d
        # deq/dx
        F[2, 0] = -tstep/T1d * (xd - x1d) * V * np.sin(states_plus[0]) / x1d
        F[2, 2] = 1 - tstep/T1d * (1 + (xd - x1d)/x1d)
        # ded/dx
        F[3, 0] = tstep/T1q * (xq - x1q) * V * np.cos(states_plus[0]) / xq
        F[3, 3] = 1 - tstep/T1q
        P = F @ P @ np.transpose(F) + Q;
        
        # Store computed state variables
        states[:, i - 1] = states_plus.reshape(modelorder)
    
    return states;