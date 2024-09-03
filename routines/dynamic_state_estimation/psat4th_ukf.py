import numpy as np
from scipy.linalg import sqrtm
import math

def psat4th_ukf(settings, ppc, simdata):

    wn = math.pi * 2 * ppc["fn"]
    ws = 1;
    noSG = settings["nSG"];
    modelorder = 4;
    nm = 2;
    t = simdata["t[s]"]
    tstep = settings["tstep"]
    tsteps = len(t) - 1  

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
    states = np.empty((modelorder, tsteps))
    states_plus = np.array([0.1634, 1, 0, 1.1394]);

    sigma_xp = np.empty((modelorder, 2 * modelorder + 1))
    sigma_fx = np.empty((modelorder, 2 * modelorder + 1))
    sigma_xm = np.empty((modelorder, 2 * modelorder + 1))
    sigma_hx = np.empty((nm, 2 * modelorder + 1)) 

    # UKF parameters
    alfa = 0.15;
    beta = 2;
    kappa = 0;
    d = alfa **2 * (modelorder + kappa) - modelorder

    
    # covariance matrices
    sigmaw = 0.5
    sigmav = 0.1
    Q = sigmaw **2 * np.eye(modelorder);
    R = sigmav**2 * np.eye(nm);
    P = 1 * np.eye(modelorder);

    for i in range(1, tsteps + 1):
        # PMU measurements
        Pe = simdata["pg" + str(noSG)][i] + np.random.normal(scale=0.02);
        Qe = simdata["qg" + str(noSG)][i] + np.random.normal(scale=0.02);
        z = np.array([Pe, Qe]);
        V = simdata["V" + str(noSG)][i] #+ np.random.normal(scale=0.005);
        theta = simdata["theta" + str(noSG)][i] #+ np.random.normal(scale=0.005);

        # Input paramters:
        Pm = simdata["pm" + str(noSG)][i] #+ np.random.normal(scale=0.005);
        Vf = simdata["vf" + str(noSG)][i] #+ np.random.normal(scale=0.005);
        
        #CORRECTION STEP
        if i != 1:
            sqrtP = np.abs(sqrtm(P))
            sigma_xm[:, 0] = states_minus;
            hx = np.zeros(nm)
            G = np.zeros((nm, nm));
            L = np.zeros((modelorder, nm));
            # Compute sigma points
            for j in range(1, modelorder + 1):
                sigma_xm[:, j] = states_minus + math.sqrt(modelorder + d) * sqrtP[:, j - 1];
                sigma_xm[:, j + modelorder] = states_minus - math.sqrt(modelorder + d) * sqrtP[:, j - 1];
            # Compute measurement values
            for j in range(2 * modelorder + 1):
                id = -V * (c1 * math.sin(sigma_xm[0, j] - theta) + c3 * math.cos(sigma_xm[0, j] - theta)) +  c1 * sigma_xm[2, j] + c3 * sigma_xm[3, j];
                iq = V * (c2 * math.sin(sigma_xm[0, j] - theta) - c1 * math.cos(sigma_xm[0, j] - theta)) - c2 * sigma_xm[2, j] + c1 * sigma_xm[3, j];
                # Active generated power
                sigma_hx[0, j] = V * (id * math.sin(sigma_xm[0, j] - theta) + iq * math.cos(sigma_xm[0, j] - theta))
                # Rective generated power
                sigma_hx[1, j] = V * (id * math.cos(sigma_xm[0, j] - theta) - iq * math.sin(sigma_xm[0, j] - theta))

                # Compute expected measurement values
                if j == 0:
                    hx += d/(modelorder + d) * sigma_hx[:, j];
                else:
                    hx += 1/(2 * modelorder + 2 * d) * sigma_hx[:, j];

            for j in range(2 * modelorder + 1):
                diffx = sigma_xm[:, j] - states_minus
                diffx.shape = (modelorder, 1)               
                diffhx = sigma_hx[:, j] - hx
                diffhx.shape = (nm, 1)     
                if j == 0:
                    G += (d/(modelorder + d) + (1 - alfa**2 + beta)) * ((diffhx) @ np.transpose(diffhx));
                    L += (d/(modelorder + d) + (1 - alfa**2 + beta)) * ((diffx) @ np.transpose(diffhx));
                else:
                    G += 1/(2 * modelorder + 2 * d) * ((diffhx) @ np.transpose(diffhx));
                    L += 1/(2 * modelorder + 2 * d) * ((diffx) @ np.transpose(diffhx));

            G = G + R;

            # Kalman gain
            K = L @ np.linalg.inv(G);
            # Computed state values
            states_plus = states_minus + (K @ (z - hx))
            P = P - K @ G @ np.transpose(K);

        # PREDICTION STEP
        sqrtP = np.abs(sqrtm(P))
        sigma_xp[:, 0] = states_plus;
        states_minus = np.zeros(modelorder);
        P = np.zeros((modelorder, modelorder));
        # Compute sigma points
        for j in range(1, modelorder + 1):
            sigma_xp[:, j] = states_plus + math.sqrt(modelorder + d) * sqrtP[:, j - 1];
            sigma_xp[:, j + modelorder] = states_plus - math.sqrt(modelorder + d) * sqrtP[:, j - 1];

        # Compute nonlinear function values
        for j in range(2 * modelorder + 1):
            # Delta function
            sigma_fx[0, j] = sigma_xp[0, j] + tstep * wn * (sigma_xp[1, j] - ws)
            # w function
            id = -V * (c1 * math.sin(sigma_xp[0, j] - theta) + c3 * math.cos(sigma_xp[0, j] - theta)) + c1 * sigma_xp[2, j] + c3 * sigma_xp[3, j]; 
            iq = V * (c2 * math.sin(sigma_xp[0, j] - theta) - c1 * math.cos(sigma_xp[0, j] - theta)) - c2 * sigma_xp[2, j] + c1 * sigma_xp[3, j];
            sigma_fx[1, j] = sigma_xp[1, j] + tstep/M * (Pm - V * (id * math.sin(sigma_xp[0, j] - theta) + iq * math.cos(sigma_xp[0, j] - theta)) - D*(sigma_xp[1, j] - ws))
            # - ra * (id**2 + iq**2)

            # e1d function
            sigma_fx[2, j] = sigma_xp[2, j] + tstep/T1q * (-sigma_xp[2, j] + (xq - x1q) * iq);
            # e1q function
            sigma_fx[3, j] = sigma_xp[3, j] + tstep/T1d * (-sigma_xp[3, j] + (xd - x1d) * id + Vf);

            # Compute predicted states and its covariance
            if j == 0:
                states_minus += d/(modelorder + d) * sigma_fx[:, j];
            else:
                states_minus += 1/(2 * modelorder + 2 * d) * sigma_fx[:, j];
        
        for j in range(2 * modelorder + 1):                                                      
            diff = sigma_fx[:, j] - states_minus;
            diff.shape = (modelorder, 1)
            if j == 0:
                P += (d/(modelorder + d) + (1 - alfa**2 + beta)) * (diff) @ np.transpose(diff);
            else:
                P += 1/(2 * modelorder + 2 * d) * (diff) @ np.transpose(diff);

        P = P + Q;

        # Store state estimates:
        states[:, i - 1] = states_plus;

    results = {"states": states}
    results["snames"] = vars

    return results;

                                                                    
            
        