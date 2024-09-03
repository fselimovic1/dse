import numpy as np
import ybus
import math
from scipy.linalg import sqrtm


def redy_ukf(settings, ppc, simdata): 

    t = simdata["t[s]"]
    tstep = settings["estep"]
    tsteps = len(t) - 1 

    ws = 1
    wn = 2 * math.pi * ppc["fn"]

    nb = ppc["bus"].shape[0]
    ng = ppc["gen"].shape[0]


    vars = [];

    # SG data
    M = ppc["sg"][:, 17];
    D = ppc["sg"][:, 18];

    # define states' names
    for i in range(ng):
        vars.append("w" + str(i + 1));
    for i in range(ng):
        vars.append("delta" + str(i + 1));

    ns = 2 * ng;
    nm = 2* ng + 2* nb;
    modelorder = 2 *ng;

    Yr, Ynn, Yss, Yns = ybus.reduced_ybus(ppc)
    YV = -np.linalg.inv(Yss) @ np.transpose(Yns);
    G = np.real(Yr);
    B = np.imag(Yr);

    Pm = np.zeros(ng);
    Eg = np.zeros(ng);

    hx = np.zeros(2 * ng + 2 * nb)
    z_pg = np.zeros(ng);
    z_qg = np.zeros(ng);
    z_vm = np.zeros(nb);
    z_va = np.zeros(nb);

    # KALMAN FILTER PARAMETERS
    states = np.empty((ns, tsteps))
    states_plus = np.hstack((np.ones(ng), np.zeros(ng)))
    for k in range(ng):
        states_plus[k + ng] = simdata["delta" + str(k + 1)][1]
    states_minus = np.zeros(ns);

    # covariance matrices
    sigma_w = 1e-3
    sigma_v = 1e-2
    Q = sigma_w **2 *  np.eye(ns);
    R = sigma_v **2 * np.eye(nm);
    P = 1e-4 * np.eye(ns);

    # UKF parameters
    alfa = 0.25;
    beta = 2;
    kappa = 0;
    d = alfa **2 * (modelorder + kappa) - modelorder

    # Constants in this simulation
    for i in range(ng):
         Eg[i] = simdata["E" + str(i + 1)];
         Pm[i] = simdata["pm" + str(i + 1)];

    sigma_xp = np.empty((modelorder, 2 * modelorder + 1))
    sigma_fx = np.empty((modelorder, 2 * modelorder + 1))
    sigma_xm = np.empty((modelorder, 2 * modelorder + 1))
    sigma_hx = np.empty((nm, 2 * modelorder + 1)) 

    for i in range(tsteps):
        # Measurements and inputs (Tm, Ef, ...)
        for j in range(ng):
            z_pg[j] = simdata["pg" + str(j + 1)][i] * (1 + np.random.normal() * 1e-3);
            z_qg[j] = simdata["qg" + str(j + 1)][i] * (1 + np.random.normal() * 1e-3);
        
        for j in range(nb):
            z_vm[j] = simdata["V" + str(j + 1)][i] * (1 + np.random.normal() * 1e-3);
            z_va[j] = simdata["theta" + str(j + 1)][i] * (1 + np.random.normal() * 1e-3);
        
        z = np.hstack((z_pg, z_qg, z_vm, z_va))
        
        if i != 0:
            # CORRECTION STEP
            sqrtP = sqrtm(P)
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
                E = Eg * np.exp(1j * sigma_xm[ng:2 * ng, j]);
                V = YV @ E;
                Vm = np.abs(V);
                Va = np.angle(V)
                Sg = E * np.conjugate(Yr @ E);
                # expected measurement values
                sigma_hx[0:ng, j] = np.real(Sg);
                sigma_hx[ng:2*ng, j] = np.imag(Sg);
                sigma_hx[2*ng:2*ng + nb, j] = Vm;
                sigma_hx[2*ng + nb:2*ng + 2*nb, j] = Va;

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
            # angles from -pi to pi
            states_plus[ng:2 * ng] = np.angle(np.exp(1j * states_plus[ng:2 * ng]));
            P = P - K @ G @ np.transpose(K);
        
        # PREDICTION STEP
        # Predicted states
        sqrtP = sqrtm(P)
        sigma_xp[:, 0] = states_plus;
        states_minus = np.zeros(modelorder);
        P = np.zeros((modelorder, modelorder));
        # Compute sigma points
        for j in range(1, modelorder + 1):
            sigma_xp[:, j] = states_plus + math.sqrt(modelorder + d) * sqrtP[:, j - 1];
            sigma_xp[:, j + modelorder] = states_plus - math.sqrt(modelorder + d) * sqrtP[:, j - 1];
        # Compute nonlinear function values
        for j in range(2 * modelorder + 1):
            E = Eg * np.exp(1j * sigma_xp[ng: 2*ng, j]);
            Ig = Yr @ E;
            Pg = np.real(E * np.conjugate(Ig));
            sigma_fx[0:ng, j] = sigma_xp[0:ng, j] + (tstep/M) * (Pm - Pg - D * (sigma_xp[0:ng, j] - ws));
            sigma_fx[ng:2*ng, j] = np.angle(np.exp(1j * (sigma_xp[ng:2*ng, j] + tstep * wn * (sigma_xp[0:ng, j] - ws))));

            # Compute predicted states
            if j == 0:
                states_minus += d/(modelorder + d) * sigma_fx[:, j];
            else:
                states_minus += 1/(2 * modelorder + 2 * d) * sigma_fx[:, j];
        
        # Prediction error covariance matrix
        for j in range(2 * modelorder + 1):                                                      
            diff = sigma_fx[:, j] - states_minus;
            diff.shape = (modelorder, 1)
            if j == 0:
                P += (d/(modelorder + d) + (1 - alfa**2 + beta)) * (diff) @ np.transpose(diff);
            else:
                P += 1/(2 * modelorder + 2 * d) * (diff) @ np.transpose(diff);

        P = P + Q;
        states[:, i] = states_plus

    results = {"states": states}
    results["snames"] = vars

    return results;