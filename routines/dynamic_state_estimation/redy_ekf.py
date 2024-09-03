import numpy as np
import ybus
import math

def redy_ekf(settings, ppc, simdata): 

    dT = int(settings["estep"]/settings["tstep"]);
    tsteps_sim = len(simdata["t[s]"]);
    tsteps_est = math.floor((tsteps_sim - 1)/dT);
    tstep = settings["estep"];

    ws = 1
    wn = 2 * math.pi * ppc["fn"]

    nb = ppc["bus"].shape[0]
    ng = ppc["gen"].shape[0]

    # SG data
    M = ppc["sg"][:, 17];
    D = ppc["sg"][:, 18];

    vars = [];

    # define states' names
    for i in range(ng):
        vars.append("w" + str(i + 1));
    for i in range(ng):
        vars.append("delta" + str(i + 1));

    ns = 2 * ng;
    nm = 2* ng + 2* nb;

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
    states = np.empty((ns, tsteps_est))
    states_plus = np.hstack((np.ones(ng), np.zeros(ng)))
    for k in range(ng):
        states_plus[k + ng] = simdata["delta" + str(k + 1)][1]
    states_minus = np.zeros(ns);

    # covariance matrices
    sigma_w = 1e-2
    sigma_v = 1e-3
    Q = sigma_w **2 *  np.eye(ns);
    R = sigma_v **2 * np.eye(nm);
    P = 1e-3 * np.eye(ns);
    
    # Constants in this simulation
    for i in range(ng):
         Eg[i] = simdata["E" + str(i + 1)][0];
         Pm[i] = simdata["Pm" + str(i + 1)][0];

    max_idx = 0;
    l = 0;
    for i in range(1, tsteps_sim, dT):
        # Measurements and inputs (Tm, Ef, ...)
        for j in range(ng):
            z_pg[j] = simdata["pg" + str(j + 1)][i] * (1 + np.random.normal() * 1e-3);
            z_qg[j] = simdata["qg" + str(j + 1)][i] * (1 + np.random.normal() * 1e-3);
        
        for j in range(nb):
            z_vm[j] = simdata["V" + str(j + 1)][i] * (1 + np.random.normal() * 1e-3);
            z_va[j] = simdata["theta" + str(j + 1)][i] * (1 + np.random.normal() * 1e-3);
        
        z = np.hstack((z_pg, z_qg, z_vm, z_va))
        
        if i != 1:
            # CORRECTION STEP
            dpg = np.zeros((ng, ns), dtype=float);
            dqg = np.zeros((ng, ns),  dtype=float);
            dvm = np.zeros((nb, ns));
            dva = np.hstack((np.zeros((nb, ng)), np.ones((nb, ng))));

            # expected measurement values
            E = Eg * np.exp(1j * states_minus[ng:2 * ng]);
            V = YV @ E;
            Vm = np.abs(V);
            Va = np.angle(V)

            Sg = E * np.conjugate(Yr @ E);
            hx[0:ng] = np.real(Sg);
            hx[ng:2*ng] = np.imag(Sg);
            hx[2*ng:2*ng + nb] = Vm;
            hx[2*ng + nb:2*ng + 2*nb] = Va;

            # measurements jacobian matrix
            for j in range(ng):
                for k in range(ns):
                    kk = np.mod(k, ng);
                    if j == k:
                        continue;
                    elif j + ng == k:
                        dpg[j, k] = np.sum(Eg[j] * Eg * (B[j, :] * np.cos(states_minus[j + ng] - states_minus[ng:2*ng]) - G[j, :] * np.sin(states_minus[j + ng] - states_minus[ng:2*ng]))) - Eg[j] ** 2 * B[j, j];
                        dqg[j, k] = np.sum(Eg[j] * Eg * (G[j, :] * np.cos(states_minus[j + ng] - states_minus[ng:2*ng]) + B[j, :] * np.sin(states_minus[j + ng] - states_minus[ng:2*ng]))) - Eg[j] ** 2 * G[j, j];
                    elif k >= ng:
                        dpg[j, k] = Eg[j] * Eg[kk] * (G[j, kk] * np.sin(states_minus[j] - states_minus[kk]) - B[j, kk] * np.cos(states_minus[j] - states_minus[kk])) 
                        dqg[j, k] = Eg[j] * Eg[kk] * (-G[j, kk] * np.cos(states_minus[j] - states_minus[kk]) - B[j, kk] * np.sin(states_minus[j] - states_minus[kk])) 
            H = np.vstack((dpg, dqg, dvm, dva));
            # Kalman gain
            K = P @ np.transpose(H) @ np.linalg.inv((H @ P @ np.transpose(H) + R))
            states_plus = states_minus + K @ (z - hx)
            states_plus[ng:2*ng] = np.angle(np.exp(1j * states_plus[ng:2*ng]));
            P = (np.eye(ns) - K @ H) @ P;
        
        # PREDICTION STEP
        # Predicted states
        E = Eg * np.exp(1j * states_plus[ng: 2*ng]);
        Ig = Yr @ E;
        Pg = np.real(E * np.conjugate(Ig));
        states_minus[0:ng] = states_plus[0:ng] + (tstep/M) * (Pm - Pg - D * (states_plus[0:ng] - ws));
        states_minus[ng:2*ng] = np.angle(np.exp(1j * (states_plus[ng: 2 *ng] + tstep * wn * (states_plus[0:ng] - ws))));

        # Prediction error covariance matrix
        dwF = np.zeros((ng, ns), dtype=float);
        ddeltaF = np.zeros((ng, ns), dtype=float);
        for j in range(ng):
            for k in range(ns):
                kk = np.mod(k, ng);
                if j == k:
                    dwF[j, k] = -(tstep/M[j]) * D[j] + 1
                    ddeltaF[j, k] = tstep * wn;
                elif j + ng == k:
                    dwF[j, k] = -tstep/M[j] * (np.sum(Eg[j] * Eg * (B[j, :] * np.cos(states_plus[j + ng] - states_plus[ng:2*ng]) - G[j, :] * np.sin(states_plus[j + ng] - states_plus[ng:2*ng]))) - Eg[j] ** 2 * B[j, j]);
                    ddeltaF[j, k] = 1;
                elif k >= ng:
                    dwF[j, k] = -tstep/M[j] * Eg[j] * Eg[kk] * (G[j, kk] * np.sin(states_plus[j] - states_plus[kk]) - B[j, kk] * np.cos(states_plus[j] - states_plus[kk])) 
        F = np.vstack((dwF, ddeltaF));
        P = F @ P @ np.transpose(F) + Q;
        states[:, l] = states_plus
        l = l + 1;
        max_idx = i;

    results = {"states": states}
    results["snames"] = vars
    # Save results
    results = {"vars" : vars};
    for i in range(len(vars)):
        results[vars[i]] = states[i, :].tolist();
    results["max_idx"] = int(max_idx);
    results["dT"] = dT;
    return results;