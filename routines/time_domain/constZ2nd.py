import numpy as np
import ybus
import sympy as sp
import math
import random

rf = 1e-2;

def new_step_calculation(E, Y, YV, delta):
        Eg = E * np.exp(1j * delta)
        Ig = Y @ Eg
        Sg = Eg * np.conj(Ig);
        Pg = np.real(Sg);
        Qg = np.imag(Sg);
        v = YV @ Eg
        V = np.abs(v)
        theta = np.angle(v)

        return Eg, Ig, Sg, Pg, Qg, V, theta;

def constZ2nd(settings, ppc):
    # Syn Gen params
    Xd = ppc["sg"][:, 8]
    M = ppc["sg"][:, 17]
    D = ppc["sg"][:, 18]
    ws = 1
    wn = 2 * math.pi * ppc["fn"]

    nb = ppc["bus"].shape[0]
    ng = ppc["gen"].shape[0]
    genIdx = (ppc["gen"][:, 0] - 1).astype(int)

    simdata = {"t[s]": (np.arange(0, settings["time"], settings["estep"])).tolist()};
    tsteps = int(settings["time"] /settings["tstep"]) + 1
    tstep = settings["tstep"]

    # Pre-event matrices
    Yr_p, Ynn_p, Yss_p, Yns_p = ybus.reduced_ybus(ppc)
    Yr = Yr_p;
    YV_p = -np.linalg.inv(Yss_p) @ np.transpose(Yns_p);
    YV = YV_p
    G_p = np.real(Yr_p);
    G = G_p;
    B_p = np.imag(Yr_p);
    B = B_p;

    #### EVENT SETUP ####
    etype = settings["event"]["etype"];
    # During event matrices
    if etype == "loadOn":
        loadpower = settings["event"]["power"];
        totalload = sum(ppc["bus"][:, 2])
        load_toadd = loadpower * totalload/100
        loadbus = np.where(ppc["bus"][:, 2] != 0)[0]
        bus_toadd = loadbus[random.randint(0, len(loadbus) - 1)];
        ppc["bus"][bus_toadd, 2] += load_toadd;

        # compute matrices
        Yr_e, Ynn_e, Yss_e, Yns_e = ybus.reduced_ybus(ppc)
        YV_e = -np.linalg.inv(Yss_e) @ np.transpose(Yns_e);
        G_e = np.real(Yr_e);
        B_e = np.imag(Yr_e);
        ppc["bus"][bus_toadd, 2] -= load_toadd;

        # Post event matrices
        Yr_c = Yr_p;
        YV_c = YV_p;
        G_c = G_p;
        B_c = B_p;
    elif etype == "bbfault":
        # event matrices
        nbus = int(settings["event"]["noBus"] - 1);
        Yss_e = np.copy(Yss_p);
        Yss_e[nbus, :] = np.zeros(nb);
        Yss_e[:, nbus] = np.zeros(nb);
        Yr_e = Ynn_p - Yns_p @ np.linalg.pinv(Yss_e) @ np.transpose(Yns_p);
        YV_e = -np.linalg.pinv(Yss_e) @ np.transpose(Yns_p);
        G_e = np.real(YV_e);
        B_e = np.imag(YV_e);
        # post event matrices
        Yr_c = Yr_p;
        YV_c = YV_p;
        G_c = G_p;
        B_c = B_p;
    elif etype == "lfault":
        nline = int(settings["event"]["noLine"] - 1);
        frombus = int(ppc["branch"][nline, 0] - 1);
        tobus = int(ppc["branch"][nline, 1] - 1);
        rl, xl, bs = ppc["branch"][nline, [2, 3, 4]];
        zl = rl + 1j * xl;
        zd_a = (zl/2 * zl/2 + zl * rf)/rf;
        zd_s  = 2 * (zl/2 * zl/2 + zl * rf)/zl;
        
        backupBranch = ppc["branch"];
        ppc["branch"] = np.delete(ppc["branch"], nline, axis=0)

        # Post event matrices
        Yr_c, Ynn_c, Yss_c, Yns_c = ybus.reduced_ybus(ppc)
        YV_c = -np.linalg.inv(Yss_c) @ np.transpose(Yns_c);
        G_c = np.real(Yr_c);
        B_c = np.imag(Yr_c);


        # EVENT matrices
        Yss_e = np.copy(Yss_c);
        Yss_e[frombus, tobus] -=  1/(zd_a);
        Yss_e[tobus, frombus] -=  1/(zd_a);
        Yss_e[frombus, frombus] += (1/(zd_s) + 1j * bs/2);
        Yss_e[tobus, tobus] += (1/(zd_s) + 1j * bs/2);
        Yr_e = Ynn_c - Yns_c @ np.linalg.inv(Yss_e) @ np.transpose(Yns_c);
        YV_e = -np.linalg.inv(Yss_e) @ np.transpose(Yns_c);
        G_e = np.real(Yr_e);
        B_e = np.imag(Yr_e);

        ppc["branch"] = backupBranch
       

    V = np.zeros((nb, tsteps));
    theta = np.zeros((nb, tsteps));
    Pg = np.zeros((ng, tsteps));
    Qg = np.zeros((ng, tsteps));
    w = np.ones((ng, tsteps));
    delta = np.zeros((ng, tsteps));

    V[:, 0] = ppc["bus"][:, 7];
    theta[:, 0] = ppc["bus"][:, 8];
    Pg[:, 0] = ppc["gen"][:, 1];
    Qg[:, 0] = ppc["gen"][:, 2];
    Eg = V[genIdx, 0] + Qg[:, 0] * Xd / V[genIdx, 0] + 1j * Pg[:, 0] * Xd / V[genIdx, 0];
    delta[:, 0] = theta[genIdx, 0] + np.angle(Eg);
    E = abs(Eg)
    Eg = E * np.exp(1j * delta[:, 0])

    Ig = Yr_p @ Eg
    Sg = Eg * np.conj(Ig);
    Pg[:, 0] = np.real(Sg);
    Qg[:, 0] = np.imag(Sg);

    # Mechanical input power
    Pm = np.real(Sg);

    # Thorugh time steps
    for i in range(1, tsteps):

        # FORWARD EULER METHOD
        if settings["nimethod"] == "forwardEuler":
            delta[:, i] = delta[:, i - 1] + tstep * wn * (w[:, i - 1] - ws);
            w[:, i] = w[:, i - 1] + tstep/M * (Pm - Pg[:, i - 1] - D * (w[:, i - 1] - ws)); 
        
        # RK2 METHOD
        if settings["nimethod"] == "RK2":
            # K1
            delta_k1 = tstep * wn * (w[:, i - 1] - ws);
            w_k1 = tstep/M * (Pm - Pg[:, i - 1] - D * (w[:, i - 1] - ws)); 

            # K2
            _, _, _, Pgk, _, _, _ = new_step_calculation(E, Yr, YV, delta[:, i - 1] + delta_k1);
            delta_k2 = tstep * wn * ((w[:, i - 1] + w_k1) - ws);
            w_k2 = tstep/M * (Pm - Pgk - D * ((w[:, i - 1] + w_k1) - ws)); 

            delta[:, i] =  delta[:, i - 1] + (delta_k1 + delta_k2)/2;
            w[:, i] = w[:, i - 1] + (w_k1 + w_k2)/2;
        
        # RK4 METHOD
        if settings["nimethod"] == "RK4":
            # K1
            delta_k1 = tstep * wn * (w[:, i - 1] - ws);
            w_k1 = tstep/M * (Pm - Pg[:, i - 1] - D * (w[:, i - 1] - ws)); 

            # K2
            _, _, _, Pgk, _, _, _ = new_step_calculation(E, Yr, YV, delta[:, i - 1] + delta_k1/2);
            delta_k2 = tstep * wn * ((w[:, i - 1] + w_k1/2) - ws);
            w_k2 = tstep/M * (Pm - Pgk - D * ((w[:, i - 1] + w_k1/2) - ws)); 

            # K3
            _, _, _, Pgk, _, _, _ = new_step_calculation(E, Yr, YV, delta[:, i - 1] + delta_k2/2);
            delta_k3 = tstep * wn * ((w[:, i - 1] + w_k2/2) - ws);
            w_k3 = tstep/M * (Pm - Pgk - D * ((w[:, i - 1] + w_k2/2) - ws));

            # K4
            _, _, _, Pgk, _, _, _ = new_step_calculation(E, Yr, YV, delta[:, i - 1] + delta_k3);
            delta_k4 = tstep * wn * ((w[:, i - 1] + w_k3) - ws);
            w_k4 = tstep/M * (Pm - Pgk - D * ((w[:, i - 1] + w_k3) - ws));

            delta[:, i] = delta[:, i - 1] + (delta_k1 + 2 * delta_k2 + 2 * delta_k3 + delta_k4)/6;
            w[:, i] = w[:, i - 1] + (w_k1 + 2 * w_k2 + 2 * w_k3 + w_k4)/6;
        
        # Trapezoidal method
        if settings["nimethod"] == "trapezoidal":
            f_delta = wn * (w[:, i - 1] - ws);
            f_w = 1/M * (Pm - Pg[:, i - 1] - D * (w[:, i - 1] - ws));
            # Prediction error covariance matrix
            dwF = np.zeros((ng, 2 * ng));
            ddeltaF = np.zeros((ng, 2 * ng));
            for j in range(ng):
                for k in range(2 * ng):
                    if j == k:
                        dwF[j, k] = -1/M[j] * (np.sum(E[j] * E * (B[j, :] * np.cos(delta[j, i - 1] - delta[:, i - 1]) - G[j, :] * np.sin(delta[j, i - 1] - delta[:, i - 1]))) - E[j] ** 2 * B[j, j])
                    elif k < ng:
                        dwF[j, k] = -1/M[j] * E[j] * E[k] * (G[j, k] * np.sin(delta[j, i - 1] - delta[k, i - 1]) - B[j, k] * np.cos(delta[j, i - 1] - delta[k, i - 1]))  
                    elif j + ng == k:
                        dwF[j, k] = -(1/M[j]) * D[j] + 1
                        ddeltaF[j, k] = wn;

            F = np.vstack((ddeltaF, dwF));
            Ac = np.eye(2 * ng) - 1/2 * tstep * F
            dx = - np.linalg.inv(Ac) @ (tstep * -np.hstack((f_delta, f_w)))
            delta[:, i] = delta[:, i - 1] + dx[0:ng].reshape(-1)
            w[:, i] = w[:, i - 1] + dx[ng: 2 * ng].reshape(-1)

        if (tstep * i) < settings["etime_s"]:
            Yr = Yr_p;
            YV = YV_p;
            G = G_p;
            B = B_p;
        elif (tstep * i) > settings["etime_e"]:
            Yr = Yr_c;
            YV = YV_c;
            G = G_c;
            B = B_c;
        else:
            Yr = Yr_e;
            YV = YV_e;
            G = G_e;
            B = B_e;
        # Update values
        Eg, Ig, Sg, Pg[:, i], Qg[:, i], V[:, i], theta[:, i] = new_step_calculation(E, Yr, YV, delta[:, i]); 
        #delta[:, i] = np.angle(Eg)

    # Collect simulation data
    estimationMoments = np.arange(0, tsteps - 1, int(settings["estep"]/tstep))
    for i in range(nb):
        simdata["V" + str(i + 1)] = V[i, estimationMoments].tolist();
        simdata["theta" + str(i + 1)] = theta[i, estimationMoments].tolist();
    for i in range(ng):
        simdata["pg" + str(i + 1)] = Pg[i, estimationMoments].tolist();
        simdata["qg" + str(i + 1)] = Qg[i, estimationMoments].tolist();
        simdata["w" + str(i + 1)] = w[i, estimationMoments].tolist();
        simdata["delta" + str(i + 1)] = delta[i, estimationMoments].tolist();
        simdata["Pm"  + str(i + 1)] = Pm[i];
        simdata["E"  + str(i + 1)] = E[i];

    return simdata;
