import numpy as np
import ybus
import sympy as sp
import math, random

rf = 1e8;

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

def constZ2ndS(settings, ppc):
    # Syn Gen params
    Xd = ppc["sg"][:, 8]
    M = ppc["sg"][:, 17]
    D = ppc["sg"][:, 18]
    ws = 1
    wn = 2 * math.pi * ppc["fn"]

    nb = ppc["bus"].shape[0]
    ng = ppc["gen"].shape[0]
    genIdx = (ppc["gen"][:, 0] - 1).astype(int)

    simdata = {"t[s]": np.arange(0, settings["time"] + settings["tstep"], settings["tstep"])}
    tsteps = simdata["t[s]"].size
    tstep = settings["tstep"]

    # Pre-event matrices
    Yr_p, Ynn_p, Yss_p, Yns_p = ybus.reduced_ybus(ppc)
    YV_p = -np.linalg.inv(Yss_p) @ np.transpose(Yns_p);
    G_p = np.real(Yr_p);
    B_p = np.imag(Yr_p);

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
        Yss_e[frombus, tobus] -= - 1/(zd_a);
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


    # INITIALIZATION
    # state variables
    x_delta = sp.Matrix(sp.symbols(f'x_delta1:{ng + 1}'))
    x_w = sp.Matrix(sp.symbols(f'x_w1:{ng + 1}'))
    # state equations
    fdelta = sp.zeros(ng, 1)
    fw = sp.zeros(ng, 1)

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

    for i in range(1, tsteps):
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
        
        point = {}
        for k in range(ng):
            # current values of variables
            point[x_delta[k]] = delta[k, i - 1]
            point[x_w[k]] = w[k, i - 1]  
            # functions
            fdelta[k] = wn * (x_w[k] - ws)
            fw[k] = 1/M[k] * (Pm[k] - Pg[k, i - 1] - D[k] * (x_w[k] - ws))  
        f =  sp.Matrix.vstack(fdelta, fw);
        x = sp.Matrix.vstack(x_delta, x_w);
        f_subs = f.subs(point)
        f_evaluated = f_subs.evalf()
        fvals = np.array(f_evaluated).astype(np.float64)
        # Explicit forward euler
        if settings["nimethod"] == "forwardEuler":
            delta[:, i] = delta[:, i - 1] + tstep * fvals[0:ng].reshape(-1)
            w[:, i] = w[:, i - 1] + tstep * fvals[ng: 2 * ng].reshape(-1)
        # Trapezoidal method
        elif settings["nimethod"] == "trapezoidal":
            derivative = f.jacobian(x)
            Jx = np.array((derivative.subs(point)).evalf()).astype(float)
            Ac = np.eye(2 * ng) - 1/2 * tstep * Jx
            dx = - np.linalg.inv(Ac) @ (tstep * -np.array(fvals))
            delta[:, i] = delta[:, i - 1] + dx[0:ng].reshape(-1)
            w[:, i] = w[:, i - 1] + dx[ng: 2 * ng].reshape(-1)

        Eg, Ig, Sg, Pg[:, i], Qg[:, i], V[:, i], theta[:, i] = new_step_calculation(E, Yr, YV, delta[:, i]); 
    
    # Collect simulation data
    for i in range(nb):
        simdata["V" + str(i + 1)] = V[i, :];
        simdata["theta" + str(i + 1)] = theta[i, :];
    
    for i in range(ng):
        simdata["pg" + str(i + 1)] = Pg[i, :];
        simdata["qg" + str(i + 1)] = Qg[i, :];
        simdata["w" + str(i + 1)] = w[i, :];
        simdata["delta" + str(i + 1)] = delta[i, :];
        simdata["pm"  + str(i + 1)] = Pm[i]
        simdata["E"  + str(i + 1)] = E[i]

    return simdata;
