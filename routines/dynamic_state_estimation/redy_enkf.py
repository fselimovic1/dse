import numpy as np
import ybus
import math
from scipy.linalg import sqrtm


def redy_enkf(settings, ppc, simdata): 

    t = simdata["t[s]"]
    tstep = settings["estep"]
    tsteps = len(t) - 1 

    ws = 1
    wn = 2 * math.pi * ppc["fn"]

    nb = ppc["bus"].shape[0]
    ng = ppc["gen"].shape[0]


    vars = [];
    
    # define states' names
    for i in range(ng):
        vars.append("w" + str(i + 1));
    for i in range(ng):
        vars.append("delta" + str(i + 1));

    ns = 2 * ng;
    nm = 2* ng + 2* nb;
    modelorder = 2 * ng;

    Yr, Ynn, Yss, Yns = ybus.reduced_ybus(ppc)
    YV = -np.linalg.inv(Yss) @ np.transpose(Yns);

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

    # covariance matrices - EnEKF parameters
    Nen = 100
    sigma_v = 1e-2
    sigma_w = 1e-3
    R = sigma_v **2 * np.eye(nm);

     # SG data
    M = ppc["sg"][:, 17];
    D = ppc["sg"][:, 18];

    # Constants in this simulation
    for i in range(ng):
         Eg[i] = simdata["E" + str(i + 1)];
         Pm[i] = simdata["pm" + str(i + 1)];
    
    x_m = np.empty((modelorder, Nen))
    x_p = np.empty((modelorder, Nen))
    zx = np.empty((nm, Nen))

    for i in range(tsteps):
        # Measurements and inputs (Tm, Ef, ...)
        for j in range(ng):
            z_pg[j] = simdata["pg" + str(j + 1)][i];
            z_qg[j] = simdata["qg" + str(j + 1)][i];
        
        for j in range(nb):
            z_vm[j] = simdata["V" + str(j + 1)][i];
            z_va[j] = simdata["theta" + str(j + 1)][i];
        
        z = np.hstack((z_pg, z_qg, z_vm, z_va))
        
        if i != 0:
            # CORRECTION STEP
            K = G @ np.linalg.inv(L + R);
            for j in range(Nen):
                x_p[:, j] = x_m[:, j] +  K @ (z * (1 + np.random.randn(nm) * sigma_v) - zx[:, j])

            states_plus = np.mean(x_p, axis = 1)
        else:
            for j in range(Nen):
                x_p[:, j] = states_plus + np.random.randn(modelorder) * sigma_w;

        # PREDICTION STEP
        # Generate ensambles
        for j in range(Nen):
            # Generate ensambles
            E = Eg * np.exp(1j * x_p[ng:2 * ng, j]);
            V = YV @ E;
            Vm = np.abs(V);
            Va = np.angle(V)
            Sg = E * np.conjugate(Yr @ E);
            Pg = np.real(Sg);
            Qg = np.imag(Sg);

            # f values
            x_m[0:ng, j] = x_p[0:ng, j] + (tstep/M) * (Pm - Pg - D * (x_p[0:ng, j] - ws)) + np.random.randn(ng) * sigma_w;
            x_m[ng:2*ng, j] = np.angle(np.exp(1j * (x_p[ng:2*ng, j] + tstep * wn * (x_p[0:ng, j] - ws)))) + np.random.randn(ng) * sigma_w;
        
        states_minus = np.mean(x_m, axis = 1);

        # z values
        for j in range(Nen):
            E = Eg * np.exp(1j * x_m[ng:2 * ng, j]);
            V = YV @ E;
            Vm = np.abs(V);
            Va = np.angle(V)
            Sg = E * np.conjugate(Yr @ E);
            Pg = np.real(Sg);
            Qg = np.imag(Sg);

            zx[0:ng, j] = Pg;
            zx[ng:2*ng, j] = Qg;
            zx[2*ng:2*ng + nb, j] = Vm;
            zx[2*ng + nb:2*ng + 2*nb, j] = Va;

        h = np.mean(zx, axis = 1);


        # Compute elements for Kalman gain
        G = np.zeros((modelorder, nm));
        L = np.zeros((nm, nm));
        for j in range(Nen):
            diff_x = x_m[:, j] - states_minus;
            diff_x.shape = (modelorder, 1) 
            diff_h = zx[:, j] - h;
            diff_h.shape = (nm, 1) 
            G += diff_x @ np.transpose(diff_h);
            L += diff_h @ np.transpose(diff_h);
        G = 1/Nen * G;
        L = 1/Nen * L;

        states[:, i] = states_plus

    results = {"states": states}
    results["snames"] = vars

    return results;