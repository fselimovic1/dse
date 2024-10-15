import numpy as np

def ybus(ppc):
    bus = ppc["bus"]
    branch = ppc["branch"]
    Y = np.zeros((bus.shape[0], bus.shape[0]), dtype=complex) 

    admittance = np.ones((branch.shape[0]))/(branch[:, 2] + 1j*branch[:, 3])
    transratio = np.zeros((branch.shape[0]), dtype=complex)
    idxline = np.nonzero(branch[:,8] == 0)
    idxtrans = np.nonzero(branch[:,8])
    transratio[idxline] =  np.exp(1j * branch[idxline, 9])
    transratio[idxtrans] = branch[idxtrans, 8] * np.exp(1j * branch[idxtrans, 9])

    transratio_conj = np.conj(transratio)

    toto = admittance + 1j/2 * branch[:, 4]
    fromfrom = toto / (transratio * transratio_conj)
    fromto = -admittance / transratio_conj
    tofrom = -admittance/transratio

    # Branch parameters
    for i in range(branch.shape[0]):
        Y[branch[i, 0].astype(int) - 1, branch[i, 0].astype(int) - 1] += fromfrom[i]
        Y[branch[i, 1].astype(int) - 1, branch[i, 1].astype(int) - 1] += toto[i]
        Y[branch[i, 0].astype(int) - 1, branch[i, 1].astype(int) - 1] += fromto[i]
        Y[branch[i, 1].astype(int) - 1, branch[i, 0].astype(int) - 1] += tofrom[i]

    # Shunt admittances and susceptances
    Y[bus[:,0].astype(int) - 1, bus[:,0].astype(int) - 1] += (bus[:,4] + 1j * bus[:,5])
    return Y


def l2z_ybus(ppc):

    Y = ybus(ppc)
        
    bus = ppc["bus"];
    for i in range(0, bus.shape[0]):
        Y[i, i] += (bus[i, 2] - 1j * bus[i, 3])/(bus[i, 7]**2);
    return Y;

def reduced_ybus(ppc):
    Yss = l2z_ybus(ppc);
    
    sg = ppc["sg"]
    nb = ppc["bus"].shape[0]
    ng = sg.shape[0]
     

    Ynn = np.zeros((ng, ng), dtype=complex)
    Yns = np.zeros((ng, nb), dtype=complex);

    for i in range(ng):
        no = sg[i, 0].astype(int) - 1
        Yg = 1/(1j * sg[i, 8]);

        Yss[no, no] += Yg
        # Additions to Ybus
        # Additional matrices
        Ynn[i, i] = Yg;
        Yns[i, no] = -Yg;

    Ysn = np.transpose(Yns);

    Yr = Ynn - Yns @ np.linalg.inv(Yss) @ Ysn;
    
    return Yr, Ynn, Yss, Yns;







