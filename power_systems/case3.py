"""Data for IEEE 3 bus test case.
"""

from numpy import array

def case3():
    ppc = {"version": '2'}
    ppc["processed"] = True
    ##-----  Power Flow Data  -----##
    ## system MVA base
    ppc["baseMVA"] = 100.0

    ## system nominal frequency [Hz]
    ppc["fn"] = 60

    ## bus data
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    ppc["bus"] = array([
        [1,  3,  0,    0,   0, 0,  1, 1,    0,    6.9, 1, 1.06, 0.94],
        [2,  2,  0,   0,    0, 0,  1, 0.9900,   0.0269, 6.9, 1, 1.06, 0.94],
        [3,  1, 0.5556, 0.2778, 0, 0,  1, 0.9837, -0.0088, 6.9, 1, 1.06, 0.94]
    ])

    ## generator data
    # bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2,
    # Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf
    ppc["gen"] = array([
        [1, 0.3386, 0.3222, 1, -0.5, 1,  100, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [2,  0.2222, -0.0285, 1, -0.5, 0.99, 100, 1, 5,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ])

    ## branch data
    # fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    ppc["branch"] = array([
        [1,   3, 0.025, 0.075, 0, 9900, 0, 0, 0,     0, 1, -360, 360],
        [1,   3, 0.025, 0.075, 0, 9900, 0, 0, 0,     0, 1, -360, 360],
        [2,   3, 0.050, 0.150, 0, 9900, 0, 0, 0,     0, 1, -360, 360]
    ])

    ## synchronous generator data
    # bus_i, baseMVA, vn, fn, modelType, xl, ra, xd, x1d, x2d, T1d, T2d, xq, x1q, x2q, T1q, T2q,
    # M, D, conn_status
    ppc["sg"] = array([
        [ 1, 100, 1, 60, 2, 0, 0.002, 0.911, 0.408, 0, 4.2, 0, 0.58, 0.58, 0, 3, 0, 13.3617, 1, 1 ], # 13.3617
        [ 2, 100, 1, 60, 2, 0, 0.002, 0.911, 0.408, 0, 4.2, 0, 0.58, 0.58, 0, 3, 0, 20.3617, 1, 1 ],
    ])

    return ppc