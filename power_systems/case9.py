# Copyright (c) 1996-2015 PSERC. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Power flow data for 9 bus, 3 generator case.
"""

from numpy import array

def case9():
    """Power flow data for 9 bus, 3 generator case.
    Please see L{caseformat} for details on the case file format.

    Based on data from Joe H. Chow's book, p. 70.

    @return: Power flow data for 9 bus, 3 generator case.
    """
    
    ppc = {"version": '2'}

    ppc["processed"] = True

    ##-----  Power Flow Data  -----##
    ## system MVA base
    ppc["baseMVA"] = 100.0

    ## nominal system frequency
    ppc["fn"] = 60

    ## bus data
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    ppc["bus"] = array([
        [1, 3, 0,    0, 0, 0, 1, 1.04, 0, 345, 1, 1.1, 0.9],
        [2, 2, 0,    0, 0, 0, 1, 1.025, 0.1620, 345, 1, 1.1, 0.9],
        [3, 2, 0,    0, 0, 0, 1, 1.025, 0.0814, 345, 1, 1.1, 0.9],
        [4, 1, 0,    0, 0, 0, 1, 1.0258, -0.0387, 345, 1, 1.1, 0.9],
        [5, 1, 0.9,  0.3, 0, 0, 1, 1.0127, -0.0644, 345, 1, 1.1, 0.9],
        [6, 1, 0,    0, 0, 0, 1, 1.0324,  0.0343, 345, 1, 1.1, 0.9],
        [7, 1, 1, 0.35, 0, 0, 1, 1.0159,  0.0127, 345, 1, 1.1, 0.9],
        [8, 1, 0,  0, 0, 0, 1,  1.0258,  0.0649, 345, 1, 1.1, 0.9],
        [9, 1, 1.25, 0.5, 0, 0, 1, 0.9956, -0.0696, 345, 1, 1.1, 0.9]
    ])

    ## generator data
    # bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2,
    # Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf
    ppc["gen"] = array([
        [1, 0.71,   0.2705, 300, -300, 1, 100, 1, 250, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [2, 1.63, 0.0665, 0, -300, 1, 100, 1, 300, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [3, 0.850, -0.1086, 0, -300, 1, 100, 1, 270, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ])

    ## branch data
    # fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    ppc["branch"] = array([
        [1, 4, 0,      0.0576, 0,     250, 250, 250, 0, 0, 1, -360, 360],
        [4, 5, 0.017,  0.092,  0.158, 250, 250, 250, 0, 0, 1, -360, 360],
        [5, 6, 0.039,  0.17,   0.358, 150, 150, 150, 0, 0, 1, -360, 360],
        [3, 6, 0,      0.0586, 0,     300, 300, 300, 0, 0, 1, -360, 360],
        [6, 7, 0.0119, 0.1008, 0.209, 150, 150, 150, 0, 0, 1, -360, 360],
        [7, 8, 0.0085, 0.072,  0.149, 250, 250, 250, 0, 0, 1, -360, 360],
        [8, 2, 0,      0.0625, 0,     250, 250, 250, 0, 0, 1, -360, 360],
        [8, 9, 0.032,  0.161,  0.306, 250, 250, 250, 0, 0, 1, -360, 360],
        [9, 4, 0.01,   0.085,  0.176, 250, 250, 250, 0, 0, 1, -360, 360]
    ])

    ## synchronous generator data
    # bus_i, baseMVA, vn, fn, modelType, xl, ra, xd, x1d, x2d, T1d, T2d, xq, x1q, x2q, T1q, T2q,
    # M, D, conn_status
    ppc["sg"] = array([
        [ 1, 100, 16.5, 60, 4, 0, 0.002, 0.1460, 0.0608, 0, 8.96, 0, 0.0969, 0.0969, 0, 0.310, 0, 47.28, 4.7, 1 ], 
        [ 2, 100, 18.0, 60, 4, 0, 0.002, 0.8958, 0.1198, 0, 6.00, 0, 0.8645, 0.1969, 0, 0.535, 0, 12.80, 2.5, 1 ],
        [ 3, 100, 13.8, 60, 4, 0, 0.002, 1.3125, 0.1813, 0, 5.89, 0, 1.2578, 0.2500, 0, 0.600, 0,  6.02, 1.8, 1 ],
    ])

    return ppc