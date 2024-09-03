import numpy as np
import math

def process_data(ppc):
    if ppc["processed"]:
        return ppc;

    # process bus data
    bus = ppc["bus"]
    bus[:, 2:5] /= ppc["baseMVA"];
    bus[:, 8] *= 2 * math.pi/360;
    ppc["bus"] = bus;

    # process generator data
    gen = ppc["gen"]
    gen[:, [1, 2, 3, 4, 8, 9]] /= ppc["baseMVA"];  
    ppc["gen"] = gen

    # process branch data
    branch = ppc["branch"]
    branch[:, 9] /= 2 * math.pi/360;
    ppc["branch"] = branch
    
    
    ppc["processed"] = True
    
    return ppc;




