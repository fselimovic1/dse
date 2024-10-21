import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
import numpy as np
import math
import ybus
import random

EPS = 1e-6;
ms = 1e-3;
rf = 1e-2;

def busesconn(ppc, busi, busj):
    for i in range(ppc["branch"].shape[0]):
        if (ppc["branch"][i, 0] == busi and ppc["branch"][i, 1] == busj)\
            or  (ppc["branch"][i, 1] == busi and ppc["branch"][i, 0] == busj):
            return True
    return False

def xml4PF(settings, ppc):

    # Basic Power System params
    nb = ppc["bus"].shape[0]
    ng = ppc["gen"].shape[0]
    genIdx = (ppc["gen"][:, 0] - 1).astype(int)

    # Compute the admittance matrix
    Y = ybus.ybus(ppc);
    #Y[8, 8] = 5.3300 - 1j * 24.09;
    # MODEL SOLVER
    model = ET.Element("Model", attrib={"type": "NR", "domain": "real", "name": settings["filename"], "eps": str(settings["eps"])});

    # VARIABLES
    vars = ET.SubElement(model, "Vars", attrib= {"out": "true"})
    for i in range(nb):
        if ppc["bus"][i, 1] > 1:
            # magnitdue
            ET.SubElement(vars, "Var", attrib={"name": "V" + str(i + 1), "val": str(ppc["bus"][i, 7])});
        else: 
            # magnitdue
            ET.SubElement(vars, "Var", attrib={"name": "V" + str(i + 1), "val": str(1)});
        # phase
        ET.SubElement(vars, "Var", attrib={"name": "theta" + str(i + 1), "val": str(0)});

    # PARAMETERS
    params = ET.SubElement(model, "Params")#, attrib= {"dT": str(ms)})

    # bus injected powers
    for i in range(nb):
        # active power
        pi = 0;
        for j in range(ng):
            if ppc["gen"][j, 0] == i + 1:
                pi = pi + ppc["gen"][j, 1];
        pi = pi - ppc["bus"][i, 2];
        ET.SubElement(params, "Param", attrib={"name": "p" + str(i + 1) + "_0", "val": str(pi)}); 
        # reactive power
        ET.SubElement(params, "Param", attrib={"name": "q" + str(i + 1) + "_0", "val": str(-ppc["bus"][i, 3])});
    # voltage data
    for i in range(nb):
        # magnitude
        ET.SubElement(params, "Param", attrib={"name": "V" + str(i + 1) + "_0", "val": str(ppc["bus"][i, 7])});
        ET.SubElement(params, "Param", attrib={"name": "theta" + str(i + 1) + "_0", "val": str(str(ppc["bus"][i, 8]))});
    # admittance matrix elements
    for i in range(nb):
        for j in range(nb):
            # conductance
            ET.SubElement(params, "Param", attrib={"name": "G_" + str(i + 1) + "_" + str(j + 1), "val": str(np.real(Y[i, j]))}); 
            # sucseptanse
            ET.SubElement(params, "Param", attrib={"name": "B_" + str(i + 1) + "_" + str(j + 1), "val": str(np.imag(Y[i, j]))}); 
    
    # Nonlinear equations
    nleqs = ET.SubElement(model, "NLEqs")
    for i in range(nb):
        # equations for slack bus
        if ppc["bus"][i][1] == 3:
            # magnitude
            mEqStr = "V" + str(i + 1) + " - V" + str(i + 1) + "_0 = 0";
            ET.SubElement(nleqs, "Eq", attrib={"fx": mEqStr});
            # angle
            aEqStr = "theta" + str(i + 1) + " - theta" + str(i + 1) + "_0 = 0";
            ET.SubElement(nleqs, "Eq", attrib={"fx": aEqStr});
        # equations for pv bus
        elif ppc["bus"][i][1] == 2:
            # magnitude
            mEqStr = "V" + str(i + 1) + " - V" + str(i + 1) + "_0 = 0";
            ET.SubElement(nleqs, "Eq", attrib={"fx": mEqStr});
            # active power
            pEqStr = "V" + str(i + 1) + "^2*G_" + str(i + 1) + "_" + str(i + 1); 
            for j in range(nb):
                if j != i and busesconn(ppc, i + 1, j + 1):
                    pEqStr = pEqStr + "+ V" + str(i + 1) + "*V" + str(j + 1) + "*(G_" + str(i + 1) + "_" + str(j + 1) + \
                    "*cos(theta" + str(i + 1) + "- theta" + str(j + 1) + ") + B_" + str(i + 1) + "_" + str(j + 1) + "*sin(theta" + \
                    str(i + 1) + " - theta" + str(j + 1) + "))";
            pEqStr = pEqStr + " - p" + str(i + 1) + "_0"; 
            ET.SubElement(nleqs, "Eq", attrib={"fx": pEqStr});
        # equations for pq bus
        else:
            # active power
            pEqStr = "V" + str(i + 1) + "^2*G_" + str(i + 1) + "_" + str(i + 1); 
            for j in range(nb):
                if j != i and busesconn(ppc, i + 1, j + 1):
                    pEqStr = pEqStr + "+ V" + str(i + 1) + "*V" + str(j + 1) + "*(G_" + str(i + 1) + "_" + str(j + 1) + \
                    "*cos(theta" + str(i + 1) + "- theta" + str(j + 1) + ") + B_" + str(i + 1) + "_" + str(j + 1) + "*sin(theta" + \
                    str(i + 1) + " - theta" + str(j + 1) + "))";
            pEqStr = pEqStr + " - p" + str(i + 1) + "_0"; 
            ET.SubElement(nleqs, "Eq", attrib={"fx": pEqStr}); 
           # reactive power
            qEqStr = "-V" + str(i + 1) + "^2*B_" + str(i + 1) + "_" + str(i + 1); 
            for j in range(nb):
                if j != i and busesconn(ppc, i + 1, j + 1):
                    qEqStr = qEqStr + "+ V" + str(i + 1) + "*V" + str(j + 1) + "*(G_" + str(i + 1) + "_" + str(j + 1) + \
                    "*sin(theta" + str(i + 1) + "- theta" + str(j + 1) + ") - B_" + str(i + 1) + "_" + str(j + 1) + "*cos(theta" + \
                    str(i + 1) + " - theta" + str(j + 1) + "))";
            qEqStr = qEqStr + " - q" + str(i + 1) + "_0"; 
            ET.SubElement(nleqs, "Eq", attrib={"fx": qEqStr});

    # Convert the ElementTree to a string
    xml_str = ET.tostring(model, encoding='utf-8', method='xml')

    # Use minidom to pretty-print the XML string
    pretty_xml = minidom.parseString(xml_str).toprettyxml(indent="    ")

    with open("modelSolver/real/" + settings["filename"] + ".xml", "w", encoding='utf-8') as files:
        files.write('<?xml version="1.0" encoding="UTF-8"?>\n' + pretty_xml.split('<?xml version="1.0" ?>\n', 1)[1])