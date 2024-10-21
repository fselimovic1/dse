import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
import numpy as np
import math
import ybus
import random

EPS = 1e-6
ms = 1e-3

rf = 1e-2;

def busesconn(ppc, busi, busj):
    for i in range(ppc["branch"].shape[0]):
        if (ppc["branch"][i, 0] == busi and ppc["branch"][i, 1] == busj)\
            or  (ppc["branch"][i, 1] == busi and ppc["branch"][i, 0] == busj):
            return True
    return False

def power_network(settings, ppc, dict4xml):
    
    nb = ppc["bus"].shape[0];
    ng = ppc["gen"].shape[0];

    # Compute the admittance matrix
    Y = ybus.ybus(ppc);

    # VARIABLES
    for i in range(nb):
        # magnitude
        dict4xml["vars"] = np.append(dict4xml["vars"], {"name": "V" + str(i + 1), "val": str(1)});
        # angle
        dict4xml["vars"] = np.append(dict4xml["vars"], {"name": "theta" + str(i + 1), "val": str(0)});

    # PARAMETERS
    # ybus
    for i in range(nb):
        for j in range(nb):
            # G
            dict4xml["params"] = np.append(dict4xml["params"], {"name": "G" + str(i + 1) + "_" + str(j + 1), "val": str(np.real(Y[i, j])) });
            # B
            dict4xml["params"] = np.append(dict4xml["params"], {"name": "B" + str(i + 1) + "_" + str(j + 1), "val": str(np.imag(Y[i, j])) });
    # load powers
    for i in range(nb):
        # active
        dict4xml["params"] = np.append(dict4xml["params"], {"name": "pl" + str(i + 1), "val": str(ppc["bus"][i, 2])});
        # reactive
        dict4xml["params"] = np.append(dict4xml["params"], {"name": "ql" + str(i + 1), "val": str(ppc["bus"][i, 3])});

    # NLEqs
    for i in range(nb):
        if np.any(ppc["sg"][:, 0] == i + 1):
            continue
        # active power
        pEqStr = "V" + str(i + 1) + "^2*G" + str(i + 1) + "_" + str(i + 1); 
        for j in range(nb):
            if j != i and busesconn(ppc, i + 1, j + 1):
                pEqStr = pEqStr + "+ V" + str(i + 1) + "*V" + str(j + 1) + "*(G" + str(i + 1) + "_" + str(j + 1) + \
                "*cos(theta" + str(i + 1) + "- theta" + str(j + 1) + ") + B" + str(i + 1) + "_" + str(j + 1) + "*sin(theta" + \
                str(i + 1) + " - theta" + str(j + 1) + "))";
        pEqStr = pEqStr + " + pl" + str(i + 1); 
        dict4xml["nleqs"] = np.append(dict4xml["nleqs"], {"fx": pEqStr});
        # reactive power
        qEqStr = "-V" + str(i + 1) + "^2*B" + str(i + 1) + "_" + str(i + 1); 
        for j in range(nb):
            if j != i and busesconn(ppc, i + 1, j + 1):
                qEqStr = qEqStr + "+ V" + str(i + 1) + "*V" + str(j + 1) + "*(G" + str(i + 1) + "_" + str(j + 1) + \
                    "*sin(theta" + str(i + 1) + "- theta" + str(j + 1) + ") - B" + str(i + 1) + "_" + str(j + 1) + "*cos(theta" + \
                    str(i + 1) + " - theta" + str(j + 1) + "))";
        qEqStr = qEqStr + " + ql" + str(i + 1);
        dict4xml["nleqs"] = np.append(dict4xml["nleqs"], {"fx": qEqStr});
    
    ################################################## INITIALIZATION ##################################################
    # VARIABLES
    for i in range(nb):
        # magnitude
        dict4xml["init"]["vars"] = np.append(dict4xml["init"]["vars"], {"name": "V" + str(i + 1), "val": str(1)});
        # angle
        dict4xml["init"]["vars"] = np.append(dict4xml["init"]["vars"], {"name": "theta" + str(i + 1), "val": str(0)}); 
    # PARAMETERS  
    # bus injected powers
    for i in range(nb):
        # active power
        pi = 0;
        for j in range(ng):
            if ppc["gen"][j, 0] == i + 1:
                pi = pi + ppc["gen"][j, 1];
        pi = pi - ppc["bus"][i, 2];
        dict4xml["init"]["params"] = np.append(dict4xml["init"]["params"], {"name": "p" + str(i + 1) + "_0", "val": str(pi)});
        # reactive power
        dict4xml["init"]["params"] = np.append(dict4xml["init"]["params"], {"name": "q" + str(i + 1) + "_0", "val": str(-ppc["bus"][i, 3])});
    # voltage data
    for i in range(nb):
        # magnitude
        dict4xml["init"]["params"] = np.append(dict4xml["init"]["params"], {"name": "V" + str(i + 1) + "_0", "val": str(ppc["bus"][i, 7])});
        #angle
        dict4xml["init"]["params"] = np.append(dict4xml["init"]["params"], {"name": "theta" + str(i + 1) + "_0", "val": str(ppc["bus"][i, 8])});
    # additonal parameters/injected currents
    for i in range(nb):
        # real injected current
        dict4xml["init"]["params"] = np.append(dict4xml["init"]["params"], {"name": "I" + str(i + 1) + "r", "val": str(0)});
        # imag injected current
        dict4xml["init"]["params"] = np.append(dict4xml["init"]["params"], {"name": "I" + str(i + 1) + "im", "val": str(0)});

    # NLEqs
    for i in range(nb):
        # equations for slack bus
        if ppc["bus"][i][1] == 3:
            # magnitude
            mEqStr = "V" + str(i + 1) + " - V" + str(i + 1) + "_0 = 0";
            dict4xml["init"]["nleqs"] = np.append(dict4xml["init"]["nleqs"], {"fx": mEqStr});
            # angle
            aEqStr = "theta" + str(i + 1) + " - theta" + str(i + 1) + "_0 = 0";
            dict4xml["init"]["nleqs"] = np.append(dict4xml["init"]["nleqs"], {"fx": aEqStr});
        # equations for pv bus
        elif ppc["bus"][i][1] == 2:
            # magnitude
            mEqStr = "V" + str(i + 1) + " - V" + str(i + 1) + "_0 = 0";
            dict4xml["init"]["nleqs"] = np.append(dict4xml["init"]["nleqs"], {"fx": mEqStr});
            # active power
            pEqStr = "V" + str(i + 1) + "^2*G" + str(i + 1) + "_" + str(i + 1); 
            for j in range(nb):
                if j != i and busesconn(ppc, i + 1, j + 1):
                    pEqStr = pEqStr + "+ V" + str(i + 1) + "*V" + str(j + 1) + "*(G" + str(i + 1) + "_" + str(j + 1) + \
                    "*cos(theta" + str(i + 1) + "- theta" + str(j + 1) + ") + B" + str(i + 1) + "_" + str(j + 1) + "*sin(theta" + \
                    str(i + 1) + " - theta" + str(j + 1) + "))";
            pEqStr = pEqStr + " - p" + str(i + 1) + "_0"; 
            dict4xml["init"]["nleqs"] = np.append(dict4xml["init"]["nleqs"], {"fx": pEqStr});
        # equations for pq bus
        else:
            # active power
            pEqStr = "V" + str(i + 1) + "^2*G" + str(i + 1) + "_" + str(i + 1); 
            for j in range(nb):
                if j != i and busesconn(ppc, i + 1, j + 1):
                    pEqStr = pEqStr + "+ V" + str(i + 1) + "*V" + str(j + 1) + "*(G" + str(i + 1) + "_" + str(j + 1) + \
                    "*cos(theta" + str(i + 1) + "- theta" + str(j + 1) + ") + B" + str(i + 1) + "_" + str(j + 1) + "*sin(theta" + \
                    str(i + 1) + " - theta" + str(j + 1) + "))";
            pEqStr = pEqStr + " - p" + str(i + 1) + "_0"; 
            dict4xml["init"]["nleqs"] = np.append(dict4xml["init"]["nleqs"], {"fx": pEqStr});
           # reactive power
            qEqStr = "-V" + str(i + 1) + "^2*B" + str(i + 1) + "_" + str(i + 1); 
            for j in range(nb):
                if j != i and busesconn(ppc, i + 1, j + 1):
                    qEqStr = qEqStr + "+ V" + str(i + 1) + "*V" + str(j + 1) + "*(G" + str(i + 1) + "_" + str(j + 1) + \
                    "*sin(theta" + str(i + 1) + "- theta" + str(j + 1) + ") - B" + str(i + 1) + "_" + str(j + 1) + "*cos(theta" + \
                    str(i + 1) + " - theta" + str(j + 1) + "))";
            qEqStr = qEqStr + " - q" + str(i + 1) + "_0"; 
            dict4xml["init"]["nleqs"] = np.append(dict4xml["init"]["nleqs"], {"fx": qEqStr});
    # PostProcessing
    for i in range(nb):
        # magnitude
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": "base.V" + str(i + 1) + " = V" + str(i + 1)});
        # angle
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": "base.theta" + str(i + 1) + " = theta" + str(i + 1)});
    # compute injected currents
    for i in range(nb):
        iReStr = "I" + str(i + 1) + "r = ";
        iImStr = "I" + str(i + 1) + "im = ";
        for j in range(nb):
            # real
            iReStr = iReStr + "+ V" + str(j + 1) + " * (" + "G" + str(i + 1) + "_" + str(j + 1) + "*cos(theta" + str(j + 1)\
                + ") - B" + str(i + 1) + "_" + str(j + 1) + "*sin(theta" + str(j + 1) + "))";
            # imag
            iImStr = iImStr + "+ V" + str(j + 1) + " * (" + "G" + str(i + 1) + "_" + str(j + 1) + "*sin(theta" + str(j + 1)\
                + ") + B" + str(i + 1) + "_" + str(j + 1) + "*cos(theta" + str(j + 1) + "))";
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": iReStr});
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": iImStr});  
    ####################################################################################################################

def synchronous_generator(settings, ppc, dict4xml):
    
    nb = ppc["bus"].shape[0];
    nsg = ppc["sg"].shape[0];
    statevars = [ 'delta', 'w', 'e1q', 'e1d', ];
    statevals = ["0", "1"]
    # wn & ws
    dict4xml["params"] = np.append(dict4xml["params"], {"name": "wn", "val": str(2 * math.pi * 60)});
    dict4xml["params"] = np.append(dict4xml["params"], {"name": "ws", "val": str(1)});
    dict4xml["params"] = np.append(dict4xml["params"], {"name": "pi", "val": str(math.pi)});

    for i in range(nsg):
        modelOrder = ppc["sg"][i, 4];
        ibus = int(ppc["sg"][i, 0]);
        # VARIABLES
        # states
        for j in range(int(modelOrder)):
            dict4xml["vars"] = np.append(dict4xml["vars"], {"name": statevars[j] + str(ibus), "val": statevals[j]});
        # algebraic
        # id
        dict4xml["vars"] = np.append(dict4xml["vars"], {"name": "id" + str(ibus), "val": str(0)});
        # iq
        dict4xml["vars"] = np.append(dict4xml["vars"], {"name": "iq" + str(ibus), "val": str(1)});

        # PARAMETERS:
        # w_0
        dict4xml["params"] = np.append(dict4xml["params"], {"name": "w" + str(ibus) + "_0", "val": "ws"});
        # ra
        dict4xml["params"] = np.append(dict4xml["params"], {"name": "ra" + str(ibus), "val": str(ppc["sg"][i, 6])});
        # x1d
        dict4xml["params"] = np.append(dict4xml["params"], {"name": "x1d" + str(ibus), "val": str(ppc["sg"][i, 7])});
        # M
        dict4xml["params"] = np.append(dict4xml["params"], {"name": "M" + str(ibus), "val": str(ppc["sg"][i, 17])});
        # D
        dict4xml["params"] = np.append(dict4xml["params"], {"name": "D" + str(ibus), "val": str(ppc["sg"][i, 18])});
        # Pm_0
        dict4xml["params"] = np.append(dict4xml["params"], {"name": "Pm" + str(ibus) + "_0", "val": str(1), "out": "true"});
        # E
        dict4xml["params"] = np.append(dict4xml["params"], {"name": "E" + str(ibus), "val": str(1), "out": "true"});
        if modelOrder > 2:
            pass
        
        # mechanical torque/power
        if modelOrder == 2:
            pmEq = "E" + str(ibus) + " * iq" + str(ibus); 
            #+ " * cos(delta" + str(ibus) + "_0) - id" + str(ibus) + " * sin(delta" + str(ibus) + "_0))";
        else: 
            pass   

        # active & reactive power
        if modelOrder == 2:
            pgEq = "V" + str(ibus) + " * (iq" + str(ibus) + " * cos(theta" + str(ibus) + " - delta" + str(ibus) + ")"\
            " - id" + str(ibus) + " * sin(theta" + str(ibus) + " - delta" + str(ibus) + "))";
            qgEq = "V" + str(ibus) + " * (id" + str(ibus) + " * cos(theta" + str(ibus) + " - delta" + str(ibus) + ")"\
            " + iq" + str(ibus) + " * sin(theta" + str(ibus) + " - delta" + str(ibus) + "))";
        else:
            pass;      

        # ODEqs
        # delta Eq
        dEq = "delta" + str(ibus) + "' = wn * (" + "w" + str(ibus) + " - ws)";
        dict4xml["odes"] = np.append(dict4xml["odes"], {"fx": dEq});
        # w Eq
        wEq = "w" + str(ibus) + "' = (1/M" + str(ibus) + ") * (Pm" + str(ibus) + "_0" +  " -(" + pmEq + ") - D" + str(ibus) + "* (w" + str(ibus) + " - ws))";
        dict4xml["odes"] = np.append(dict4xml["odes"], {"fx": wEq});
        if modelOrder > 2:
            pass

        # NLEqs
        # real and imaginary KVL equations:
        reKVLEq = "V" + str(ibus) + "*cos(theta" + str(ibus) + ") - E" + str(ibus) + " * cos(delta" + str(ibus) +\
            ") - ra" + str(ibus) + "*iq" + str(ibus) + "*cos(delta" + str(ibus) + ")"\
            " - ra" + str(ibus) + "*id" + str(ibus) + "*sin(delta" + str(ibus) + ")"\
            " + x1d" + str(ibus) + "*iq" + str(ibus) + "*sin(delta" + str(ibus) + ")"\
            " - x1d" + str(ibus) + "*id" + str(ibus) + "*cos(delta" + str(ibus) + ") = 0";
        dict4xml["nleqs"] = np.append(dict4xml["nleqs"], {"fx": reKVLEq});
        imKVLEq = "V" + str(ibus) + "*sin(theta" + str(ibus) + ") - E" + str(ibus) + " * sin(delta" + str(ibus) +\
            ") - ra" + str(ibus) + "*iq" + str(ibus) + "*sin(delta" + str(ibus) + ")"\
            " + ra" + str(ibus) + "*id" + str(ibus) + "*cos(delta" + str(ibus) + ")"\
            " - x1d" + str(ibus) + "*iq" + str(ibus) + "*cos(delta" + str(ibus) + ")"\
            " - x1d" + str(ibus) + "*id" + str(ibus) + "*sin(delta" + str(ibus) + ") = 0";
        dict4xml["nleqs"] = np.append(dict4xml["nleqs"], {"fx": imKVLEq});
        # active power injected into ibus node
        pEqStr = "V" + str(ibus) + "^2*G" + str(ibus) + "_" + str(ibus); 
        for j in range(nb):
            if j != i and busesconn(ppc, i + 1, j + 1):
                pEqStr = pEqStr + "+ V" + str(ibus) + "*V" + str(j + 1) + "*(G" + str(ibus) + "_" + str(j + 1) + \
                "*cos(theta" + str(ibus) + "- theta" + str(j + 1) + ") + B" + str(ibus) + "_" + str(j + 1) + "*sin(theta" + \
                str(ibus) + " - theta" + str(j + 1) + "))";
        pEqStr = pEqStr + " = " + pgEq  + " - pl" + str(ibus); 
        dict4xml["nleqs"] = np.append(dict4xml["nleqs"], {"fx": pEqStr});
        # reactive power injected into ibus node
        qEqStr = "-V" + str(ibus) + "^2*B" + str(ibus) + "_" + str(ibus); 
        for j in range(nb):
            if j != (ibus - 1) and busesconn(ppc, i + 1, j + 1):
                qEqStr = qEqStr + "+ V" + str(ibus) + "*V" + str(j + 1) + "*(G" + str(ibus) + "_" + str(j + 1) + \
                    "*sin(theta" + str(ibus) + "- theta" + str(j + 1) + ") - B" + str(ibus) + "_" + str(j + 1) + "*cos(theta" + \
                    str(i + 1) + " - theta" + str(j + 1) + "))";
        qEqStr = qEqStr + " = " + qgEq +  "- ql" + str(ibus);
        dict4xml["nleqs"] = np.append(dict4xml["nleqs"], {"fx": qEqStr});

        ####################################### INITIALIZATION #########################################################
        # Parameters
        # imag injected current
        dict4xml["init"]["params"] = np.append(dict4xml["init"]["params"], {"name": "E" + str(ibus) + "r", "val": str(0)});
        dict4xml["init"]["params"] = np.append(dict4xml["init"]["params"], {"name": "E" + str(ibus) + "im", "val": str(0)});
        # delta_0
        dict4xml["init"]["params"] = np.append(dict4xml["init"]["params"], {"name": "delta" + str(ibus) + "_0", "val": str(0)});

        # PostProcessing
        # base w
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": "base.w" + str(ibus) + " = ws"});
        # Ere
        eReStr = "E" + str(ibus) + "r = V" + str(ibus) +" * cos(theta" + str(ibus) + ") + ra" + str(ibus) + \
        "* I" + str(ibus) + "r - x1d" + str(ibus) + "* I" + str(ibus) + "im"; 
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": eReStr});
        # Eim
        eImStr = "E" + str(ibus) + "im = V" + str(ibus) +" * sin(theta" + str(ibus) + ") + ra" + str(ibus) + \
        "* I" + str(ibus) + "im + x1d" + str(ibus) + "* I" + str(ibus) + "r"; 
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": eImStr});
        # base E0
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": "base.E" + str(ibus) + " = sqrt(E" + \
                                                                   str(ibus) + "r^2 + E" + str(ibus) + "im^2)"});

        # base delta
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": "delta" + str(ibus) + "_0 = atg(E" + \
                                                                    str(ibus) + "im/E" + str(ibus) + "r)"});
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": "base.delta" + str(ibus) + " = delta" + str(ibus) + "_0"});
        # id
        idStr = "I" + str(ibus) + "r * cos(pi/2 - delta" + str(ibus) + "_0) - I" + str(ibus) + "im * sin(pi/2 - delta" + str(ibus) + "_0)"
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": "base.id" + str(ibus) + " = " + idStr});
        # iq 
        iqStr =  "I" + str(ibus) + "r * sin(pi/2 - delta" + str(ibus) + "_0) + I" + str(ibus) + "im * cos(pi/2 - delta" + str(ibus) + "_0)"
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": "base.iq" + str(ibus) + " = " + iqStr});
        # base pm
        dict4xml["init"]["pproc"] = np.append(dict4xml["init"]["pproc"], {"fx": "base.Pm" + str(ibus) + "_0 = E" + str(ibus)\
                                                                           + " * (" + str(iqStr) + ")"});
        ################################################################################################################

def xml4DAE(settings, ppc):

    # Define dictioneries for nodes in XML
    dict4xml = { 
                "vars": np.empty(0, dtype=object), "params": np.empty(0, dtype=object), 
                #################################### initialization ####################################
                "init": 
                { 
                "vars": np.empty(0, dtype=object), "params": np.empty(0, dtype=object),\
                "nleqs": np.empty(0, dtype=object), "pproc": np.empty(0, dtype=object)
                },
                ########################################################################################
                "odes": np.empty(0, dtype=object), "nleqs": np.empty(0, dtype=object)
                };
    
    # POWER NETWORK XML configuration
    power_network(settings, ppc, dict4xml); 

    # SYNCHRONOUS GENERATOR XML configuration
    synchronous_generator(settings, ppc, dict4xml);

    # MODEL SOLVER
    model = ET.Element("Model", attrib={"type": "DAE", "domain": "real", "method": settings["nimethod"], "eps": str(settings["EPS"]), "name": settings["filename"]});

    # VARIABLES
    vars = ET.SubElement(model, "Vars", attrib= {"out": "true"})
    nvars = dict4xml["vars"].shape[0];
    for i in range(nvars):
        ET.SubElement(vars, "Var", attrib= dict4xml["vars"][i]);
    
    # PARAMETERS
    params = ET.SubElement(model, "Params")
    nparams = dict4xml["params"].shape[0];
    for i in range(nparams):
        ET.SubElement(params, "Param", attrib= dict4xml["params"][i]);
    
    # INITIALIZATION
    init = ET.SubElement(model, "Init");
    imodel = ET.SubElement(init, "Model", attrib={"type": "NR", "domain": "real", "eps" : str(1e-6), "name": "PF Subproblem for DAE"});

    # I_Vars
    ivars = ET.SubElement(imodel, "Vars")
    nivars = dict4xml["init"]["vars"].shape[0];
    for i in range(nivars):
        ET.SubElement(ivars, "Var", attrib= dict4xml["init"]["vars"][i]); 
    
    # I_PARAMETERS
    iparams = ET.SubElement(imodel, "Params")
    niparams = dict4xml["init"]["params"].shape[0];
    for i in range(niparams):
        ET.SubElement(iparams, "Param", attrib= dict4xml["init"]["params"][i]);
    
    # I_NLeqs
    inleqs = ET.SubElement(imodel, "NLEqs")
    ninleqs = dict4xml["init"]["nleqs"].shape[0];
    for i in range(ninleqs):
        ET.SubElement(inleqs, "Eq", attrib= dict4xml["init"]["nleqs"][i]);
    
    # I_Pproc
    ipproc = ET.SubElement(imodel, "PostProc")
    nipproc = dict4xml["init"]["pproc"].shape[0];
    for i in range(nipproc):
        ET.SubElement(ipproc, "Eq", attrib= dict4xml["init"]["pproc"][i]);

    # ODEqs
    odes = ET.SubElement(model, "ODEqs")
    nodes = dict4xml["odes"].shape[0];
    for i in range(nodes):
        ET.SubElement(odes, "Eq", attrib= dict4xml["odes"][i]);
    
    # NLEqs
    nleqs = ET.SubElement(model, "NLEqs")
    nnleqs = dict4xml["nleqs"].shape[0];
    for i in range(nnleqs):
        ET.SubElement(nleqs, "Eq", attrib= dict4xml["nleqs"][i]);

    # Convert the ElementTree to a string
    xml_str = ET.tostring(model, encoding='utf-8', method='xml')

    # Use minidom to pretty-print the XML string
    pretty_xml = minidom.parseString(xml_str).toprettyxml(indent="    ")

    with open("modelSolver/real/" + settings["filename"] + ".xml", "w", encoding='utf-8') as files:
        files.write('<?xml version="1.0" encoding="UTF-8"?>\n' + pretty_xml.split('<?xml version="1.0" ?>\n', 1)[1])