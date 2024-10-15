import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
import numpy as np
import math
import ybus
import random

EPS = 1e-6
ms = 1e-3

rf = 1e-2;

def xml4DAE_CMNNE(settings, ppc):
    
    # Power System params
    nb = ppc["bus"].shape[0]
    ng = ppc["gen"].shape[0]
    genIdx = (ppc["gen"][:, 0] - 1).astype(int)
    
    # Syn Gen params
    Xd = ppc["sg"][:, 8]
    M = ppc["sg"][:, 17]
    D = ppc["sg"][:, 18]
    ws = 1
    wn = 2 * math.pi * ppc["fn"]

    # Power Network
    Yr, Ynn, Yss, Yns = ybus.reduced_ybus(ppc);
    Gr = np.real(Yr);
    Br = np.imag(Yr);
    Yv = -np.linalg.inv(Yss) @ np.transpose(Yns);
    Gv = np.real(Yv);
    Bv = np.imag(Yv);

    # INITIALITATION PHASE
    # Voltages
    V = ppc["bus"][:, 7];
    theta = ppc["bus"][:, 8];
    v = V * np.exp(1j * theta);
    Vr = np.real(v);
    Vim = np.imag(v);

    # Powers
    Pg = ppc["gen"][:, 1];
    Qg = ppc["gen"][:, 2];
    Eg = V[genIdx] + Qg * Xd / V[genIdx] + 1j * Pg * Xd / V[genIdx];
    delta = theta[genIdx] + np.angle(Eg);
    E = abs(Eg)
    Eg = E * np.exp(1j * delta)
    Ig = Yr @ Eg
    Sg = Eg * np.conj(Ig);
    Pg = np.real(Sg);
    Qg = np.imag(Sg);
    # Mechanical input power
    Pm = Pg;

    # MODEL SOLVER
    model = ET.Element("Model", attrib={"type": "DAE", "domain": "real", "method": settings["nimethod"], "eps": str(settings["EPS"]), "name": settings["filename"]});

    # VARIABLES
    vars = ET.SubElement(model, "Vars", attrib= {"out": "true"})
    for i in range(ng):
        ET.SubElement(vars, "Var", attrib={"name": "delta" + str(i + 1), "val": str(delta[i])})
        ET.SubElement(vars, "Var", attrib={"name": "w" + str(i + 1), "val": str(ws)})

    if settings["intType"] == "simultaneous":
        for i in range(nb):
            ET.SubElement(vars, "Var", attrib={"name": "Vr" + str(i + 1), "val": str(Vr[i])});
            ET.SubElement(vars, "Var", attrib={"name": "Vim" + str(i + 1), "val": str(Vim[i])}); 

    # PARAMETERS
    params = ET.SubElement(model, "Params")#, attrib= {"dT": str(ms)})
    ET.SubElement(params, "Param", attrib={"name": "ws", "val": str(ws)});
    ET.SubElement(params, "Param", attrib={"name": "wn", "val": str(wn)});

    for i in range(ng):
        ET.SubElement(params, "Param", attrib={"name": "x1d" + str(i + 1), "val": str(Xd[i])}); 
        ET.SubElement(params, "Param", attrib={"name": "M" + str(i + 1), "val": str(M[i])});
        ET.SubElement(params, "Param", attrib={"name": "D" + str(i + 1), "val": str(D[i])}); 
        ET.SubElement(params, "Param", attrib={"name": "Pm" + str(i + 1), "val": str(Pm[i]), "out": "true"});
        ET.SubElement(params, "Param", attrib={"name": "E" + str(i + 1), "val": str(E[i]), "out": "true"});    
        ET.SubElement(params, "Param", attrib={"name": "pg" + str(i + 1), "val": str(0), "out": "true"});
        ET.SubElement(params, "Param", attrib={"name": "qg" + str(i + 1), "val": str(0), "out": "true"});  

    # Yr 
    for i in range(ng):
        for j in range(ng):
            ET.SubElement(params, "Param", attrib={"name": "Gr" + str(i + 1) + "_"  + str(j + 1), "val": str(Gr[i][j])}); 
            ET.SubElement(params, "Param", attrib={"name": "Br" + str(i + 1) + "_" + str(j + 1), "val": str(Br[i][j])}); 
    # Yv    
    for i in range(nb):
        for j in range(ng):
            ET.SubElement(params, "Param", attrib={"name": "Gv" + str(i + 1) + "_" +  str(j + 1), "val": str(Gv[i][j])}); 
            ET.SubElement(params, "Param", attrib={"name": "Bv" + str(i + 1) + "_" + str(j + 1), "val": str(Bv[i][j])}); 
    
    if settings["intType"] == "partitioned":
        for i in range(nb):
            ET.SubElement(params, "Param", attrib={"name": "Vr" + str(i + 1), "val": str(Vr[i])});
            ET.SubElement(params, "Param", attrib={"name": "Vim" + str(i + 1), "val": str(Vim[i])}); 
    
    for i in range(nb):
        ET.SubElement(params, "Param", attrib={"name": "V" + str(i + 1), "val": str(V[i]), "out": "true"});
        ET.SubElement(params, "Param", attrib={"name": "theta" + str(i + 1), "val": str(theta[i]), "out": "true"});
    
    # Ordinary differential equations
    odeqs = ET.SubElement(model, "ODEqs")
    for i in range(ng):
        wODE = "w" + str(i + 1) + "' = 1/M" + str(i + 1) + " * (Pm" + str(i + 1) +  " - E"\
                + str(i + 1) + "^2 * Gr" + str(i + 1) + "_" + str(i + 1) + " - D" + str(i+1) + " * (w" + str(i + 1)\
                + " - ws)";
        for j in range(ng):
            if i == j:
                continue
            wODE += " - " + "E" + str(i + 1) + " * E" + str(j + 1) + " * (Br" + str(i + 1) + "_" + str(j + 1) \
            + " * sin(delta" + str(i + 1) + " - delta" + str(j + 1) + ") + Gr" + str(i + 1) + "_" + str(j + 1) \
            + " * cos(delta" + str(i + 1) + " - delta" + str(j + 1) + "))" 
        wODE += ")";
        dODE = "delta" + str(i + 1) + "' = wn * (w" + str(i + 1) + " - ws)";
        ET.SubElement(odeqs, "Eq", attrib={"fx": wODE});  
        ET.SubElement(odeqs, "Eq", attrib={"fx": dODE}); 
    
    # Nonlinear equations
    nleqs = ET.SubElement(model, "NLEqs")
    if settings["intType"] == "simultaneous": 
        for i in range(nb):
            vrNLE = "Vr" + str(i + 1) + " =";
            vimNLE = "Vim" + str(i + 1) + " =";
            for j in range(ng):
                vrNLE += " + E" + str(j + 1) + " * (Gv" + str(i + 1) + "_" + str(j + 1) + " * cos(delta" \
                    + str(j + 1) + ") - Bv" + str(i + 1) + "_" + str(j + 1) + " * sin(delta" + str(j + 1) + "))";
                vimNLE += " + E" + str(j + 1) + " * (Gv" + str(i + 1) + "_" + str(j + 1) + " * sin(delta" \
                    + str(j + 1) + ") + Bv" + str(i + 1) + "_" + str(j + 1) + " * cos(delta" + str(j + 1) + "))";

            ET.SubElement(nleqs, "Eq", attrib={"fx": vrNLE});  
            ET.SubElement(nleqs, "Eq", attrib={"fx": vimNLE}); 

    # Post-processing -> disturbance generation
    # POST-PROCESSING
    pproc = ET.SubElement(model, "PostProc")
    if settings["intType"] == "partitioned":
        for i in range(nb):
            vrNLE = "Vr" + str(i + 1) + " =";
            vimNLE = "Vim" + str(i + 1) + " =";
            for j in range(ng):
                vrNLE += " + E" + str(j + 1) + " * (Gv" + str(i + 1) + "_" + str(j + 1) + " * cos(delta" \
                    + str(j + 1) + ") - Bv" + str(i + 1) + "_" + str(j + 1) + " * sin(delta" + str(j + 1) + "))";
                vimNLE += " + E" + str(j + 1) + " * (Gv" + str(i + 1) + "_" + str(j + 1) + " * sin(delta" \
                    + str(j + 1) + ") + Bv" + str(i + 1) + "_" + str(j + 1) + " * cos(delta" + str(j + 1) + "))";

            ET.SubElement(pproc, "Eq", attrib={"fx": vrNLE});  
            ET.SubElement(pproc, "Eq", attrib={"fx": vimNLE});

    # Compute bus voltages in polar coordinates    
    for i in range(nb):
        ET.SubElement(pproc, "Eq", attrib={"fx": "V" + str(i + 1) + " = sqrt(Vr" + str(i + 1) + "^2 + " + "Vim" + str(i + 1) + "^2)" });
        ET.SubElement(pproc, "Eq", attrib={"fx": "theta" + str(i + 1) + " = atg(Vim" + str(i + 1) +"/" + "Vr" + str(i + 1) + ")"});
    
    # Compute active and reactive generator powers
    for i in range(ng):
        eqPg = "pg" + str(i + 1) + "= E" + str(i + 1) + "^2 * Gr" + str(i + 1) + "_" + str(i + 1);
        eqQg = "qg" + str(i + 1) + "= -E" + str(i + 1) + "^2 * Br" + str(i + 1) + "_" + str(i + 1);    
        for j in range(ng):
            if i == j:
                continue;
            eqPg += "+" +  "E" + str(i + 1) + " * E" + str(j + 1) + " * (Br" + str(i + 1) + "_" + str(j + 1) \
            + " * sin(delta" + str(i + 1) + " - delta" + str(j + 1) + ") + Gr" + str(i + 1) + "_" + str(j + 1) \
            + " * cos(delta" + str(i + 1) + " - delta" + str(j + 1) + "))" 
            eqQg += "+" +  "E" + str(i + 1) + " * E" + str(j + 1) + " * (-Br" + str(i + 1) + "_" + str(j + 1) \
            + " * cos(delta" + str(i + 1) + " - delta" + str(j + 1) + ") + Gr" + str(i + 1) + "_" + str(j + 1) \
            + " * sin(delta" + str(i + 1) + " - delta" + str(j + 1) + "))" 
        ET.SubElement(pproc, "Eq", attrib={"fx": eqPg });
        ET.SubElement(pproc, "Eq", attrib={"fx": eqQg });

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

        # Event network model
        Yr_e, Ynn_e, Yss_e, Yns_e = ybus.reduced_ybus(ppc)
        Gr_e = np.real(Yr_e);
        Br_e = np.imag(Yr_e);
        Yv_e = -np.linalg.inv(Yss_e) @ np.transpose(Yns_e);
        Gv_e = np.real(Yv_e);
        Bv_e = np.imag(Yv_e);
        ppc["bus"][bus_toadd, 2] -= load_toadd;
    
        # post event matrices
        Gr_c = Gr;
        Br_c = Br;
        Gv_c = Gv;
        Bv_c = Bv;
    
    elif etype == "bbfault":
        # event matrices
        nbus = int(settings["event"]["noBus"] - 1);
        Yss_e = np.copy(Yss);
        Yss_e[nbus, :] = np.zeros(nb);
        Yss_e[:, nbus] = np.zeros(nb);
        Yr_e = Ynn - Yns @ np.linalg.pinv(Yss_e) @ np.transpose(Yns);
        Gr_e = np.real(Yr_e);
        Br_e = np.imag(Yr_e);
        Yv_e = -np.linalg.pinv(Yss_e) @ np.transpose(Yns);
        Gv_e = np.real(Yv_e);
        Bv_e = np.imag(Yv_e);
        # post event matrices
        Gr_c = Gr;
        Br_c = Br;
        Gv_c = Gv;
        Bv_c = Bv;
    
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
        Yr_c, Ynn_c, Yss_c, Yns_c = ybus.reduced_ybus(ppc);
        Gr_c = np.real(Yr_c);
        Br_c = np.imag(Yr_c);
        Yv_c = -np.linalg.inv(Yss_c) @ np.transpose(Yns_c);
        Gv_c = np.real(Yv_c);
        Bv_c = np.imag(Yv_c);


        # EVENT matrices
        Yss_e = np.copy(Yss_c);
        Yss_e[frombus, tobus] -=  1/(zd_a);
        Yss_e[tobus, frombus] -=  1/(zd_a);
        Yss_e[frombus, frombus] += (1/(zd_s) + 1j * bs/2);
        Yss_e[tobus, tobus] += (1/(zd_s) + 1j * bs/2);
        Yr_e = Ynn_c - Yns_c @ np.linalg.inv(Yss_e) @ np.transpose(Yns_c);
        Gr_e = np.real(Yr_e);
        Br_e = np.imag(Yr_e);
        Yv_e = -np.linalg.inv(Yss_e) @ np.transpose(Yns_c);
        Gv_e = np.real(Yv_e);
        Bv_e = np.imag(Yv_e);

        ppc["branch"] = backupBranch
    
    # Event generation
    # Yr  & Yv matrices
    for i in range(ng):
        for j in range(ng):
            cond = ET.SubElement(pproc, "Eq", attrib={"cond": "t > " + str(settings["etime_s"])});
            ET.SubElement(cond, "Then", attrib={"fx": "Gr" + str(i + 1) + "_" + str(j + 1) + " = " + str(Gr_e[i][j])});
            cond = ET.SubElement(pproc, "Eq", attrib={"cond": "t > " + str(settings["etime_s"])}); 
            ET.SubElement(cond, "Then", attrib={"fx": "Br" + str(i + 1) + "_" + str(j + 1)  + " = " + str(Br_e[i][j])});  
        
    for i in range(nb):
        for j in range(ng):
            cond = ET.SubElement(pproc, "Eq", attrib={"cond": "t > " + str(settings["etime_s"])});
            ET.SubElement(cond, "Then", attrib={"fx": "Gv" + str(i + 1) + "_" + str(j + 1) + " = " + str(Gv_e[i, j])}); 
            cond = ET.SubElement(pproc, "Eq", attrib={"cond": "t > " + str(settings["etime_s"])});
            ET.SubElement(cond, "Then", attrib={"fx": "Bv" + str(i + 1) + "_" + str(j + 1) + " = " + str(Bv_e[i, j])});
    
    # Yr  & Yv matrices
    for i in range(ng):
        for j in range(ng):
            cond = ET.SubElement(pproc, "Eq", attrib={"cond": "t > " + str(settings["etime_e"])});
            ET.SubElement(cond, "Then", attrib={"fx": "Gr" + str(i + 1) + "_" + str(j + 1) + " = " + str(Gr_c[i][j])}); 
            cond = ET.SubElement(pproc, "Eq", attrib={"cond": "t > " + str(settings["etime_e"])});
            ET.SubElement(cond, "Then", attrib={"fx": "Br" + str(i + 1) + "_" + str(j + 1)  + " = " + str(Br_c[i][j])});  
    for i in range(nb):
        for j in range(ng):
            cond = ET.SubElement(pproc, "Eq", attrib={"cond": "t > " + str(settings["etime_e"])});
            ET.SubElement(cond, "Then", attrib={"fx": "Gv" + str(i + 1) + "_" + str(j + 1) + " = " + str(Gv_c[i, j])});
            cond = ET.SubElement(pproc, "Eq", attrib={"cond": "t > " + str(settings["etime_e"])});
            ET.SubElement(cond, "Then", attrib={"fx": "Bv" + str(i + 1) + "_" + str(j + 1) + " = " + str(Bv_c[i, j])});
    
    # Convert the ElementTree to a string
    xml_str = ET.tostring(model, encoding='utf-8', method='xml')

    # Use minidom to pretty-print the XML string
    pretty_xml = minidom.parseString(xml_str).toprettyxml(indent="    ")

    with open("modelSolver/real/" + settings["filename"] + ".xml", "w", encoding='utf-8') as files:
        files.write('<?xml version="1.0" encoding="UTF-8"?>\n' + pretty_xml.split('<?xml version="1.0" ?>\n', 1)[1])