import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
import numpy as np
import math
import ybus
import random

EPS = 1e-6
ms = 1e-3

rf = 1e-2;

def xml4DAE(settings, ppc):

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
    
    # Convert the ElementTree to a string
    xml_str = ET.tostring(model, encoding='utf-8', method='xml')

    # Use minidom to pretty-print the XML string
    pretty_xml = minidom.parseString(xml_str).toprettyxml(indent="    ")

    with open("modelSolver/real/" + settings["filename"] + ".xml", "w", encoding='utf-8') as files:
        files.write('<?xml version="1.0" encoding="UTF-8"?>\n' + pretty_xml.split('<?xml version="1.0" ?>\n', 1)[1])