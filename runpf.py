import os
from routines.power_flows.xml4PF import xml4PF
from process_data import process_data
import numpy as np

from power_systems import case118
ppc = process_data(case118.case118())

outputFile = "ieee_pf_results"
# POWER FLOWS settings
settings = {"filename": "ieee118_pf"}
settings["eps"] = 1e-10;

# generate xml file for PF analysis
xml4PF(settings, ppc) 

# Run power flow computations using modelSolver
os.chdir("modelSolver");
msCommand = "modelSolver NR real real\\" + settings["filename"] + ".xml" + " ..\\" + outputFile;
os.system(msCommand);