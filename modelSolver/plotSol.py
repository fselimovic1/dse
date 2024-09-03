#!/usr/bin/env python3
#(c) Izudin Dzafic, 2023
import sys
import numpy as np
import matplotlib.pyplot as plt

# varsToPlot = ["w1", "w2", "theta1", "theta2", "theta3", "V3" ]
varsToPlot = [ "w1", "w2" ]
#varsToPlot = [ "w1", "w2", "w3", "w4", "w5", "w6", "w7", "w8", "w9", "w10" ]
# varsToPlot = [ "w1", "w2", "w3" ]
args = len(sys.argv)-1
if args != 1:
   print('Missing input file name!')
   sys.exit()
    
# set dark background
#plt.style.use('dark_background')
# Read the file and extract the relevant lines
lines = []
header = []

with open(sys.argv[1], 'r') as file:
    foundSolutionData = False
    headerLoaded = False
    firstHeaderLine = False
    secondHeaderLine = False
    for line in file:
        if foundSolutionData and firstHeaderLine and not headerLoaded:
            header = line.strip().split()
            headerLoaded = True
            print('Header is loaded')
        elif foundSolutionData:
            if not firstHeaderLine:
                if line.startswith('------'):
                    firstHeaderLine = True
                    continue
                else:
                    print('File is not formated properly! First line (--------) is not available!')
                    sys.exit()
            if not secondHeaderLine:       
                if line.startswith('------'):
                    secondHeaderLine = True
                    continue
                else:
                    print('File is not formated properly! Second line (--------) is not available!')
                    sys.exit()             
            if line.startswith('----'):
                print('Detected end of SOLUTION DATA')
                break
            lines.append(line)
        elif line.strip() == 'SOLUTION DATA':
            foundSolutionData = True
            print('Detected SOLUTION DATA')

# Parse the data
data = np.loadtxt(lines, dtype=float)

# Check if data is empty
if data.size == 0:
    print("No data found in the file.")
else:
    # Extract the columns
    t = data[:, 0]
    columns = data[:, 1:].T
    
    # Plot the data
    for i, column in enumerate(columns):
        if header[i+1] in varsToPlot:
            plt.plot(t, column, label=header[i+1])

    # Add labels and legend
    if len(header) <= 11:
        xLbl = ''
        yLbl = ''
        k=0
        for lbl in header:
            if k>1:
                yLbl += ', ' + lbl.strip()
            elif k>0:
                yLbl = lbl.strip()
            else:
                xLbl = lbl.strip()
            k=k+1
        plt.xlabel(xLbl)
        plt.ylabel(yLbl)     
    else:
        plt.xlabel('Time (t) [s]')
        plt.ylabel('Multiple Values')
        
    plt.legend()
    
    # show grid
    plt.grid(linestyle='dotted', linewidth=1)
    
    # Display the plot
    plt.show()