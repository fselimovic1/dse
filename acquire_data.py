import numpy as np

def acquire_data(fileName):
    solution_data_found = False
    count = 0
    first_row = True
    is_data_line = False
    with open(fileName, 'r') as file:
        for line in file:
            if solution_data_found and count  == 1 and line.strip() != "---------------------------":
                vars = line.split()
            if solution_data_found and count == 2:
                is_data_line = True
            else:
                is_data_line = False
            if is_data_line and line.strip() != "------------------------":
                if first_row:
                    data = np.array(line.split(), dtype=float);
                    first_row = False;
                else:
                    data = np.vstack((data, np.array(line.split(), dtype=float)));
            if line.strip() == 'SOLUTION DATA':
                solution_data_found = True;
            if solution_data_found and (line.strip() == "---------------------------" or line.strip() == "------------------------"):
                count += 1;
    sim_data = dict();
    for i in range(0,len(vars)):
        sim_data[vars[i]] = data[:, i].tolist()
    return sim_data

