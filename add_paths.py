import os
import sys

# Get the current working directory
current_dir = os.getcwd()

# Add the current directory to sys.path
sys.path.append(current_dir)

# Walk through the directory tree and add each directory to sys.path
for root, dirs, files in os.walk(current_dir):
    for dir in dirs:
        dir_path = os.path.join(root, dir)
        if dir_path not in sys.path:
            sys.path.append(dir_path)
