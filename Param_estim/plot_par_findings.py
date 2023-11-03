import matplotlib.pyplot as plt
import pandas as pd
import csv
from instructions import *

param_values = []
stor_fold = storage_folder()
file = os.path.join(stor_fold, 'log_book.txt')
with open(file, "r") as f:
    reader = csv.reader(f, delimiter='\t')
    param_names = next(reader)  # Read the first row as parameter names
    for row in reader:
        param_values.append([float(x) for x in row])

# Transpose the list so that each parameter has its own list of values
param_values = list(map(list, zip(*param_values)))

# Plot each parameter's values over time
import matplotlib.pyplot as plt

fig, axs = plt.subplots(len(param_names), 1, sharex=True)

for i, param in enumerate(param_values):
    axs[i].plot(param)
    axs[i].set_ylabel(str(param_names[i]))

plt.xlabel("Iteration")
plt.show()
