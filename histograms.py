# Code for calculating mean and variance from energy file

# Import the required libraries
import numpy as np
import matplotlib.pyplot as plt

# Import the data file
energy_data = np.loadtxt("energy.txt")
bins_range = np.linspace(np.min(energy_data), np.max(energy_data), 10)
# Plot the data
plt.hist(energy_data, bins = bins_range, color = "black", edgecolor = "grey")
plt.xlabel("Energy")
plt.ylabel("Frequency")
plt.show()
