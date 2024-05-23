import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Temperature values
temp = [1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7, 4.0]

# Pressure values corresponding to different densities
pressure = [0.100000, 0.188889, 0.277778, 0.366667, 0.455556, 0.544444, 0.633333, 0.722222, 0.811111, 0.900000]

# Energy values for each density
energy = [
    [-0.786604, -0.716732, -0.673604, -0.642987, -0.619113, -0.599275, -0.582061, -0.566665, -0.552600, -0.539552],
    [-1.441956, -1.317690, -1.242492, -1.189154, -1.147193, -1.111900, -1.080903, -1.052878, -1.027047, -1.002903],
    [-2.026193, -1.874499, -1.781042, -1.712504, -1.656688, -1.608299, -1.564778, -1.524716, -1.487288, -1.451951],
    [-2.542881, -2.393244, -2.294372, -2.216326, -2.148956, -2.088059, -2.031682, -1.978777, -1.928685, -1.880983],
    [-3.031578, -2.897083, -2.795373, -2.706692, -2.625160, -2.548610, -2.476165, -2.407258, -2.341516, -2.278658],
    [-3.534010, -3.405193, -3.291146, -3.183119, -3.079609, -2.980506, -2.885792, -2.795388, -2.709118, -2.626726],
    [-4.055449, -3.910736, -3.769117, -3.629386, -3.493571, -3.362950, -3.238217, -3.119455, -3.006525, -2.899091],
    [-4.559224, -4.379319, -4.196528, -4.013845, -3.835825, -3.664821, -3.501993, -3.347552, -3.201299, -3.062828],
    [-4.988417, -4.760701, -4.528033, -4.294359, -4.066096, -3.846757, -3.638309, -3.441082, -3.254831, -3.079013],
    [-5.282570, -4.997567, -4.709290, -4.418953, -4.134236, -3.860305, -3.599602, -3.352981, -3.120308, -2.900944]
]

# Convert lists to numpy arrays
temp = np.array(temp)
pressure = np.array(pressure)
energy = np.array(energy)

# Create a meshgrid for temperature and pressure
T, P = np.meshgrid(temp, pressure)

# Create the figure and a 3D Axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
surf = ax.plot_surface(T, P, energy, cmap='viridis')

# Add labels
ax.set_xlabel('Temperature')
ax.set_ylabel('Pressure')
ax.set_zlabel('Energy')

# Add a color bar which maps values to colors
fig.colorbar(surf)

# Show plot
plt.show()

