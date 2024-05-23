# Lennard-Jones Potential
# Created by Alex O'Hare, 28/04/2024

# Import the required libraries
import matplotlib.pyplot as plt
import numpy as np
import random as r

# Define the Lennard-Jones potential equation
def LJ(r, epsilon, sigma):
	return (4 * epsilon) * (((sigma / r) ** 12) - ((sigma / r) ** 6))

def rep(r, epsilon, sigma):
	return (4 * epsilon) * ((sigma / r) ** 12)

def att(r, epsilon, sigma):
	return - (4 * epsilon) * ((sigma / r) ** 6)
	
epsilon = 0.4396
sigma = 0.3851

# Evaluate the LJ potential for a range of distances 0 - 4
distances = np.linspace(0.9 * sigma, 3 * sigma, 1000000)
phi = [LJ(dis, epsilon, sigma) for dis in distances]

# Plot the LJ potential 
plt.plot(distances, phi, color = "black", linewidth = 0.9)
plt.xlabel("r ($\AA$)")
plt.ylabel("$\phi$(r) (eV)")
plt.text(0.38, 2, r"$\phi_{LJ}(r) = 4\varepsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]$")

# Plot the details
## r_min
plt.scatter(distances[phi.index(min(phi))], 0, color = "black", s = 15	)
plt.text(distances[phi.index(min(phi))], 0, "$r_{min}$", 
		va = "bottom", color = "black")

## r = sigma
zero_crossing = np.interp(0, np.flip(phi), np.flip(distances))
plt.scatter(zero_crossing, 0, color = "black", s = 15)
plt.text(zero_crossing-0.08, 0, "$r = \sigma$", va = "bottom", 
		color = "black")
		
## epsilon
plt.plot([distances[phi.index(min(phi))], distances[phi.index(min(phi))]], 
		[0, min(phi)], linestyle = "--", color = "blue", linewidth = 0.9)
plt.text(0.44, -0.2, r"$\varepsilon$")

## repulsion
A = [rep(dis, epsilon, sigma) for dis in distances]
plt.plot(distances, A, color = "red", linewidth = 0.9, linestyle = "--")
plt.text(0.41, 0.9, r"$+4\varepsilon\left(\frac{\sigma}{r}\right)^{12}$")
## attraction
B = [att(dis, epsilon, sigma) for dis in distances]
plt.plot(distances, B, color = "green", linewidth = 0.9, linestyle = "--")
plt.text(0.52, -0.5, r"$-4\varepsilon\left(\frac{\sigma}{r}\right)^6$")
## x axis
plt.axhline(0, color = "black", linewidth = 0.9, linestyle = "--")

## plot settings
plt.ylim(-0.7, 3.0)
plt.xlim(0.2, 1.1)
plt.savefig("carbon-carbon_LJ.png")

