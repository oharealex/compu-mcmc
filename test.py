import matplotlib.pyplot as plt
import numpy as np
import random as r

ydata = []
for i in range (100):
	ydata.append(r.randint(0, 50))

#plt.scatter(np.arange(len(ydata)), ydata, color = "red")
plt.errorbar(np.arange(len(ydata)), ydata, fmt = "-o", yerr = 5, color = "black", markersize = 5, ecolor = "red", capsize = 2)
plt.show()
