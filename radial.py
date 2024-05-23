import matplotlib.pyplot as plt
import numpy as np

rho = 0.1
def N(rho, r):
	return 4 * np.pi * rho * r * r
	
values = np.arange(1, 3, (3 - 1) / 100)
answers = np.arange(100)
for val in range(len(values)):
	answers[val] = N(rho, val)

plt.plot(values, answers, color = "black")
plt.xlabel("r")
plt.ylabel("N(r)")
#plt.show()
plt.close()

def g(rho, val):
	return (val
