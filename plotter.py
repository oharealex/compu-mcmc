import matplotlib.pyplot as plt
import numpy as np


data = np.load("energy_his_001.npy")

plt.figure(figsize = (12, 4))
plt.plot(data, color = "black")
plt.xlabel("iterations")
plt.ylabel(r"$\langle e \rangle$")

print("mean across iterations: " + str(np.mean(data)))
print("max across iterations: " + str(np.max(data)))
print("min across iterations: " + str(np.min(data)))
print("error: " + str((np.max(data) - np.min(data)) / 2))

#plt.savefig("sample_energy.png")
#plt.show()
plt.close()

data2 = []
errdata = []
for i in range(len(data)):
	if i <= (len(data) - 500):
		if i % 500 == 0:
			dataset = [data[j] for j in range(i, i + 499)]
			dataval = sum(dataset) / len(dataset)
			data2.append(dataval)
			errdata.append(max(dataset) - min(dataset))

plt.figure(figsize = (12, 4))
plt.errorbar(np.arange(len(data2)), data2, yerr = np.mean(errdata), fmt = "-o",
		color = "black", ecolor = "red", markersize = 5, capsize = 2)
plt.xlabel("iterations")
plt.ylabel(r"$\langle e \rangle$")
#plt.savefig("sample_energy_error.png")
#plt.show()
plt.close()

bins = 25
prob_dens = [[] for _ in range(bins)]
energies = np.arange(np.min(data), np.max(data), (np.max(data) - np.min(data)) / bins)
for i in data:
	closest = np.argmin(np.abs(energies - i))
	prob_dens[closest].append(i)

ydata = np.zeros(bins)
errdata = np.zeros(bins)
for i in range(len(prob_dens)):
	ydata[i] = len(prob_dens[i])
	if sum(prob_dens[i]) > 0:
		errdata[i] = abs((max(prob_dens[i]) - min(prob_dens[i])) / 2)
	else:
		errdata[i] = 0

plt.bar(energies, ydata, color = "black", width = (energies[1] - energies[0]) * 0.8 )
plt.xlabel(r"$\langle e \rangle$")
plt.ylabel(r"$P(\langle e \rangle)$")
plt.show()
