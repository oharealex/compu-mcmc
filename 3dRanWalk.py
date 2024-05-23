# 3-D random walk
# Created by Alex O'Hare 29/04/2024

# Import the required libraries
import matplotlib.pyplot as plt
import numpy as np
import random as r
import copy

# System settings
ID = "001"
kB = 1.380649 * (10 ** (-23))
temp = 1.3
beta = 1 / (kB * temp)
box_length = 12
box_volume = box_length * box_length * box_length
density = 0.1
num_particles = int(density * box_volume)
step_mag = 2
r_cut = 2.5
print("num_particles "+str(num_particles))

# Particle matrix
part_mat = np.array([[r.randint(0, box_length) for a in range(3)] 
		for b in range(num_particles)])

# Energy matrix
energy_mat = np.zeros(num_particles)
	
# Continuous box
def cont_box(particle):
	modified = False
	for i in range(3):
		if part_mat[particle][i] > box_length:
			part_mat[particle][i] = 0
			modified = True
		elif part_mat[particle][i] < 0:
			part_mat[particle][i] = box_length
			modified = True
	return modified

def random_walk(particle, step_mag):

	dir_axis = r.randint(0, 2)  # Choose a random dimension
	dir_sign = r.choice([-1, 1])  # Choose a random direction
	dir_mag = 1 #r.randint(1, step_mag)
	new_position = part_mat[particle][dir_axis] + dir_sign * dir_mag
	if new_position == 13:
		new_position = 0
	if new_position == 14:
		new_position = 1
	if new_position == 15:
		new_position = 2
	if new_position == -1:
		new_position = 12
	if new_position == -2:
		new_position = 11
	if new_position == -3:
		new_position = 10
	for i in range(num_particles):
		if i != particle:
			x0, y0, z0 = part_mat[particle]
			x1, y1, z1 = part_mat[i]
			dis = distance_calc(x0, x1, y0, y1, z0, z1)
		#	dis = abs(new_position - part_mat[i][dir_axis])
			if dis <= 0: # r_cut:
				pass
			else:
				part_mat[particle][dir_axis] = new_position
					
					
# PDF data storage
pos_1dx_PDF = np.zeros(box_length + 1)
pos_1dy_PDF = np.zeros(box_length + 1)
pos_1dz_PDF = np.zeros(box_length + 1)

# 1D PDF
def PDF(particle):
	x = part_mat[particle][0]
	y = part_mat[particle][1]
	z = part_mat[particle][2]
	pos_1dx_PDF[x] += 1
	pos_1dy_PDF[y] += 1
	pos_1dz_PDF[z] += 1

# Particle-particle distance calculation
def distance_calc(x0, x1, y0, y1, z0, z1):
	return float(np.sqrt(((x0 - x1) * (x0 - x1)) + ((y0 - y1) * (y0 - y1)) + 
			((z0 - z1) * (z0 - z1))))

# Lennard-Jones Potential
def lennard_jones(epsilon, sigma, R):
	if r_cut > R > 0:
		return 4 * epsilon * (((sigma / R) ** 12) - ((sigma / R) ** 6))

# Neighbours list
def old_neighbours(particle):
	neighbours = []
	part0 = particle
	for part1 in range(len(part_mat)):
		if part1 != particle:
			x0, y0, z0 = part_mat[part0]
			x1, y1, z1 = part_mat[part1]
			R = abs(distance_calc(x0, x1, y0, y1, z0, z1))
			if R <= r_cut:
				neighbours.append(part1)
	return neighbours
	
# Energy function
def get_energy_init():
	#energy_mat.fill(0)
	for part0 in range(len(part_mat)):
		for part1 in range(part0 + 1, len(part_mat)):
			x0, y0, z0 = part_mat[part0]
			x1, y1, z1 = part_mat[part1]
			R = abs(distance_calc(x0, x1, y0, y1, z0, z1))
			if r_cut > R > 0:
				E = lennard_jones(1, 1, R)
				energy_mat[part0] += E
				energy_mat[part1] += E

# Energy function
def get_energy(particle, oldies):
	#energy_mat.fill(0)
	part0 = particle
	old_neighbours = oldies
	new_neighbours = []
	for part1 in range(len(part_mat)):
		if part1 != part0:
			x0, y0, z0 = part_mat[part0]
			x1, y1, z1 = part_mat[part1]
			R = abs(distance_calc(x0, x1, y0, y1, z0, z1))
			if r_cut > R > 0:
				new_neighbours.append(part1)
				energy_mat[part0] = 0
				#energy_mat[part1] = 0
				E = lennard_jones(1, 1, R)
				energy_mat[part0] += E
				#energy_mat[part1] += E
	neighbours = list(set(old_neighbours + new_neighbours))
	return neighbours

# Update neighbour energies
def neighbour_energy(neighbours):
	for nei in neighbours:
		part0 = nei
		energy_mat[part0] = 0
		for part1 in range(part0 + 1, len(part_mat)):
			x0, y0, z0 = part_mat[part0]
			x1, y1, z1 = part_mat[part1]
			R = abs(distance_calc(x0, x1, y0, y1, z0, z1))
			if r_cut > R > 0:
				E = lennard_jones(1, 1, R)
				energy_mat[part0] += E
				energy_mat[part1] += E

# Metropolis-Hastings 
def metropolis(U, dU):
	delta = dU - U 
	P_accept = min(1, np.exp(-beta * delta))
	#print("P_accept: " + str(P_accept))
	P_random = r.randint(0, 1)
	if P_random <= P_accept:
		return True
	else:
		return False

# Monte Carlo loop
t = 0
T = 50000

# Energy history
energy_his = np.zeros(T)
#for i in range(num_particles):
get_energy_init()
energy_his[0] = np.mean(energy_mat)
U = energy_his[0]
acceptance_rate = 0

while t < T:
	if t % 100 == 0:
		print(t)
		#print(energy_his[t])
	
#	# Implement random walks
#	particle = r.randint(0, (num_particles - 1))
#	random_walk(particle, step_mag)
#		
	for part in range(len(part_mat) - 1):
		PDF(part)
#	
#	# Implement energy function
#	get_energy()
#	energy_his[t] = np.mean(energy_mat)
#	print(r"$\delta U$" + str(energy_his[t]))
#	test = (3/2) * 1.5 * (1.38 * 10 ** (-23)) + energy_his[t] * num_particles
#	print(test)
	
	particle = r.randint(0, (num_particles - 1))
	old_positions = copy.deepcopy(part_mat)
	oldies = old_neighbours(particle)
	random_walk(particle, step_mag)
	old_energy = copy.deepcopy(energy_mat)
	neighbours = get_energy(particle, oldies)
	neighbour_energy(neighbours)
	dU = np.mean(energy_mat)
	check = metropolis(U, dU)
	if check == True:
		energy_his[t] = np.mean(energy_mat) + (((3 / 2) * temp * 1.38 * 10 ** (-23)) )#* num_particles)
		U = dU
		#print(energy_his[t])
		acceptance_rate += 1
	else:
		energy_mat = old_energy
		part_mat = old_positions
		energy_his[t] = np.mean(energy_mat) + (((3 / 2) * temp * 1.38 * 10 ** (-23)) )#* num_particles) 
		#print("not accepted")
	
		
	
#		energy_mat.fill(0)
#		for part0 in range(len(part_mat)):
#			for part1 in range(part0 + 1, len(part_mat)):
#				x0, y0, z0 = part_mat[part0]
#				x1, y1, z1 = part_mat[part1]
#				R = abs(distance_calc(x0, x1, y0, y1, z0, z1))
#				if r_cut > R > 0:
#					E = lennard_jones(1, 1, R)
#					energy_mat[part0] += E
#					energy_mat[part1] += E
	#print(energy_mat)
	#print(part_mat)
	
	# Reset the energy matrix
#	energy_mat = np.zeros(num_particles)
	
#	# Plot the data
#	## 3D scatter plots
#	fig = plt.figure(facecolor = "black")
#	ax = fig.add_subplot(111, projection = "3d", facecolor = "black")
#	ax.set_xlim([0, box_length])
#	ax.set_ylim([0, box_length])
#	ax.set_zlim([0, box_length])
#	ax.grid(False)
#	ax.axis("off")
#	ax.scatter(part_mat[:, 0], part_mat[:, 1], part_mat[:, 2], 
#			color = "red", s = 10, alpha = 0.6)
	#plt.show()
	#plt.close()
	
#	plt.savefig("random_walk_images/" + str(t) + ".png")
#	plt.close()
	
#	## 2D scatter plots (x, y)
#	fig = plt.figure(facecolor = "black")
#	ax = fig.add_subplot(facecolor = "black")
#	ax.set_xlim([0, box_length])
#	ax.set_ylim([0, box_length])
#	ax.grid(False)
#	ax.axis("off")
#	ax.scatter(part_mat[:, 0], part_mat[:, 1], color = "red", s = 10, 
#			alpha = 0.6)
#	plt.savefig("random_walk_images/" + "xy_" + str(t) + ".png")
#	plt.close()
#	
#	## 2D scatter plots (x, z)
#	fig = plt.figure(facecolor = "black")
#	ax = fig.add_subplot(facecolor = "black")
#	ax.set_xlim([0, box_length])
#	ax.set_ylim([0, box_length])
#	ax.grid(False)
#	ax.axis("off")
#	ax.scatter(part_mat[:, 0], part_mat[:, 2], color = "red", s = 10, 
#			alpha = 0.6)
#	plt.savefig("random_walk_images/" + "xz_" + str(t) + ".png")
#	plt.close()
#	
#	## 2D scatter plots (y, z)
#	fig = plt.figure(facecolor = "black")
#	ax = fig.add_subplot(facecolor = "black")
#	ax.set_xlim([0, box_length])
#	ax.set_ylim([0, box_length])
#	ax.grid(False)
#	ax.axis("off")
#	ax.scatter(part_mat[:, 1], part_mat[:, 2], color = "red", s = 10, 
#			alpha = 0.6)
#	plt.savefig("random_walk_images/" + "yz_" + str(t) + ".png")
#	plt.close()
	t += 1

#	# 2D scatter plots (x, y)
#	fig = plt.figure(figsize = (15, 15))#facecolor = "black")
#	ax = fig.add_subplot()#facecolor = "black")
#	ax.set_xlim([0, box_length])
#	ax.set_ylim([0, box_length])
#	#ax.grid(False)
#	#ax.axis("off")
#	ax.set_xticks(np.arange(0, box_length, 1))
#	ax.set_yticks(np.arange(0, box_length, 1))
#	ax.grid(color="black", linestyle="-", linewidth=1)
#	ax.scatter(part_mat[:, 0], part_mat[:, 1], color = "red", s = 100, 
#			alpha = 0.6)
#	plt.savefig("random_walk_images/" + "xy_" + str(t) + str(IT) + ".png")
#	plt.close()

#	## 2D scatter plots (x, z)
#	fig = plt.figure(figsize = (15, 15))#facecolor = "black")
#	ax = fig.add_subplot()#facecolor = "black")
#	ax.set_xlim([0, box_length])
#	ax.set_ylim([0, box_length])
#	#ax.grid(False)
#	#ax.axis("off")
#	ax.set_xticks(np.arange(0, box_length, 1))
#	ax.set_yticks(np.arange(0, box_length, 1))
#	ax.grid(color="black", linestyle="-", linewidth=1)
#	ax.scatter(part_mat[:, 0], part_mat[:, 2], color = "red", s = 100, 
#			alpha = 0.6)
#	plt.savefig("random_walk_images/" + "xz_" + str(t)  +str(IT) + ".png")
#	plt.close()

#	## 2D scatter plots (y, z)
#	fig = plt.figure(figsize = (15, 15))#facecolor = "black")
#	ax = fig.add_subplot()#facecolor = "black")
#	ax.set_xlim([0, box_length])
#	ax.set_ylim([0, box_length])
#	#ax.grid(False)
#	#ax.axis("off")
#	ax.set_xticks(np.arange(0, box_length, 1))
#	ax.set_yticks(np.arange(0, box_length, 1))
#	ax.grid(color="black", linestyle="-", linewidth=1)
#	ax.scatter(part_mat[:, 1], part_mat[:, 2], color = "red", s = 100, 
#			alpha = 0.6)
#	plt.savefig("random_walk_images/" + "yz_" + str(t) + str(IT) + ".png")
#	plt.close()
print("Acceptance rate: " + str((acceptance_rate / T) * 100) + "%")
# Plot probability density (x-axis)
#pos_1dx_PDF_norm = pos_1dx_PDF / (T * num_particles)
#plt.bar(np.arange(box_length + 1), pos_1dx_PDF_norm, color = "black")
#plt.xlabel("x-axis position")
#plt.ylabel("probability density")
#plt.savefig("prob_dens_images/pdf_x" + ".png")
#plt.close()

## Plot probability density (y-axis)
#pos_1dy_PDF_norm = pos_1dy_PDF / (T * num_particles)
#plt.bar(np.arange(box_length + 1), pos_1dy_PDF_norm, color = "black")
#plt.xlabel("y-axis position")
#plt.ylabel("probability density")
#plt.savefig("prob_dens_images/pdf_y" + ".png")
#plt.close()

## Plot probability density (z-axis)
#pos_1dz_PDF_norm = pos_1dz_PDF / (T * num_particles)
#plt.bar(np.arange(box_length + 1), pos_1dz_PDF_norm, color = "black")
#plt.xlabel("z-axis position")
#plt.ylabel("probability density")
#plt.savefig("prob_dens_images/pdf_z" + ".png")
#plt.close()

## Save probability density (x-axis)
#np.save("xax", pos_1dx_PDF_norm)
#np.save("yax", pos_1dy_PDF_norm)
#np.save("zax", pos_1dz_PDF_norm)

# Plot the energy over time
np.save("energy_his_" + ID, energy_his)
plt.plot(np.arange(T), energy_his, color = "black", linewidth = 1)
plt.xlabel("iterations")
plt.ylabel(r"$\langle e \rangle$")
plt.savefig("energy_history_2" + ".png")
plt.close()
new_list = []
for _ in range(len(energy_his)):
	if _ < len(enegy_his) - 5:
		new_list.append(np.mean(energy_his[_] + energy_his[_+1] + energy_his[_+2] +
				energy_his[_+3] + energy_his[_+4]))

plt.plot(new_list)
plt.xlabel("iterations")
plt.ylabel(r"$\langle e \rangle$")
plt.show()
	
	


