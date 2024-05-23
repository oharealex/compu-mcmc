# 3-D random walk
# Created by Alex O'Hare 29/04/2024

# Import the required libraries
import matplotlib.pyplot as plt
import numpy as np
import random as r


for IT in range(1):
	# System settings
	box_length = 12
	box_volume = box_length * box_length * box_length
	density = 0.1
	num_particles = int(density * box_volume)
	step_mag = 3
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
		dir_mag = r.randint(1, step_mag)
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
				if dis < r_cut:
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
	
	# Energy function
	def get_energy():
		energy_mat.fill(0)
		for part0 in range(len(part_mat)):
			for part1 in range(part0 + 1, len(part_mat)):
				x0, y0, z0 = part_mat[part0]
				x1, y1, z1 = part_mat[part1]
				R = abs(distance_calc(x0, x1, y0, y1, z0, z1))
				if r_cut > R > 0:
					E = lennard_jones(1, 1, R)
					energy_mat[part0] += E
					energy_mat[part1] += E
		
	# Monte Carlo loop
	t = 0
	T = 10000
	
	# Energy history
	energy_his = np.zeros(T)
	
	while t < T:
		print(t)
		
		# Implement random walks
		for particle in range(num_particles):
			random_walk(particle, step_mag)
			
		for part in range(len(part_mat) - 1):
			PDF(part)
		
		# Implement energy function
		get_energy()
		energy_his[t] = np.mean(energy_mat)
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

	# Plot probability density (x-axis)
	pos_1dx_PDF_norm = pos_1dx_PDF / (T * num_particles)
	plt.bar(np.arange(box_length + 1), pos_1dx_PDF_norm, color = "black")
	plt.xlabel("x-axis position")
	plt.ylabel("probability density")
	plt.savefig("prob_dens_images/pdf_x" + str(IT) + ".png")
	plt.close()

	# Plot probability density (y-axis)
	pos_1dy_PDF_norm = pos_1dy_PDF / (T * num_particles)
	plt.bar(np.arange(box_length + 1), pos_1dy_PDF_norm, color = "black")
	plt.xlabel("y-axis position")
	plt.ylabel("probability density")
	plt.savefig("prob_dens_images/pdf_y" + str(IT) + ".png")
	plt.close()

	# Plot probability density (z-axis)
	pos_1dz_PDF_norm = pos_1dz_PDF / (T * num_particles)
	plt.bar(np.arange(box_length + 1), pos_1dz_PDF_norm, color = "black")
	plt.xlabel("z-axis position")
	plt.ylabel("probability density")
	plt.savefig("prob_dens_images/pdf_z" + str(IT) + ".png")
	plt.close()

	# Save probability density (x-axis)
	np.save(str(IT)+"xax", pos_1dx_PDF_norm)
	np.save(str(IT)+"yax", pos_1dy_PDF_norm)
	np.save(str(IT)+"zax", pos_1dz_PDF_norm)
	
	# Plot the energy over time
	plt.plot(np.arange(T), energy_his, color = "black", linewidth = 1)
	plt.xlabel("iterations")
	plt.ylabel(r"$<\phi(r)>$")
	plt.savefig("energy_history_" + str(IT) + ".png")
	


