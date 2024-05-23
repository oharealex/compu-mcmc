# 3-D random walk
# Created by Alex O'Hare 29/04/2024

# Import the required libraries
import matplotlib.pyplot as plt
import numpy as np
import random as r
import copy

simulation_number = 0
density_values = [0.1, 0.17, 0.24, 0.31, 0.38, 0.45, 0.52, 0.59, 0.66, 0.73, 0.8]
temperature_values = [1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4]
history_set = np.array()
for dv in density_values:
	for tv in temperature_values:
	
		simulation_number += 1
		print("simulation number #" + str(simulation_number))
		
		# Save the data
		# temperature, pressure, number of particles, mean energy across iterations,
		# min energy across iterations, max energy across iterations, error.
		history = [tv, dv] 
		# System settings
		ID = "001"
		kB = 1.380649 * (10 ** (-23))
		temp = tv
		beta = 1 / (kB * temp)
		box_length = 12
		box_volume = box_length * box_length * box_length
		density = dv
		num_particles = int(density * box_volume)
		history.append(num_particles)
		step_mag = 2
		r_cut = 2.5
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
					if dis <= 0: # r_cut:
						pass
					else:
						part_mat[particle][dir_axis] = new_position

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
		get_energy_init()
		energy_his[0] = np.mean(energy_mat)
		U = energy_his[0]
		acceptance_rate = 0

		while t < T:
			
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
				energy_his[t] = np.mean(energy_mat) + (((3 / 2) * temp * 1.38 * 10 ** (-23)) )
				U = dU
				acceptance_rate += 1
			else:
				energy_mat = old_energy
				part_mat = old_positions
				energy_his[t] = np.mean(energy_mat) + (((3 / 2) * temp * 1.38 * 10 ** (-23)) )

		# Save the data
		np.save("energy_his_" + ID + str(nv) + str(tv), energy_his)
		history.append(np.mean(energy_his))
		history.append(np.min(energy_his))
		history.append(np.max(energy_his))
		history.append((np.max(energy_his) - np.min(energy_his)) / 2)
		history_set.append(history)
		



