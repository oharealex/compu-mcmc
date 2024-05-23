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
