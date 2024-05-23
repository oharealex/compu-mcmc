import matplotlib.pyplot as plt
import numpy as np

mat0x = np.load("0xax.npy")
mat1x = np.load("1xax.npy")
mat2x = np.load("2xax.npy")
mat3x = np.load("3xax.npy")
mat4x = np.load("4xax.npy")
mat5x = np.load("5xax.npy")
mat6x = np.load("6xax.npy")
mat7x = np.load("7xax.npy")
mat8x = np.load("8xax.npy")
mat9x = np.load("9xax.npy")

mat0y = np.load("0yax.npy")
mat1y = np.load("1yax.npy")
mat2y = np.load("2yax.npy")
mat3y = np.load("3yax.npy")
mat4y = np.load("4yax.npy")
mat5y = np.load("5yax.npy")
mat6y = np.load("6yax.npy")
mat7y = np.load("7yax.npy")
mat8y = np.load("8yax.npy")
mat9y = np.load("9yax.npy")

mat0z = np.load("0zax.npy")
mat1z = np.load("1zax.npy")
mat2z = np.load("2zax.npy")
mat3z = np.load("3zax.npy")
mat4z = np.load("4zax.npy")
mat5z = np.load("5zax.npy")
mat6z = np.load("6zax.npy")
mat7z = np.load("7zax.npy")
mat8z = np.load("8zax.npy")
mat9z = np.load("9zax.npy")

collectionx = np.stack((mat0x, mat1x, mat2x, mat3x, mat4x, mat5x, mat6x, mat7x,
		mat8x, mat9x), axis = 0)

collectiony = np.stack((mat0y, mat1y, mat2y, mat3y, mat4y, mat5y, mat6y, mat7y,
		mat8y, mat9y), axis = 0)

collectionz = np.stack((mat0z, mat1z, mat2z, mat3z, mat4z, mat5z, mat6z, mat7z,
		mat8z, mat9z), axis = 0)

mean_valsx = np.mean(collectionx, axis = 0)
mean_arrayx = np.array(mean_valsx)

mean_valsy = np.mean(collectiony, axis = 0)
mean_arrayy = np.array(mean_valsy)

mean_valsz = np.mean(collectionz, axis = 0)
mean_arrayz = np.array(mean_valsz)

plt.bar(np.arange(13), mean_arrayx, color = "black")
plt.xlabel("x-axis position")
plt.ylabel("mean probability density")
plt.savefig("meanx.png")
plt.close()
plt.bar(np.arange(13), mean_arrayy, color = "black")
plt.xlabel("y-axis position")
plt.ylabel("mean probability density")
plt.savefig("meany.png")
plt.close()
plt.bar(np.arange(13), mean_arrayz, color = "black")
plt.xlabel("z-axis position")
plt.ylabel("mean probability density")
plt.savefig("meanz.png")
plt.close()
