// Monte Carlo for Lennard-Jones fluids
// Created by Alex O'Hare, 28/04/2024

// Import the required libraries
#include <iostream>
#include <cmath>
#include <fstream>
#include <random>

using namespace std;

// Minimum image convention
void minimum_image_convention(float * vec, float box_length){
	for (int i = 0; i < 3; ++i){
		vec[i] -= floorf(0.5 + vec[i] / box_length) * box_length;
		}
	return;
}

float energy(float position[], float N_O_Pos[], int itag);

// Define global variables
float sigma = 1.0;
float epsilon = 1.0;
float box_length = 12.0;
float density_init = 0.1;
float volume = pow(box_length, 3);
int num_particles = (int)(density_init * volume);
int dim = 3;
float temperature = 1.5;
float r_cut = 2.5 * sigma;

// Main loop
int main (void){
	
	// Initialise and seed the random number generator
	random_device rand_dev;
	mt19937 generator(rand_dev());
	
	// Calculate the density
	float density = num_particles / volume;
	
	// Print the data
	cout << "Density" << density << endl;
	cout << "Temperature" << temperature << endl;
	
	// System variables
	float position[dim * num_particles];
	
	// Monte Carlo variables
	int num_steps = 50000;
	int energy_step = 10;
	int pos_step = 100;
	int naccept = 0;
	float max_disp = 0.01;
	
	// Data records
	ofstream save_energy, save_position;
	save_energy.open("energy.txt");
	save_position.open("position.txt");
	
	// Define random distributions
	uniform_real_distribution<> dis1(-1.0, 1.0);
	uniform_real_distribution<> dis2(0.0, 1.0);
	uniform_int_distribution<> dis3(0, num_particles - 1);
	
	// Generate initial positions
	for (int i = 0; i < 3 * num_particles; ++i){
		position[i] = dis1(generator) * box_length / 2.0;
		}
	
	// Monte Carlo loop
	for (int istep = 0; istep < num_steps; ++istep){
	
		// Select a particle at random
		int itag = dis3(generator);
		
		//Define variables for new pos., current pos., 
		float new_pos[dim], old_pos[dim];
		
		// Define variables for new energy, current energy, prob. of acceptance
		float new_ene, old_ene, prob;
		
		// Move particle
		for (int k = 0; k < dim; ++k){
			old_pos[k] = position[dim * itag + k];
			new_pos[k] = position[dim * itag + k] + max_disp * dis1(generator);
			}
		old_ene = energy(position, old_pos, itag);
		new_ene = energy(position, new_pos, itag);
		float sample_ene = old_ene;
		
		// Probability ratio for acceptance
		prob = exp(-(new_ene - old_ene) / temperature);
		float xi = dis2(generator);
		if (prob > xi){
			for (int k = 0; k < dim; ++k) position[itag * dim + k] = new_pos[k];
			naccept = naccept + 1;
			sample_ene = new_ene;
			}
		
		if (istep % energy_step == 0) save_energy << 0.5 * sample_ene << endl;
		}
	
	save_energy.close();
	save_position.close();
	// Print the data
	cout << "Density" << density << endl;
	cout << "Temperature" << temperature << endl;

	return 0;
}

// Energy calculation function
float energy(float position[], float pos_itag[], int itag){
	float phi_total = 0.0;
	float r_cut_sq = r_cut * r_cut;
	
	for (int i = 0; i < num_particles; ++i){
		if (i == itag) continue;
		float vec_dist[dim];
		float phi_ij = 0.0;
		
		for (int k = 0; k < dim; ++k){
			vec_dist[k] = (pos_itag[k] - position[i * dim + k]);
			}
			
		minimum_image_convention(vec_dist, box_length);
		
		float r_ij_sq = 0.0;
		for (int k = 0; k < dim; ++k){
			r_ij_sq += pow(vec_dist[k], 2);
			}
		
		// Lennard-Jones potential
		if (r_ij_sq < r_cut_sq){
			float r_mod = sqrt(r_ij_sq);
			float r12 = pow((sigma/r_mod), 12.0);
			float r6 = pow((sigma/r_mod), 6.0);
			float rc12 = pow((sigma/r_cut), 12.0);
			float rc6 = pow((sigma/r_cut), 6.0);
			phi_ij = 4 * epsilon * ((r12 - r6) - (rc12 - rc6));
			}
		phi_total += phi_ij;
		}
	return phi_total;
	}

		
