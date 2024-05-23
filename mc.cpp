#include <iostream>
#include <cmath>
#include <fstream>
#include <random>

using namespace std;

// Prototipos de las funciones

/* Minimum image convention */
void mic(float * vec, float Lbox){

  for (int i = 0; i < 3 ; ++i) {
    vec[i] -= floorf(0.5 + vec[i]/Lbox)*Lbox;
  }
  
  return;
}

float Energy(float Position[], float N_O_Pos[], int itag);

// Global variable

float sigma = 1.0;
float epsilon = 1.0;
float Lbox = 15.0;
float densinput = 0.1;
float volume = Lbox*Lbox*Lbox;
int Npart=(int)(densinput*volume);
int dim = 3;
float Temp=1.5;
float r_cut = 2.5 * sigma;

int main (void)
{
    random_device rand_dev;
    mt19937 generator(rand_dev());

    float dens = Npart/volume;
    // LOG of the run
    cout <<"Density "<< dens <<endl;
    cout <<"Temperature "<< Temp <<endl;
    // system variables

    float Position[dim*Npart];

    // Montecarlo variable

    int Nstep = 50000;
    int nsamp_ener = 10;
    int nsamp_pos = 100;
    int naccept = 0;
    float deltaR = 0.1;

    // Data files

    ofstream fich_ener, fich_posi;

    fich_ener.open("energy.txt");
    fich_posi.open("position.txt");


    uniform_real_distribution<> dis1(-1.0, 1.0);
    uniform_real_distribution<> dis2(0.0, 1.0);
    uniform_int_distribution<> dist(0, Npart-1);

    // Initial positions

    for (int i=0; i<3*Npart; ++i) Position[i] = dis1(generator)*Lbox/2.;

    // Monte carlo loop

    for (int istep = 0; istep<Nstep; ++istep)
    {

        int itag = dist(generator);

	
        float PosNew[dim], PosOld[dim];
        float Enew, Eold, prob;

        // Move particle itag

        for (int k=0; k<dim; ++k)
        {
            PosOld[k] = Position[dim*itag + k];
            PosNew[k] = Position[dim*itag + k] + deltaR*dis1(generator);
        } 

        Eold = Energy(Position, PosOld, itag);
        Enew = Energy(Position, PosNew, itag);

	/*	cout<<Eold<<" "<<Enew<<" "<<deltaR<<endl;*/
        float Esample = Eold;

        // Probabiliy ratio

        prob = exp(-(Enew-Eold)/Temp);

        float xi = dis2(generator);

        if (prob > xi) // If accept
        {
            for (int k=0; k<dim; ++k) Position[itag*dim + k] = PosNew[k];

            naccept = naccept + 1;

            Esample = Enew;
        }

       if (istep % nsamp_ener == 0) fich_ener << 0.5*Esample << endl;
       /* if (istep % nsamp_pos == 0) 
	        {
        
		}*/
    }


    fich_ener.close();
    fich_posi.close();

    return 0.;
}


// Funcion que calcula la energia

float Energy(float Position[], float Pos_itag[], int itag)
{
  float U_tot = 0.0;
  float r_2cut = r_cut * r_cut;

  for ( int i = 0; i < Npart; ++i)  {
    if (i == itag) continue;
    float vec_dist[dim];
    float u_Ij = 0.0;
    
    for (int k = 0; k < dim; ++k)
      vec_dist[k] = (Pos_itag[k] - Position[i * dim + k]);
    
    mic(vec_dist, Lbox);

    /*    cout << vec_dist[0]<<" "<<vec_dist[1]<<" "<<vec_dist[2]<<endl;*/
      
    float r2_Ij=0.0;
    for (int k = 0; k < dim; ++k)
      { r2_Ij += pow(vec_dist[k], 2);   
	/*	cout << r2_Ij<<endl;*/
      }
    //lennard-jones potential
	
    if (r2_Ij < r_2cut)
      {
	float r_mod = sqrt(r2_Ij);
	float r12 = pow((sigma/r_mod), 12.0);
	float r6 = pow((sigma/r_mod), 6.0);
	float rc12 = pow((sigma/r_cut), 12.0);
	float rc6 = pow((sigma/r_cut), 6.0);
	
	u_Ij = 4 * epsilon * ((r12-r6) - (rc12-rc6));   
      }
    U_tot += u_Ij;
      }
  return U_tot;  
}

