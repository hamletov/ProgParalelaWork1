#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <iostream>

typedef double vect_t[2];  
const double G = 6.673e-11;  
int my_rank, comm_sz;
MPI_Comm comm;
MPI_Datatype vect_mpi_t;
vect_t* vel;

void get_arg(int* n_p, int* n_steps_p, double* delta_t_p);
void Gen_init_cond(double masses[], vect_t pos[],  vect_t loc_vel[], int n, int loc_n);
void Salida(double time, double masses[], vect_t pos[], vect_t loc_vel[], int n, int loc_n);
void Obtener_Fuerza(int loc_part, double masses[], vect_t loc_forces[], vect_t pos[], int n, int loc_n);
void Actualiza_PosVel(int loc_part, double masses[], vect_t loc_forces[], vect_t loc_pos[], vect_t loc_vel[], int n, int loc_n, double delta_t);

using namespace std;
int main(int argc, char* argv[]) {
	int n;int loc_n; //n total de particulas, loc_n cantidad de particulas asignada por proceso  
	int n_steps;       
	int loc_part;           
	double delta_t;                   
	double* masses;             
	vect_t* loc_pos;vect_t* pos;                
	vect_t* loc_vel; vect_t* loc_forces;         

	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &comm_sz);
	MPI_Comm_rank(comm, &my_rank);

	get_arg(&n, &n_steps, &delta_t);
	loc_n = n/comm_sz; 
	masses = (double *)malloc(n*sizeof(double));pos = (vect_t *)malloc(n*sizeof(vect_t));loc_forces = (vect_t *)malloc(loc_n*sizeof(vect_t));
   
	loc_pos = pos + my_rank*loc_n;
	loc_vel = (vect_t *)malloc(loc_n*sizeof(vect_t));
	if (my_rank == 0)vel=(vect_t *)malloc(n*sizeof(vect_t));
	MPI_Type_contiguous(2, MPI_DOUBLE, &vect_mpi_t);
	MPI_Type_commit(&vect_mpi_t);

    Gen_init_cond(masses, pos, loc_vel, n, loc_n);
	/*
	if(my_rank==0)
		for(int i=0;i<n;i++)
			cout<<pos[i][0]<<" "<<pos[i][1]<<" "<<vel[1][0]<<" "<<vel[1][1]<<endl;
	*/

    double start = MPI_Wtime();
	for(int step = 1; step <= n_steps; step++) {
		double t = step*delta_t;
		if(step==n_steps)Salida(t,masses,pos,loc_vel,n,loc_n);
		for(loc_part = 0; loc_part < loc_n; loc_part++)
			Obtener_Fuerza(loc_part, masses, loc_forces, pos, n, loc_n);
		for(loc_part = 0; loc_part < loc_n; loc_part++)
			Actualiza_PosVel(loc_part, masses, loc_forces, loc_pos, loc_vel, n, loc_n, delta_t);
		MPI_Allgather(MPI_IN_PLACE, loc_n, vect_mpi_t,pos, loc_n, vect_mpi_t, comm);
	}
   
	double finish = MPI_Wtime();
	if (my_rank == 0)printf("Elapsed time = %e seconds\n", finish-start);

	MPI_Type_free(&vect_mpi_t);
	
	free(masses);free(pos);free(loc_forces);free(loc_vel);
	if(my_rank==0)free(vel);
	MPI_Finalize();
	return 0;
}  

/*---------------------------------------------------------------------
 * Funcion:  Get_args
 * Proposito:   Obtener datos de entrada
 * Salida:
 *    n_p:             Numero de particulas
 *    n_steps_p:       Numero de timesteps
 *    delta_t_p:       tamanio de cada paso
 */
void get_arg(int* n_p, int* n_steps_p, double* delta_t_p) {
   if (my_rank == 0) {
	   *n_p=12;
	   *n_steps_p =1000;
	   *delta_t_p=1e-3;
   }
   MPI_Bcast(n_p, 1, MPI_INT, 0, comm);
   MPI_Bcast(n_steps_p, 1, MPI_INT, 0, comm);
   MPI_Bcast(delta_t_p, 1, MPI_DOUBLE, 0, comm);
}

void Gen_init_cond(double masses[], vect_t pos[], 
      vect_t loc_vel[], int n, int loc_n) {
	
	if (my_rank == 0) {
		for(int i=0;i<n;i++)
			masses[i]= (rand()/10.0);
	
		for(int i=0;i<n;i++){
			pos[i][0]=rand()%100-50;
			pos[i][1]=rand()%100-50;
			vel[i][0]=rand()%100-50;
			vel[i][1]=rand()%100-50;
		}	  
   }

	MPI_Bcast(masses, n, MPI_DOUBLE, 0, comm);
	MPI_Bcast(pos, n, vect_mpi_t, 0, comm);
	MPI_Scatter(vel, loc_n, vect_mpi_t,loc_vel, loc_n, vect_mpi_t, 0, comm);
}  

void Salida(double time, double masses[], vect_t pos[],
      vect_t loc_vel[], int n, int loc_n) {
	
	MPI_Gather(loc_vel, loc_n, vect_mpi_t, vel, loc_n, vect_mpi_t, 0, comm);
	if (my_rank == 0) {
		printf("%.2f\n", time);
		for(int part = 0; part < n; part++) {
			printf("%3d %10.3e ", part, pos[part][0]);
			printf("  %10.3e ", pos[part][1]);
			printf("  %10.3e ", vel[part][0]);
			printf("  %10.3e\n", vel[part][1]);
		}
		printf("\n");
	}
}  

void Obtener_Fuerza(int loc_part, double masses[], vect_t loc_forces[], 
      vect_t pos[], int n, int loc_n) {
	vect_t f_part_k;

	int part = my_rank*loc_n + loc_part;
	loc_forces[loc_part][0] = loc_forces[loc_part][1] = 0.0;
	for(int k = 0; k < n; k++) {
		if (k != part) {
			f_part_k[0] = pos[part][0] - pos[k][0];
			f_part_k[1] = pos[part][1] - pos[k][1];
			double len = sqrt(f_part_k[0]*f_part_k[0] + f_part_k[1]*f_part_k[1]);
			double len_3 = len*len*len;
			double mg = -G*masses[part]*masses[k];
			double fact = mg/len_3;
			f_part_k[0] *= fact;
			f_part_k[1] *= fact;
   
			loc_forces[loc_part][0] += f_part_k[0];
			loc_forces[loc_part][1] += f_part_k[1];
      }
   }
}

/*---------------------------------------------------------------------
 * Function:  Actualiza_PosVel
 * Purpose:   Actualiza la velocidad y posicion
 * In args:
 *    loc_part:    indice a actualizar
 *    masses:      masa
 *    loc_forces:  arreglo local de fuerzas
 *    n:           numero de particulas
 *    loc_n:       particulas asignadas para este proceso
 *    delta_t:     step size
 */
void Actualiza_PosVel(int loc_part, double masses[], vect_t loc_forces[], 
      vect_t loc_pos[], vect_t loc_vel[], int n, int loc_n, 
      double delta_t) {
		  
   int part = my_rank*loc_n + loc_part;
   double fact = delta_t/masses[part];
   loc_pos[loc_part][0]+= delta_t * loc_vel[loc_part][0];
   loc_pos[loc_part][1]+= delta_t * loc_vel[loc_part][1];
   loc_vel[loc_part][0]+= fact * loc_forces[loc_part][0];
   loc_vel[loc_part][1]+= fact * loc_forces[loc_part][1];
} 

