#include<iostream>
#include<cstring>
#include<cstdio>
#include<vector>
#include<cmath>
#include<sstream>
#include<algorithm>
#include<set>
#include<ctime>
#include<omp.h>
using namespace std;

double pos[101][2];
double vel[101][3];
double forces[101][2];
double forces_qk[2];
double masses[101];
int n;

int main(){
	omp_set_num_threads(4);  
		//lectura de masas
	n=36;
	int nsteps=10000;
	
	double G=6.673*1e-11;
	for(int i=0;i<n;i++){
		masses[i]= ( rand()/10.0);
	}
	
	for(int i=0;i<n;i++){
		pos[i][0]=rand()%100-50;
		pos[i][1]=rand()%100-50;
		vel[i][0]=rand()%100-50;
		vel[i][1]=rand()%100-50;
	}
	
	
	//double x_diff,y_diff,dist,dist_cubed;
	double delta_t=1e-5;
	#pragma omp parallel
	for(int timestep=0;timestep<nsteps;timestep++){
		//cout<<timestep<<" "<<omp_get_thread_num()<<endl;		
		if(timestep==nsteps-1){
			#pragma omp single
			for(int i=0;i<5;i++)
				cout<<pos[i][0]<<" "<<pos[i][1]<<" "<<vel[i][0]<<" "<<vel[i][1]<<endl;
		}
		
		#pragma omp for schedule(static, n/4)
		for(int q=0;q<n;q++){
			forces[q][0]=0,forces[q][1]=0;
			for(int k=0;k<n;k++){
				if(k==q)continue;
				double x_diff= pos[q][0]-pos[k][0];
				double y_diff= pos[q][1]-pos[k][1];
				double dist= sqrt( abs(x_diff*x_diff+ y_diff*y_diff));
				double dist_cubed=dist*dist*dist;
				forces[q][0]+=G*masses[q]*masses[k]/dist_cubed*x_diff;
				forces[q][1]+=G*masses[q]*masses[k]/dist_cubed*y_diff;
			}	
		}
		
		#pragma omp for schedule(static, n/4)
		for(int q=0;q<n;q++){
			pos[q][0]+=delta_t*vel[q][0];
			pos[q][1]+=delta_t*vel[q][1];
			vel[q][0]+=(delta_t/masses[q])*forces[q][0];
			vel[q][1]+=(delta_t/masses[q])*forces[q][1];
		}
		
	}

	
    return 0;
}

