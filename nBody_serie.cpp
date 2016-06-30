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
	
	//lectura de masas
	n=36;
	int nsteps=10000;
	
	
	double G=6.673*1e-11;
	for(int i=0;i<n;i++){
		masses[i]= (rand()/10.0);
	}
	
	for(int i=0;i<n;i++){
		pos[i][0]=rand()%100-50;
		pos[i][1]=rand()%100-50;
		vel[i][0]=rand()%100-50;
		vel[i][1]=rand()%100-50;
	}
	
	
	double x_diff,y_diff,dist,dist_cubed;
	double delta_t=1e-5;
	
	for(int timestep=0;timestep<nsteps;timestep++){
		if(timestep==nsteps-1){
			for(int i=0;i<5;i++)
				cout<<pos[i][0]<<" "<<pos[i][1]<<" "<<vel[i][0]<<" "<<vel[i][1]<<endl;
		}
		
		
		for(int i=0;i<n;i++)
			forces[i][0]=0,forces[i][1]=0;
		
		for(int q=0;q<n;q++)
			for(int k=q+1;k<n;k++){
				x_diff= pos[q][0]-pos[k][0];
				y_diff= pos[q][1]-pos[k][1];
				dist= sqrt(abs(x_diff*x_diff+ y_diff*y_diff));
				dist_cubed=dist*dist*dist;
				forces_qk[0]= G*masses[q]*masses[k]/dist_cubed*x_diff;
				forces_qk[1]= G*masses[q]*masses[k]/dist_cubed*y_diff;
				
				forces[q][0]+=forces_qk[0];
				forces[q][1]+=forces_qk[1];
				forces[k][0]-=forces_qk[0];
				forces[k][1]-=forces_qk[1];
			}
		
		for(int q=0;q<n;q++){
			pos[q][0]+=delta_t*vel[q][0];
			pos[q][1]+=delta_t*vel[q][1];
			vel[q][0]+=delta_t/masses[q]*forces[q][0];
			vel[q][1]+=delta_t/masses[q]*forces[q][1];
		}
	}
	
	
    return 0;
}

