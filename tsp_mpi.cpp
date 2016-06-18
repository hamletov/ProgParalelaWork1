#include<iostream>
#include<cstring>
#include<cstdio>
#include<stack>
#include<queue>
#include<vector>
#include<cmath>
#include<sstream>
#include<algorithm>
#include<set>
#include<ctime>
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<utility>

int myrank, NumProcs, NumCities;
int *dis;
int init_tour_count;
int n=12;
double minimo_global=1e+10;

struct path{
	int id,distance,mask;
	path(){}
	path(int _id,int _distance,int _mask){
		id=_id;distance=_distance;mask=_mask;
	}
};
path partial[1001];
double ans[25];
int *shared_id;
int *shared_distance;
int *shared_mask;

void Fill_Dist(){
	if(myrank == 0)NumCities=12;
	// global operation, all processes must call it
	MPI_Bcast( &NumCities, 1, MPI_INT, 0, MPI_COMM_WORLD);
	dis = (int*)malloc(NumCities*NumCities*sizeof(int));

	if(myrank == 0){
		for(int i=0;i<NumCities;i++)dis[i*NumCities+i]=0;
		
		for(int i=0;i<NumCities;i++)
			for(int j=i+1;j<NumCities;j++){
				dis[i*NumCities+j]=rand()%100;
				dis[j*NumCities+i]=dis[i*NumCities+j];
			}
	}
  
	//Broadcast
	MPI_Bcast(dis,NumCities*NumCities,MPI_INT, 0, MPI_COMM_WORLD); 
  
	if (myrank == 0){
		printf("Number of cities: %d\n", NumCities);
		for(int i=0;i<NumCities;i++){
			for(int j=0;j<NumCities;j++)
				printf("%5d",dis[i*NumCities+j] );
			printf("\n");
		}
	}
}

int level(int t){
	int fact=1;
	for(int i=n-1;i>0;i--){
		fact*=i;
		if(fact>=t)return n-i;
	}
	return 1;
}

int f(int t){
	int fact=1;
	for(int i=n-1;i>0;i--){
		fact*=i;
		if(fact>=t)return fact;
	}
	return 1;
}

void Build_initial_queue() {
	int val=f(NumProcs);
	shared_id=(int*)malloc(val*sizeof(int));
	shared_distance=(int*)malloc(val*sizeof(int));
	shared_mask=(int*)malloc(val*sizeof(int));
	int curr_sz = 0;
	
	if (myrank==0){
		curr_sz = 0;
		std::queue<path>Q;
		Q.push(path(0,0,1));
		int niv= level(NumProcs);
		while(!Q.empty()){
			path p =Q.front();
			Q.pop();
			int id=p.id;
			double distance=p.distance;
			int mask=p.mask;
		
			if( __builtin_popcount(mask) ==niv+1){
				partial[curr_sz++]=  path(id,distance,mask );
			}else{
				for(int i=0;i<n;i++)
					if( (mask&(1<<i))==0 )
						Q.push(path(i,distance+dis[id*NumProcs+i],mask|(1<<i)  )  );
			}
		}
		
	
		for(int i=0;i<curr_sz;i++){
			shared_id[i]=partial[i].id;
			shared_distance[i]=partial[i].distance;
			shared_mask[i]=partial[i].mask;
		}
	}

	MPI_Bcast(&curr_sz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(shared_id,curr_sz,MPI_INT, 0, MPI_COMM_WORLD);   	
	MPI_Bcast(shared_distance,curr_sz,MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(shared_mask,curr_sz,MPI_INT, 0, MPI_COMM_WORLD);

	init_tour_count = curr_sz;
}

void Set_init_tours(int* my_first_tour_p, int* my_last_tour_p) {
   int quotient, remainder, my_count;
   quotient = init_tour_count/NumProcs;
   remainder = init_tour_count%NumProcs;
   if (myrank < remainder) {
      my_count = quotient+1;
      *my_first_tour_p = myrank*my_count;
   } else {
      my_count = quotient;
      *my_first_tour_p = myrank*my_count + remainder;
   }
   *my_last_tour_p = *my_first_tour_p + my_count - 1;
} 

void Partition_tree(std::queue<std::pair<std::pair<int,double>,int> > &Q) {
   	int my_first_tour=0, my_last_tour=0, i;
	Set_init_tours(&my_first_tour, &my_last_tour);
   	
	for(i = my_last_tour; i >= my_first_tour; i--)
		Q.push(std::make_pair(std::make_pair(shared_id[i],shared_distance[i]),shared_mask[i]));	
   
}  


void Par_tree_search() {
	std::queue<std::pair<std::pair<int,double>,int> >Q;
	Partition_tree(Q);
	double minimo=1e+10;
	
	while(!Q.empty()){
		std::pair<std::pair<int,double>,int>p =Q.front();
		Q.pop();
		int id=p.first.first;
		double distance=p.first.second;
		int mask=p.second;
		
		if( __builtin_popcount(mask) ==n){
			minimo=std::min(minimo,distance+dis[id*NumProcs+0]);			
		}else{
			for(int i=0;i<n;i++)
				if( (mask&(1<<i))==0 )
					Q.push(std::make_pair(std::make_pair(i,distance+dis[id*NumProcs+i]),mask|(1<<i)  )  );
		}
	}
	
	ans[myrank]=minimo;
	MPI_Reduce(&minimo, &minimo_global, 1, MPI_DOUBLE, MPI_MIN, 0,
           MPI_COMM_WORLD);
	//minimo_global
	
	return;
}


using namespace std;
int main(int argc, char* argv[]) {
	MPI_Init (&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
	for(int i=0;i<NumProcs;i++)ans[i]=1e+10;
	double start = MPI_Wtime();
	
	Fill_Dist();  // process 0 master
	Build_initial_queue();
	Par_tree_search();
	
	double dev=1e+10;
	
	double finish = MPI_Wtime();
	
	if(myrank==0){
		cout<<minimo_global<<endl;
		printf("Elapsed time = %e seconds\n", finish-start);
	}
	

	MPI_Finalize();
	return 0;
}
