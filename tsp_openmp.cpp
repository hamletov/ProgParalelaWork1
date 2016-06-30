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
#include<omp.h>
using namespace std;
int thread_count;
bool visited[31];
double dis[31][31];
int n=12;
int init_tour_count;
pair<pair<int,double>,int>partial[1001];
double ans[25];

int f(int t){
	int fact=1;
	for(int i=n-1;i>0;i--){
		fact*=i;
		if(fact>=t)return n-i;
	}
	return 1;
}

void Set_init_tours(long long my_rank, int* my_first_tour_p, int* my_last_tour_p) {
   int quotient, remainder, my_count;
   quotient = init_tour_count/thread_count;
   remainder = init_tour_count % thread_count;
   if (my_rank < remainder) {
      my_count = quotient+1;
      *my_first_tour_p = my_rank*my_count;
   } else {
      my_count = quotient;
      *my_first_tour_p = my_rank*my_count + remainder;
   }
   *my_last_tour_p = *my_first_tour_p + my_count - 1;
} 

void Build_initial_queue(){
   	int curr_sz = 0;
   	queue<pair<pair<int,double>,int> >Q;
	Q.push(make_pair(make_pair(0,0),1));
	int niv= f(thread_count);
	
   	while(!Q.empty()){
   		pair<pair<int,double>,int>p =Q.front();
		Q.pop();
		int id=p.first.first;
		double distance=p.first.second;
		int mask=p.second;
		
		if( __builtin_popcount(mask) ==niv+1){
			partial[curr_sz++]=make_pair(make_pair(id,distance),mask );
		
		}else{
			for(int i=0;i<n;i++)
				if( (mask&(1<<i))==0 )
					Q.push(make_pair(make_pair(i,distance+dis[id][i]),mask|(1<<i)  )  );
		}
   	}
   	init_tour_count = curr_sz; 
}

void Partition_tree(int my_rank, queue<pair<pair<int,double>,int> > &Q) {
   	int my_first_tour=0, my_last_tour=0, i;
	Set_init_tours(my_rank, &my_first_tour, &my_last_tour);

	for(i = my_last_tour; i >= my_first_tour; i--){
		Q.push(partial[i]);	
    }
}  

void Par_tree_search() {
	int my_rank = omp_get_thread_num();
   	queue<pair<pair<int,double>,int> >Q;
	Partition_tree(my_rank, Q);
	
	while(!Q.empty()){
		pair<pair<int,double>,int>p =Q.front();
		Q.pop();
		int id=p.first.first;
		double distance=p.first.second;
		int mask=p.second;
		
		if( __builtin_popcount(mask) ==n){
			ans[my_rank]=min(ans[my_rank],distance+dis[id][0]);
		}else{
			for(int i=0;i<n;i++)
				if( (mask&(1<<i))==0 )
					Q.push( make_pair(make_pair(i,distance+dis[id][i]),mask|(1<<i)  )  );
		}
	}
	
	return ;
}

int main(){
	
	for(int i=0;i<n;i++)dis[i][i]=0;
	
	for(int i=0;i<n;i++)
		for(int j=i+1;j<n;j++){
			dis[i][j]=rand()%100;
			dis[j][i]=dis[i][j];
		}
	thread_count=8;
	omp_set_num_threads(thread_count);
	for(int i=0;i<thread_count;i++)ans[i]=1e+10;
	Build_initial_queue();

	#pragma omp parallel for
	for(int thread = 0; thread < thread_count; thread++)
		Par_tree_search();
	
   	double mini=1e+10;
   	for(int i=0;i<thread_count;i++)mini=min(mini,ans[i]);
   	
   	printf("%.10lf\n",mini);
   	return 0;
} 
