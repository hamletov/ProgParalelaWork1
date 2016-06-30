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

double memo[1<<20][21];
bool done[1<<20][21];

double dis[31][31];
int n;

double dp(int mask,int pos){
	if(__builtin_popcount(mask)==n)return dis[pos][0];
	if(done[mask][pos])return memo[mask][pos];
	double dev=1e+10;
	
	for(int i=0;i<n;i++)
		if( (mask&(1<<i))==0 )
			dev=min(dev, dis[pos][i]+ dp(mask|(1<<i),i) );
	
	memo[mask][pos]=dev;
	done[mask][pos]=1;
	return dev;
}

int main(){
	
	memset(done,0,sizeof(done));
	n=20;
	
	for(int i=0;i<n;i++)dis[i][i]=0;
	
	for(int i=0;i<n;i++)
		for(int j=i+1;j<n;j++){
			dis[i][j]=rand()%100;
			dis[j][i]=dis[i][j];
		}
	
	printf("%.10lf\n",dp(1,0));
	
    return 0;
}

