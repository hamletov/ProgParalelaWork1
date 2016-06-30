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
bool visited[31];

double dis[31][31];
int n;

double dfs(int pos,int cont){
	if(cont==n)return dis[pos][0];
	
	double dev=1e+10;
	for(int i=0;i<n;i++)
		if(!visited[i]){
			visited[i]=1;
			dev=min(dev, dfs(i,cont+1)+dis[pos][i] );
			visited[i]=0;
		}
	return dev;
}

int main(){
	
	n=12;
	
	for(int i=0;i<n;i++)dis[i][i]=0;
	
	for(int i=0;i<n;i++)
		for(int j=i+1;j<n;j++){
			dis[i][j]=rand()%100;
			dis[j][i]=dis[i][j];
		}
	
	memset(visited,0,sizeof(visited));
	visited[0]=1;
	printf("%.10lf\n",dfs(0,1));
	
    return 0;
}

