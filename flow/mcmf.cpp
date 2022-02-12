//Submitted to https://open.kattis.com/problems/mincostmaxflow
#include <cstdio>
#include <queue>
#include <vector>
#include <stdint.h>

const int64_t INF=1e9+7;

std::vector<int64_t> elist;
std::vector<int64_t> cap;
std::vector<int64_t> weight;
std::vector<int64_t> adj[250];

int64_t dist[250];
int64_t prev[250];


void add_he(int64_t node,int64_t C,int64_t W){
  adj[node].push_back(elist.size());
  elist.push_back(node);
  cap.push_back(C);
  weight.push_back(W);
}

int main(){
  int64_t N,M,S,T;
  scanf("%ld %ld %ld %ld",&N,&M,&S,&T);
  for(int64_t i=0;i<M;i++){
    int64_t U,V,C,W;
    scanf("%ld %ld %ld %ld",&U,&V,&C,&W);
    add_he(U,C,W);
    add_he(V,0,-W);
  }
  int64_t flow=0;
  int64_t cost=0;
  while(true){
    std::fill(dist,dist+N,INF);
    std::fill(prev,prev+N,-1);
    dist[S]=0;
    std::queue<int> frontier;
    frontier.push(S);
    while(!frontier.empty()){
      int node=frontier.front();
      frontier.pop();
      for(int64_t e:adj[node]){
	if(cap[e]){
	  int64_t U=elist[e],V=elist[e^1];
	  int64_t W=weight[e];
	  if(dist[V]>dist[U]+W){
	    dist[V]=dist[U]+W;
	    prev[V]=e;
	    frontier.push(V);
	  }
	}
      }
    }
    if(prev[T]==-1) break;
    int64_t aug=INF;
    for(int64_t node=T;node!=S;node=elist[prev[node]]){
      int64_t e=prev[node];
      aug=std::min(aug,cap[e]);
    }
    for(int64_t node=T;node!=S;node=elist[prev[node]]){
      int64_t e=prev[node];
      cap[e]-=aug;
      cap[e^1]+=aug;
      cost+=aug*weight[e];
    }
    flow+=aug;
  }
  printf("%ld %ld\n",flow,cost);
  return 0;
}
