#include <cstdio>
#include <cassert>

//Counts the number of integer pairs (x,y) that satisfy
//x>0
//y>0
//ax+by<=c
//Only works if the number of such points is finite

int lattice(int a,int b,int c){
  return (c<a+b)?0:lattice(b,a%b,c-a/b*c/a*b)+(a/b*c/a)*(c/a)-(c/a)*(c/a+1)/2*(a/b);
}

int brute(int a,int b,int c){
  int cnt=0;
  for(int i=1;a*i<=c;i++){
    for(int j=1;a*i+b*j<=c;j++){
      cnt++;
    }
  }
  return cnt;
}

int main(){
  for(int i=1;i<=200;i++){
    for(int j=1;j<=200;j++){
      for(int k=1;k<=200;k++){
	assert(lattice(i,j,k)==brute(i,j,k));
      }
    }
  }
  return 0;
}
