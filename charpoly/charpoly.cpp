#include <cstdio>
#include <algorithm>
#include <vector>

//Computes the characteristic polynomial in O(n^3) arithmetic operations
//Based on discussion ar https://codeforces.com/blog/entry/92248?#comment-809408

//Tested at https://judge.yosupo.jp/problem/characteristic_polynomial
//Submission: https://judge.yosupo.jp/submission/53682

const int MOD=998244353;

int mat[500][500];
struct Poly{
  std::vector<int> poly;
  Poly():poly({0}){
  }
  Poly(std::vector<int> poly):poly(poly){
  }
  friend Poly operator*(Poly a,Poly b){
    std::vector<int> res(a.size()+b.size()-1);
    for(size_t i=0;i<a.size();i++){
      for(size_t j=0;j<b.size();j++){
	res[i+j]=(res[i+j]+1LL*a[i]*b[j])%MOD;
      }
    }
    return Poly{res};
  }
  friend Poly operator+(Poly a,Poly b){
    std::vector<int> res(std::max(a.size(),b.size()));
    for(size_t i=0;i<a.size();i++){
      res[i]=(res[i]+a[i])%MOD;
    }
    for(size_t i=0;i<b.size();i++){
      res[i]=(res[i]+b[i])%MOD;
    }
    return Poly{res};
  }
  friend Poly operator*(int scale,Poly a){
    std::vector<int> res(a.size());
    for(size_t i=0;i<a.size();i++){
      res[i]=(1LL*a[i]*scale)%MOD;
    }
    return Poly{res};
  }
  int operator[](int index){
    return poly[index];
  }
  size_t size(){
    return poly.size();
  }
}det[500][500];//det[i][j]: determinant of matrix formed by rows 0..i (inclusive) and columns 0..i-1 and j

int N;

int modexp(int base,int exp){
  int ac=1;
  for(;exp>0;exp>>=1){
    if(exp&1) ac=1LL*ac*base%MOD;
    base=1LL*base*base%MOD;
  }
  return ac;
}

int inverse(int x){
  return modexp(x,MOD-2);
}

int main(){
  scanf("%d",&N);
  if(N==0){
    printf("1\n");
    return 0;
  }
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      scanf("%d",&mat[i][j]);
    }
  }
  //Convert to Hessenberg form
  for(int i=0;i<N;i++){
    //process ith column
    //find nonzero in column
    for(int j=i+1;j<N;j++){
      if(mat[j][i]!=0){
	for(int k=0;k<N;k++){
	  std::swap(mat[j][k],mat[i+1][k]);
	}
	for(int k=0;k<N;k++){
	  std::swap(mat[k][j],mat[k][i+1]);
	}
	break;
      }
    }
    if(mat[i+1][i]==0) continue;
    //make subdiagonal element one
    int inv=mat[i+1][i];
    int scale=inverse(inv);
    for(int k=0;k<N;k++){
      mat[i+1][k]=1LL*mat[i+1][k]*scale%MOD;
    }
    for(int k=0;k<N;k++){
      mat[k][i+1]=1LL*mat[k][i+1]*inv%MOD;
    }
    //zero everything below subdiagonal
    for(int j=i+2;j<N;j++){
      int scale=mat[j][i];
      for(int k=0;k<N;k++){
	mat[j][k]=(mat[j][k]-1LL*scale*mat[i+1][k])%MOD;
      }
      for(int k=0;k<N;k++){
	mat[k][i+1]=(mat[k][i+1]+1LL*scale*mat[k][j])%MOD;
      }
    }
  }
  //Laplace expansion by DP
  for(int i=0;i<N;i++){
    for(int j=i;j<N;j++){
      Poly p=(i==j)?Poly({-mat[i][j],1}):Poly({-mat[i][j]});
      det[i][j]=(i==0)?p:p*det[i-1][i-1]+mat[i][i-1]*det[i-1][j];
    }
  }
  for(int k=0;k<=N;k++){
    if(k) printf(" ");
    printf("%d",(det[N-1][N-1][k]+MOD)%MOD);
  }
  printf("\n");
}
