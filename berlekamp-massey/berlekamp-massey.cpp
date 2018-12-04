//Only works for prime modulo
#include <cstdio>
#include <vector>
#include <cassert>
#include <stdint.h>

const int64_t MOD=1e9+7;

int64_t modexp(int64_t base,int64_t exp){
  int64_t ac=1;
  for(;exp>0;exp>>=1){
    if(exp&1) ac=ac*base%MOD;
    base=base*base%MOD;
  }
  return ac;
}

int64_t inverse(int64_t x){
  return modexp(x,MOD-2);
}

void addmod(int64_t& x,int64_t y){
  x=(x+y)%MOD;
}

std::vector<int64_t> berlekamp_massey(const std::vector<int64_t>& seq){
  std::vector<int64_t> curr(1,1),prev(1,1);
  int64_t b=1;
  for(int64_t i=0,shift=1;i<seq.size();i++,shift++){
    int64_t d=0;
    for(int64_t j=0;j<curr.size();j++){
      addmod(d,curr[j]*seq[i-j]);
    }
    if(d==0) continue;
    std::vector<int64_t> next=(curr);
    if(prev.size()+shift>curr.size()){
      next.resize(prev.size()+shift);
    }
    int64_t scale=(MOD-d)*inverse(b)%MOD;
    for(int64_t j=0;j<prev.size();j++){
      addmod(next[j+shift],scale*prev[j]);
    }
    if(next.size()>curr.size()){
      prev=curr;
      b=d;
      shift=0;
    }
    curr=next;
  }
  return curr;
}

int main(){
  int64_t N,M;
  scanf("%ld %ld",&N,&M);
  std::vector<int64_t> seq;
  for(int64_t i=0;i<N;i++){
    int64_t x;
    scanf("%ld",&x);
    seq.push_back(x);
  }
  std::vector<int64_t> recur=berlekamp_massey(seq);
  for(int64_t i=seq.size();i<M;i++){
    int64_t ac=0;
    for(int64_t j=1;j<recur.size();j++){
      addmod(ac,recur[j]*seq[i-j]);
    }
    seq.push_back((MOD-ac)%MOD);
  }
  for(int64_t i=0;i<seq.size();i++){
    if(i) printf(" ");
    printf("%ld",seq[i]);
  }
  printf("\n");
  return 0;
}
