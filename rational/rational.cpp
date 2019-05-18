#include <cstdio>
#include <stdint.h>
#include <cassert>

struct Ratio{
  int64_t p,q;
  Ratio(int64_t p,int64_t q):p(p),q(q){
  }
  bool operator<(struct Ratio r)const{
    return p*r.q<q*r.p;
  }
};

//weighted mediant a*x+b
struct Ratio mediant(struct Ratio a,int64_t x,struct Ratio b){
  return Ratio(a.p*x+b.p,a.q*x+b.q);
}

//returns index of last 0
template<class T>
int64_t gallop(T func){
  int64_t low=0,high=1;
  while(!func(high)){
    high*=2;
  }
  while(high-low>1){
    int64_t mid=(low+high)/2;
    if(func(mid)){
      high=mid;
    }else{
      low=mid;
    }
  }
  return low;
}

//Finds simplest ratio in (low,high)
struct Ratio between(struct Ratio low,struct Ratio high){
  assert(low<high);
  struct Ratio lb(0,1),ub(1,0);
  while(true){
    ub=mediant(lb,gallop([&](int64_t x){return mediant(lb,x,ub)<high;}),ub);
    lb=mediant(ub,gallop([&](int64_t x){return low<mediant(ub,x,lb);}),lb);
    struct Ratio mt=mediant(lb,1,ub);
    if(low<mt&&mt<high){
      return mt;
    }
  }
}

int main(){
  int64_t A,B,C,D;
  scanf("%ld %ld %ld %ld",&A,&B,&C,&D);
  struct Ratio r=between(Ratio(A,B),Ratio(C,D));
  printf("%ld/%ld\n",r.p,r.q);
  return 0;
}
