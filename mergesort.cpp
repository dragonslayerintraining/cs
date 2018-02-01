#include <iostream>
#include <vector>
#include <algorithm>

//for debugging/testing purposes
void output(const std::vector<int>& vs){
  for(unsigned int i=0;i<vs.size();i++){
    if(i>0){
      std::cerr<<" ";
    }
    std::cerr<<vs[i];
  }
  std::cerr<<std::endl;
}

void mergesort(std::vector<int>& vs,int begin,int end){
  if(end-begin>1){
    int mid=(begin+end)/2;
    mergesort(vs,begin,mid);
    mergesort(vs,mid,end);
    int i=begin,j=mid;
    std::vector<int> tmp;
    while(i<mid||j<end){
      if(j==end||(i<mid&&vs[i]<=vs[j])){
	tmp.push_back(vs[i++]);
      }else{
	tmp.push_back(vs[j++]);
      }
    }
    std::copy(tmp.begin(),tmp.end(),vs.begin()+begin);
  }
}

void mergesort(std::vector<int>& vs){
  mergesort(vs,0,vs.size());
}

//tester
int main(){
  std::vector<int> vs;
  for(int i=0;i<100;i++){
    vs.push_back(i*17%100);
  }
  std::cout<<"Before: ";
  output(vs);
  mergesort(vs);
  std::cout<<"After: ";
  output(vs);
}
