#include <cstdio>
#include <chrono>
#include <random>
#include <climits>
#include <cassert>
#include <stdint.h>

const int64_t INF=1e9+7;

std::mt19937 rng(std::chrono::system_clock::now().time_since_epoch().count());

struct Treap{
  int priority;//max at root
  int64_t key,key_lazy;
  int val,val_min,val_lazy;
  int size;
  struct Treap* left,*right;
  static struct Treap* nil;
  Treap():priority(-1),
	  key(0),key_lazy(0),val(INF),val_min(INF),val_lazy(0),
	  size(0),left(nil),right(nil){
  }
  Treap(int64_t key,int val):priority(std::uniform_int_distribution<int>(0,INT_MAX)(rng)),
				    key(key),key_lazy(0),val(val),val_min(val),val_lazy(0),
				    size(1),left(nil),right(nil){
  }
  //warning: missing destructor
  void pull(){
    if(this!=nil){
      left->push();
      right->push();
      size=left->size+right->size+1;
      val_min=std::min(val,std::min(left->val_min,right->val_min));
    }
  }
  
  void push(){
    if(this!=nil){
      left->key_lazy+=key_lazy;
      right->key_lazy+=key_lazy;
      key+=key_lazy;
      key_lazy=0;
      
      left->val_lazy+=val_lazy;
      right->val_lazy+=val_lazy;
      val_min+=val_lazy;
      val+=val_lazy;
      val_lazy=0;
    }
  }
  
  static void dump(struct Treap* x){
    x->push();
    if(x!=nil){
      printf("[");
      dump(x->left);
      printf("%ld:%d",x->key,x->val);
      dump(x->right);
      printf("]");
    }
  }

  static void split(struct Treap* x,int k,struct Treap*& L,struct Treap*& R){
    assert(k>=0&&k<=x->size);
    x->push();
    L=R=nil;
    if(x==nil) return;
    if(k<=x->left->size){
      split(x->left,k,L,x->left);
      R=x;
    }else{
      split(x->right,k-x->left->size-1,x->right,R);
      L=x;
    }
    x->pull();
  }

  static struct Treap* concat(struct Treap* x,struct Treap* y){
    x->push();
    y->push();
    if(x==nil) return y;
    if(y==nil) return x;
    if(x->priority>y->priority){
      x->right=concat(x->right,y);
      x->pull();
      return x;
    }else{
      y->left=concat(x,y->left);
      y->pull();
      return y;
    }
  }

  static int lower_bound(struct Treap* x,int64_t key){
    x->push();
    if(x==nil) return 0;
    if(key<=x->key){
      return lower_bound(x->left,key);
    }else{
      return lower_bound(x->right,key)+x->left->size+1;
    }
  }

  static int upper_bound(struct Treap* x,int64_t key){
    x->push();
    if(x==nil) return 0;
    if(key<x->key){
      return upper_bound(x->left,key);
    }else{
      return upper_bound(x->right,key)+x->left->size+1;
    }
  }

  //returns index of leftmost minimum value
  static int argmin(struct Treap* x){
    if(x==nil) return 0;
    x->push();
    x->left->push();
    x->right->push();
    if(x->left!=nil&&x->val_min==x->left->val_min){
      return argmin(x->left);
    }else if(x->val_min==x->val){
      return x->left->size;
    }else if(x->right!=nil&&x->val_min==x->right->val_min){
      return argmin(x->right)+x->left->size+1;
    }else{
      assert(false);
    }
  }

  static struct Treap* kth(struct Treap* x,int k){
    assert(k>=0&&k<x->size);
    x->push();
    if(k<x->left->size){
      return kth(x->left,k);
    }else if(k>x->left->size){
      return kth(x->right,k-x->left->size-1);
    }else{
      return x;
    }
  }

  static struct Treap* extract(struct Treap*& x,int l,int r){
    struct Treap* L,*R;
    split(x,r,x,R);
    split(x,l,L,x);
    struct Treap* res=x;
    x=concat(L,R);
    return res;
  }
  static struct Treap* extract(struct Treap*& x,int l){
    return extract(x,l,l+1);
  }

  //warning: not stable
  static struct Treap* merge(struct Treap* x,struct Treap* y){
    struct Treap* res=nil;
    while(true){
      if(x==nil) return concat(res,y);
      if(y==nil) return concat(res,x);
      if(kth(x,0)->key>kth(y,0)->key){
	std::swap(x,y);
      }
      int k=upper_bound(x,kth(y,0)->key);
      assert(k>0);
      struct Treap* L;
      split(x,k,L,x);
      res=concat(res,L);
    }
  }
};
struct Treap* Treap::nil=new Treap();

int main(){
  struct Treap* root=Treap::nil;
  for(int i=0;i<10;i++){
    root=Treap::merge(root,new Treap(i,i));
  }
  Treap::dump(root);
  printf("\n");
  return 0;
}
