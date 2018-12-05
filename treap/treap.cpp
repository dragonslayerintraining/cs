//Replace the random number generator if needed
//Based on https://codeforces.com/blog/entry/11148?#comment-162254
#include <cstdio>
#include <cstdlib>
#include <cassert>

struct Treap{
  struct Treap* left,*right;
  int priority;
  int size;
  //Data:
  int val;
  int lazy;
  Treap():left(this),right(this),priority(0),size(0){
  }
  Treap(int val);
  ~Treap();
  
  static void push(struct Treap* node){
    node->val+=node->lazy;
    node->left->lazy+=node->lazy;
    node->right->lazy+=node->lazy;
    node->lazy=0;
  }
  static void pull(struct Treap* node){
    node->size=node->left->size+node->right->size+1;
  }
  
  static void split(struct Treap* node,struct Treap*& L,struct Treap*& R,int k);
  static struct Treap* merge(struct Treap* L,struct Treap* R);
  static int kth(struct Treap* node,int k);
  static void dump(struct Treap* node);
  static struct Treap* increase(struct Treap* node,int l,int r,int x);
}* nil=new Treap();

Treap::Treap(int val):left(nil),right(nil),priority(rand()),size(1),val(val),lazy(0){
}

Treap::~Treap(){
  if(left!=nil) delete left;
  if(right!=nil) delete right;
}

void Treap::split(struct Treap* node,struct Treap*& L,struct Treap*& R,int k){
  L=R=nil;
  if(node==nil) return;
  push(node);
  if(k<=node->left->size){
    split(node->left,L,node->left,k);
    R=node;
  }else{
    split(node->right,node->right,R,k-node->left->size-1);
    L=node;
  }
  pull(node);
}

struct Treap* Treap::merge(struct Treap* L,struct Treap* R){
  if(L==nil) return R;
  if(R==nil) return L;
  if(L->priority>R->priority){
    push(L);
    L->right=merge(L->right,R);
    pull(L);
    return L;
  }else{
    push(R);
    R->left=merge(L,R->left);
    pull(R);
    return R;
  }
}

int Treap::kth(struct Treap* node,int k){
  assert(k>=0&&k<node->size);
  push(node);
  if(k<node->left->size){
    return kth(node->left,k);
  }else if(k>node->left->size){
    return kth(node->right,k-node->left->size-1);
  }else{
    return node->val;
  }
}

//For debugging:
void Treap::dump(struct Treap* node){
  if(node==nil) return;
  printf("[");
  dump(node->left);
  printf("%d(%d)",node->val,node->lazy);
  dump(node->right);
  printf("]");
}

//[l,r)
struct Treap* Treap::increase(struct Treap* node,int l,int r,int x){
  struct Treap* L,*R;
  split(node,node,R,r);
  split(node,L,node,l);
  node->lazy+=x;
  return merge(L,merge(node,R));
}

void show(struct Treap* node){
  for(int i=0;i<node->size;i++){
    if(i>0) printf(" ");
    printf("%d",Treap::kth(node,i));
  }
  printf("\n");
}

int main(){
  struct Treap* root=nil;
  for(int i=0;i<10;i++){
    root=Treap::merge(root,new Treap(i));
  }
  for(int i=0;i<root->size;i+=2){
    root=Treap::increase(root,i,i+2,10);
    show(root);
  }
  show(root);
  delete root;
  root=nil;
  for(int i=0;i<1000000;i++){
    root=Treap::merge(root,new Treap(i));
  }
  for(int i=0;i<root->size;i++){
    struct Treap* left,*right;
    Treap::split(root,left,right,rand()%root->size);
    root=Treap::merge(right,left);
  }
  show(root);
  delete root;
  return 0;
}
