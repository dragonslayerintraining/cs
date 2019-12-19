//Implementation of online fully dynamic convex hull
//Currently only for one side but it should be easy to add the other side
//O(log^2n) per operation
//Based on paper "Maintenance of configurations in the plane" by Overmars and van Leeuwen
//Using join-based treaps as BST
//Assumes all points are always distinct
//Assumes coordinates are at most 10^6

//TODO: Benchmark (without asserts)
//      Is this faster than naively recomputing convex hull?
//TODO: We don't need to allocate a new Treap<Point> to delete point
//TODO: Clean up code

#include <cstdio>
#include <random>
#include <chrono>
#include <fstream>
#include <set>
#include <cassert>

std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

struct Point{
  int64_t x,y;
  Point operator-(Point p)const{
    return {x-p.x,y-p.y};
  }
  int64_t cross(Point p)const{
    return x*p.y-y*p.x;
  }
  int64_t dot(Point p)const{
    return x*p.x+y*p.y;
  }
  bool operator<(Point p)const{
    if(y!=p.y) return y<p.y;
    return x<p.x;
  }
};

int64_t cross(Point a,Point b,Point c){
  return (b-a).cross(c-a);
}

//Warning: memory leak if deleting Treap with more than one element
template<class T>
struct Treap{
  Treap* left,*right;
  Treap* leftmost;
  decltype(rng()) priority;
  int64_t size;
  T data;
  Treap(T data):left(NULL),right(NULL),leftmost(this),priority(rng()),size(1),data(data){
  }
  static int64_t get_size(Treap* x){
    return x?x->size:0;
  }
  void pull(){
    size=1+get_size(left)+get_size(right);
    leftmost=left?left->leftmost:this;
  }
  void push(){
  }
  static Treap* concat(Treap* x,Treap* y){
    if(x==NULL) return y;
    if(y==NULL) return x;
    if(x->priority<y->priority){
      y->push();
      y->left=concat(x,y->left);
      y->pull();
      return y;
    }else{
      x->push();
      x->right=concat(x->right,y);
      x->pull();
      return x;
    }
  }
  static void split(Treap* x,int64_t k,Treap*& L,Treap*& R){
    assert(k>=0&&k<=get_size(x));
    L=R=NULL;
    if(x==NULL) return;
    x->push();
    if(k<=get_size(x->left)){
      split(x->left,k,L,x->left);
      R=x;
    }else{
      split(x->right,k-get_size(x->left)-1,x->right,R);
      L=x;
    }
    x->pull();
  }
  static int64_t lower_bound(Treap* x,T key){
    if(x==NULL) return 0;
    if(!(x->data<key)){
      return lower_bound(x->left,key);
    }else{
      return lower_bound(x->right,key)+get_size(x->left)+1;
    }
  }
  static Treap* insert(Treap* x,T key){
    int64_t k=lower_bound(x,key);
    Treap* left,*right;
    split(x,k,left,right);
    return concat(concat(left,new Treap(key)),right);
  }
  static Treap* erase(Treap* x,T key){
    int64_t k=lower_bound(x,key);
    Treap* left,*right;
    split(x,k,left,x);
    assert(x);
    split(x,1,x,right);
    assert(!(x->data<key)&&!(key<x->data));
    delete x;
    return concat(left,right);
  }
  static T kth(Treap* x,int64_t k){
    assert(x);
    assert(k>=0&&k<get_size(x));
    if(k<get_size(x->left)){
      return kth(x->left,k);
    }else if(k>get_size(x->left)){
      return kth(x->right,k-get_size(x->left)-1);
    }else{
      return x->data;
    }
  }
};

struct Walker{
  Treap<Point>* m,*extra;
  int index;
  Walker(Treap<Point>* root):m(root),extra(NULL),index(0){
    assert(root);
    canonicalize();
  }
  void go_left(){
    assert(m);
    extra=m;
    m=m->left;
    canonicalize();
  }
  void go_right(){
    assert(m);
    index+=Treap<Point>::get_size(m->left)+1;
    m=m->right;
    canonicalize();
  }
  //ensure m->right or extra
  void canonicalize(){
    assert(m||extra);
    if(!extra&&!m->right){
      go_left();
      assert(extra);
    }
  }
  bool unique(){
    return !m;
  }
  Point ml(){
    if(!m) return extra->data;
    return m->data;
  }
  Point mr(){
    if(!m) return extra->data;
    assert(m->right||extra);
    return m->right?m->right->leftmost->data:extra->data;
  }
};

std::pair<int64_t,int64_t> find_bridge(Treap<Point>* left,Treap<Point>* right){
  assert(left);
  assert(right);
  assert(Treap<Point>::kth(left,Treap<Point>::get_size(left)-1)<Treap<Point>::kth(right,0));
  struct Walker w1(left),w2(right);
  int64_t split_y=Treap<Point>::kth(right,0).y;
  while(!w1.unique()||!w2.unique()){
    Point b=w1.mr(),c=w2.ml();
    //Change cross(...)>0 to cross(...)>=0 to ignore points on edges
    if(!w1.unique()&&cross(w1.ml(),b,c)>0){
      w1.go_left();
    }else if(!w2.unique()&&cross(b,c,w2.mr())>0){
      w2.go_right();
    }else if(w1.unique()){
      w2.go_left();
    }else if(w2.unique()){
      w1.go_right();
    }else{
      Point a=w1.ml(),d=w2.mr();
      int64_t s1=cross(a,b,c);
      int64_t s2=cross(b,a,d);
      assert(s1+s2>=0);
      //y=(s1*d.y+s2*c.y)/(s1+s2);
      if(s1+s2==0||s1*d.y+s2*c.y<split_y*(s1+s2)){
	w1.go_right();
      }else{
	w2.go_left();
      }
    }
  }
  int l1=w1.index,r1=w2.index;
  assert(l1==0||cross(Treap<Point>::kth(left,l1-1),Treap<Point>::kth(left,l1),Treap<Point>::kth(right,r1))<=0);
  assert(r1==right->size-1||cross(Treap<Point>::kth(left,l1),Treap<Point>::kth(right,r1),Treap<Point>::kth(right,r1+1))<=0);
  assert(l1==left->size-1||cross(Treap<Point>::kth(left,l1),Treap<Point>::kth(left,l1+1),Treap<Point>::kth(right,r1))>=0);
  assert(r1==0||cross(Treap<Point>::kth(left,l1),Treap<Point>::kth(right,r1-1),Treap<Point>::kth(right,r1))>=0);
  return {l1,r1};
}

Treap<Point>* concat_hulls(Treap<Point>* left,Treap<Point>* right,Treap<Point>*& save1,Treap<Point>*& save2,int64_t& cut){
  cut=Treap<Point>::get_size(left);
  if(left==NULL) return right;
  if(right==NULL) return left;
  std::pair<int64_t,int64_t> bridge=find_bridge(left,right);
  cut=bridge.first+1;
  Treap<Point>::split(left,bridge.first+1,left,save1);
  Treap<Point>::split(right,bridge.second,save2,right);
  return Treap<Point>::concat(left,right);
}

void unconcat_hulls(Treap<Point>* x,Treap<Point>*& left,Treap<Point>*& right,Treap<Point>* save1,Treap<Point>* save2,int64_t cut){
  Treap<Point>::split(x,cut,left,right);
  left=Treap<Point>::concat(left,save1);
  right=Treap<Point>::concat(save2,right);
}

struct HullData{
  Point p;
  Treap<Point>* hull;
  Treap<Point>* save1;//from merging p and R
  Treap<Point>* save2,*save3;//from merging L and p+R
  int64_t cut1;
  HullData(Point p):p(p),hull(new Treap<Point>(p)),save1(NULL),save2(NULL),save3(NULL){
  }
  HullData(const HullData& h):HullData(h.p){
    assert(Treap<Point>::get_size(h.hull)==1);
  }
  ~HullData(){
    assert(Treap<Point>::get_size(hull)==1);
    delete hull;
  }
  bool operator<(HullData h)const{
    return p<h.p;
  }
};

template<>
void Treap<HullData>::pull(){
  size=1+get_size(left)+get_size(right);
  if(right){
    Treap<Point>* save0;
    int64_t cut2;
    data.hull=concat_hulls(data.hull,right->data.hull,save0,data.save1,cut2);
    assert(save0==NULL);
    assert(cut2==1);
  }
  if(left){
    data.hull=concat_hulls(left->data.hull,data.hull,data.save2,data.save3,data.cut1);
  }
}

template<>
void Treap<HullData>::push(){
  if(left){
    unconcat_hulls(data.hull,left->data.hull,data.hull,data.save2,data.save3,data.cut1);
  }
  if(right){
    unconcat_hulls(data.hull,data.hull,right->data.hull,NULL,data.save1,1);
  }
}


//===Testing code begins here===
void log_point(Point p){
  std::ofstream fout;
  fout.open("geom.txt",std::ios::app);
  fout<<"P "<<p.x<<" "<<p.y<<std::endl;
}

void log_seg(Point p,Point q){
  std::ofstream fout;
  fout.open("geom.txt",std::ios::app);
  fout<<"L "<<p.x<<" "<<p.y<<" "<<q.x<<" "<<q.y<<std::endl;
}

int main(){
  Treap<HullData>* hull=NULL;
  if(0){
    std::set<Point> points;
    for(int64_t i=0;i<100000;i++){
      Point p{rng()%100000,rng()%100000};
      if(points.count(p)) continue;
      points.insert(p);
      log_point(p);
      hull=Treap<HullData>::insert(hull,HullData(p));
    }
  }else{
    for(int64_t i=0;i<10;i++){
      for(int64_t j=0;j<10;j++){
	log_point(Point{i,j});
	hull=Treap<HullData>::insert(hull,HullData(Point{i,j}));
      }
    }
  }
  while(hull){
    Treap<Point>* root=hull->data.hull;
    for(int64_t i=0;i<Treap<Point>::get_size(root)-1;i++){
      log_seg(Treap<Point>::kth(root,i),Treap<Point>::kth(root,i+1));
    }
    std::vector<Point> ps;
    for(int64_t i=0;i<Treap<Point>::get_size(root);i++){
      ps.push_back(Treap<Point>::kth(root,i));
    }
    for(Point p:ps){
      hull=Treap<HullData>::erase(hull,HullData(p));
    }
  }
}
