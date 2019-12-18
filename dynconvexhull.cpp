//Implementation of online fully dynamic convex hull
//Currently only for one side but it should be easy to add the other side
//O(log^3n) per operation
//Based on paper "Maintenance of configurations in the plane" by Overmars and van Leeuwen
//Assumes all points are always distinct
//Assumes coordinates are at most 10^6

//TODO: Benchmark (without asserts)
//      Is this faster than naively recomputing convex hull?
//TODO: One of the saves is always NULL, remove it. Similar for a cut.
//TODO: We don't need to allocate a new Treap<Point> to delete point
//TODO: Clean up code
//TODO: Remove extra log from binary search

#include <cstdio>
#include <random>
#include <chrono>
#include <cassert>
#include <fstream>
#include <set>


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

template<class T>
struct Treap{
  Treap* left,*right;
  int64_t priority;
  int64_t size;
  T data;
  Treap(T data):left(NULL),right(NULL),priority(rng()),size(1),data(data){
  }
  static int64_t get_size(Treap* x){
    return x?x->size:0;
  }
  Treap* pull(){
    size=1+get_size(left)+get_size(right);
    return this;
  }
  void push(){
  }
  static Treap* concat(Treap* x,Treap* y){
    if(x==NULL) return y;
    if(y==NULL) return x;
    if(x->priority<y->priority){
      y->push();
      y->left=concat(x,y->left);
      return y->pull();
    }else{
      x->push();
      x->right=concat(x->right,y);
      return x->pull();
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
  static std::pair<int64_t,int64_t> bridge(Treap<Point>* left,Treap<Point>* right){
    assert(left);
    assert(right);
    assert(kth(left,get_size(left)-1)<kth(right,0));
    int64_t l1=0,r1=get_size(left)-1;
    int64_t l2=0,r2=get_size(right)-1;
    int64_t split_y=kth(right,0).y;
    while(l1<r1||l2<r2){
      int64_t m1=(l1+r1+1)/2;
      int64_t m2=(l2+r2)/2;
      Point b=kth(left,m1),c=kth(right,m2);
      //Change cross(...)>0 to cross(...)>=0 to ignore points on edges
      if(l1<r1&&cross(kth(left,m1-1),b,c)>0){
	r1=m1-1;
      }else if(l2<r2&&cross(b,c,kth(right,m2+1))>0){
	l2=m2+1;
      }else if(l1==r1){
	r2=m2;
      }else if(l2==r2){
	l1=m1;
      }else{
	Point a=kth(left,m1-1),d=kth(right,m2+1);
	int64_t s1=cross(a,b,c);
	int64_t s2=cross(b,a,d);
	if(s1+s2==0){
	  l1=m1;
	  r2=m2;
	  continue;
	}
	assert(s1+s2>0);
	//y=(s1*d.y+s2*c.y)/(s1+s2);
	if(s1*d.y+s2*c.y<split_y*(s1+s2)){
	  l1=m1;
	}else{
	  r2=m2;
	}
      }
    }
    assert(l1==r1);
    assert(l2==r2);
    assert(l1==0||cross(kth(left,l1-1),kth(left,l1),kth(right,r2))<=0);
    assert(r2==right->size-1||cross(kth(left,l1),kth(right,r2),kth(right,r2+1))<=0);
    assert(l1==left->size-1||cross(kth(left,l1),kth(left,l1+1),kth(right,r2))>=0);
    assert(r2==0||cross(kth(left,l1),kth(right,r2-1),kth(right,r2))>=0);
    return {l1,r2};
  }
  static Treap<Point>* concat_hulls(Treap<Point>* left,Treap<Point>* right,Treap<Point>*& save1,Treap<Point>*& save2,int64_t& cut){
    cut=get_size(left);
    if(left==NULL) return right;
    if(right==NULL) return left;
    std::pair<int64_t,int64_t> b=bridge(left,right);
    cut=b.first+1;
    split(left,b.first+1,left,save1);
    split(right,b.second,save2,right);
    return concat(left,right);
  }
  static void unconcat_hulls(Treap<Point>* x,Treap<Point>*& left,Treap<Point>*& right,Treap<Point>* save1,Treap<Point>* save2,int64_t cut){
    split(x,cut,left,right);
    left=concat(left,save1);
    right=concat(save2,right);
  }
};

struct HullData{
  Point p;
  Treap<Point>* hull;
  Treap<Point>* save1,*save2;//from merging L and M
  Treap<Point>* save3,*save4;//from merging L+M and R
  int64_t cut1,cut2;
  HullData(Point p):p(p),hull(new Treap<Point>(p)),save1(NULL),save2(NULL),save3(NULL),save4(NULL){
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
Treap<HullData>* Treap<HullData>::pull(){
  size=1+get_size(left)+get_size(right);
  if(left){
    data.hull=Treap<Point>::concat_hulls(left->data.hull,data.hull,data.save1,data.save2,data.cut1);
  }
  if(right){
    data.hull=Treap<Point>::concat_hulls(data.hull,right->data.hull,data.save3,data.save4,data.cut2);
  }
  return this;
}

template<>
void Treap<HullData>::push(){
  if(right){
    Treap<Point>::unconcat_hulls(data.hull,data.hull,right->data.hull,data.save3,data.save4,data.cut2);
  }
  if(left){
    Treap<Point>::unconcat_hulls(data.hull,left->data.hull,data.hull,data.save1,data.save2,data.cut1);
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
  if(0){
    Treap<Point>* root1=NULL;
    for(int64_t i=0;i<5;i++){
      root1=Treap<Point>::concat(root1,new Treap<Point>({(i-2)*(i-2),i}));
    }
    Treap<Point>* root2=NULL;
    for(int64_t i=0;i<5;i++){
      root2=Treap<Point>::concat(root2,new Treap<Point>({(i-2)*(i-2),i+5}));
    }
    for(int64_t i=0;i<root1->size;i++){
      log_point(Treap<Point>::kth(root1,i));
    }
    for(int64_t i=0;i<root2->size;i++){
      log_point(Treap<Point>::kth(root2,i));
    }
    std::pair<int64_t,int64_t> bridge=Treap<Point>::bridge(root1,root2);
    log_seg(Treap<Point>::kth(root1,bridge.first),
	    Treap<Point>::kth(root2,bridge.second));
  }
  if(0){
    Treap<Point>* root1=NULL;
    for(int64_t i=0;i<5;i++){
      root1=Treap<Point>::concat(root1,new Treap<Point>({(i-4)*(i-4)+5,i}));
    }
    Treap<Point>* root2=NULL;
    for(int64_t i=0;i<5;i++){
      root2=Treap<Point>::concat(root2,new Treap<Point>({(i-0)*(i-0),i+5}));
    }
    for(int64_t i=0;i<root1->size;i++){
      log_point(Treap<Point>::kth(root1,i));
    }
    for(int64_t i=0;i<root2->size;i++){
      log_point(Treap<Point>::kth(root2,i));
    }
    std::pair<int64_t,int64_t> bridge=Treap<Point>::bridge(root1,root2);
    log_seg(Treap<Point>::kth(root1,bridge.first),
	    Treap<Point>::kth(root2,bridge.second));
  }
  if(0){
    Treap<Point>* root=new Treap<Point>(Point{-1,-1});
    Treap<Point>* ign1,*ign2;
    for(int64_t i=0;i<20;i++){
      Point p{rng()%1000,i};
      log_point(p);
      int64_t cut;
      Treap<Point>* add=new Treap<Point>({p});
      root=Treap<Point>::concat_hulls(root,add,ign1,ign2,cut);
      Treap<Point>::unconcat_hulls(root,root,add,ign1,ign2,cut);
      root=Treap<Point>::concat_hulls(root,add,ign1,ign2,cut);
    }
    for(int64_t i=0;i<root->size-1;i++){
      log_seg(Treap<Point>::kth(root,i),Treap<Point>::kth(root,i+1));
    }
  }
  if(0){
    Treap<HullData>* hull=NULL;
    for(int64_t i=0;i<20;i++){
      Point p{rng()%1000,i};
      log_point(p);
      hull=Treap<HullData>::concat(hull,new Treap<HullData>(p));
    }
    Treap<Point>* root=hull->data.hull;
    for(int64_t i=0;i<root->size-1;i++){
      log_seg(Treap<Point>::kth(root,i),Treap<Point>::kth(root,i+1));
    }
  }
  if(1){
    Treap<HullData>* hull=NULL;
    std::set<Point> points;
    for(int64_t i=0;i<100000;i++){
      Point p{rng()%1000,rng()%1000};
      if(points.count(p)) continue;
      points.insert(p);
      //log_point(p);
      hull=Treap<HullData>::insert(hull,HullData(p));
    }
    while(hull){
      Treap<Point>* root=hull->data.hull;
      for(int64_t i=0;i<Treap<Point>::get_size(root)-1;i++){
	//log_seg(Treap<Point>::kth(root,i),Treap<Point>::kth(root,i+1));
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
  if(0){
    Treap<HullData>* hull=NULL;
    for(int64_t i=0;i<10;i++){
      for(int64_t j=0;j<10;j++){
	log_point(Point{i,j});
	hull=Treap<HullData>::insert(hull,HullData(Point{i,j}));
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
}
