//Tested on https://codeforces.com/contest/997/problem/E
#include <cstdio>
#include <vector>
#include <numeric>
#include <array>

const int INF=1e9+7;

struct MinCntAccSegtree{
  struct Node{
    int min,cnt;
    long long sum;
    Node():min(INF),cnt(0),sum(0){
    }
    Node(int min,int cnt,long long sum):min(min),cnt(cnt),sum(sum){
    }
  };
private:
  int N;
  Node combine(Node a,Node b){
    int min=std::min(a.min,b.min);
    return Node(min,(a.min==min)*a.cnt+(b.min==min)*b.cnt,a.sum+b.sum);
  }
  std::vector<Node> st;
  std::vector<int> lazy_inc,lazy_acc;
  void pull(int w,int L,int R){
    st[w]=combine(st[w<<1],st[w<<1|1]);
  }
  void push2(int w,int L,int R){
    //apply inc
    if(R-L>1){
      lazy_inc[w<<1]+=lazy_inc[w];
      lazy_inc[w<<1|1]+=lazy_inc[w];
    }
    st[w].min+=lazy_inc[w];
    lazy_inc[w]=0;
  }
  void push(int w,int L,int R){
    push2(w,L,R);
    //apply acc
    if(R-L>1){
      int M=(L+R)/2;
      push2(w<<1,L,M);
      push2(w<<1|1,M,R);
      if(st[w].min==st[w<<1].min){
	lazy_acc[w<<1]+=lazy_acc[w];
      }
      if(st[w].min==st[w<<1|1].min){
	lazy_acc[w<<1|1]+=lazy_acc[w];
      }
    }
    st[w].sum+=st[w].cnt*lazy_acc[w];
    lazy_acc[w]=0;
  }
  void build(int w,int L,int R){
    if(R-L==1){
      st[w]=Node(INF,1,0);//INF means disabled
    }else{
      int M=(L+R)/2;
      build(w<<1,L,M);
      build(w<<1|1,M,R);
      pull(w,L,R);
    }
  }
  void update_inc(int w,int L,int R,int a,int b,int v){
    push(w,L,R);
    if(a>=R||b<=L) return;
    if(a<=L&&b>=R){
      lazy_inc[w]+=v;
      push(w,L,R);
    }else{
      int M=(L+R)/2;
      update_inc(w<<1,L,M,a,b,v);
      update_inc(w<<1|1,M,R,a,b,v);
      pull(w,L,R);
    }
  }
  void update_acc(int w,int L,int R,int a,int b,int v){
    push(w,L,R);
    if(a>=R||b<=L) return;
    if(a<=L&&b>=R){
      if(st[w].min==0){
	lazy_acc[w]+=v;
      }
      push(w,L,R);
    }else{
      int M=(L+R)/2;
      update_acc(w<<1,L,M,a,b,v);
      update_acc(w<<1|1,M,R,a,b,v);
      pull(w,L,R);
    }
  }
  Node query(int w,int L,int R,int a,int b){
    push(w,L,R);
    if(a>=R||b<=L) return Node();
    if(a<=L&&b>=R){
      return st[w];
    }else{
      int M=(L+R)/2;
      return combine(query(w<<1,L,M,a,b),
		     query(w<<1|1,M,R,a,b));
    }
  }
		     
public:
  MinCntAccSegtree(int N):N(N),st(N*4),lazy_inc(N*4),lazy_acc(N*4){
    build(1,0,N);
  }
  void update_inc(int a,int b,int v){
    update_inc(1,0,N,a,b,v);
  }
  void update_acc(int a,int b,int v){
    update_acc(1,0,N,a,b,v);
  }
  long long query_sum(int a,int b){
    return query(1,0,N,a,b).sum;
  }
  Node query(int i){
    return query(1,0,N,i,i+1);
  }
};
