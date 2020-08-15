#include <cstdio>
#include <vector>
#include <cassert>
#include <cstdint>
#include <functional>

//A bitset implementation thats support efficient range operations

uint64_t get_bit(const uint64_t* data,size_t index){
  return (data[index/64]>>(index%64))&1;
}

void set_bit(uint64_t* data,size_t index,uint64_t val){
  data[index/64]=(data[index/64]&~(1ULL<<(index%64)))|(val<<(index%64));
}

//not in the main loop so doesn't need to be super optimized
uint64_t get_word_part(const uint64_t* data,size_t index,size_t len){
  assert(len>=0&&len<=64);
  uint64_t word=0;
  for(size_t k=0;k<len;k++){
    word|=get_bit(data,index+k)<<k;
  }
  return word;
}

void set_word_part(uint64_t* data,size_t index,size_t len,uint64_t val){
  assert(len>=0&&len<=64);
  for(size_t k=0;k<len;k++){
    set_bit(data,index+k,(val>>k)&1);
  }
}

struct op_second{
  uint64_t operator()(uint64_t x,uint64_t y){
    return y;
  }
};

//ranges may overlap
template<class Combine=op_second>
void bitrange_assign(uint64_t* dest,size_t dest_begin,uint64_t* src,size_t src_begin,size_t len){
  Combine op;
  const size_t dest_end=dest_begin+len;
  const size_t src_end=src_begin+len;
  if(len<=64){
    set_word_part(dest,dest_begin,len,op(get_word_part(dest,dest_begin,len),get_word_part(src,src_begin,len)));
  }else{
    size_t mid1=(dest_begin+63)/64;
    size_t mid2=(dest_end)/64;
    size_t head_len=mid1*64-dest_begin;
    size_t tail_len=dest_end-mid2*64;
    assert(mid1<=mid2);

    if(src==dest&&dest_begin<src_begin){
      //copy forwards
      set_word_part(dest,dest_begin,head_len,op(get_word_part(dest,dest_begin,head_len),get_word_part(src,src_begin,head_len)));
      size_t index=src_begin+head_len;
      size_t j=index/64,shift=index%64;
      if(shift==0){
	for(size_t i=mid1;i<mid2;i++){
	  dest[i]=op(dest[i],src[j]);
	  j++;
	}
      }else{
	for(size_t i=mid1;i<mid2;i++){
	  dest[i]=op(dest[i],(src[j]>>shift)|(src[j+1]<<(64-shift)));
	  j++;
	}
      }
      set_word_part(dest,dest_end-tail_len,tail_len,op(get_word_part(dest,dest_end-tail_len,tail_len),get_word_part(src,src_end-tail_len,tail_len)));
    }else{
      //copy backwards
      set_word_part(dest,dest_end-tail_len,tail_len,op(get_word_part(dest,dest_end-tail_len,tail_len),get_word_part(src,src_end-tail_len,tail_len)));
      size_t index=src_end-tail_len;
      size_t j=index/64,shift=index%64;
      if(shift==0){
	for(size_t i=mid2;i>mid1;i--){
	  dest[i-1]=op(dest[i-1],src[j-1]);
	  j--;
	}
      }else{
	for(size_t i=mid2;i>mid1;i--){
	  dest[i-1]=op(dest[i-1],(src[j-1]>>shift)|(src[j]<<(64-shift)));
	  j--;
	}
      }
      set_word_part(dest,dest_begin,head_len,op(get_word_part(dest,dest_begin,head_len),get_word_part(src,src_begin,head_len)));
    }
  }
}

size_t bitrange_popcount(uint64_t* data,size_t begin,size_t len){
  const size_t end=begin+len;
  if(len<=64){
    return __builtin_popcountll(get_word_part(data,begin,len));
  }else{
    size_t mid1=(begin+63)/64;
    size_t mid2=(end)/64;
    size_t head_len=mid1*64-begin;
    size_t tail_len=end-mid2*64;
    assert(mid1<=mid2);

    size_t res=__builtin_popcountll(get_word_part(data,begin,head_len));
    for(size_t i=mid1;i<mid2;i++){
      res+=__builtin_popcountll(data[i]);
    }
    res+=__builtin_popcountll(get_word_part(data,end-tail_len,tail_len));
    return res;
  }
}

int main(){
  int N,Q;
  scanf("%d %d ",&N,&Q);
  std::vector<uint64_t> buffer(N/64+1);
  for(int i=0;i<N;i++){
    char c;
    scanf("%c",&c);
    set_bit(buffer.data(),i,c-'0');
  }    
  for(int i=0;i<Q;i++){
    char op;
    scanf(" %c",&op);
    if(op=='C'||op=='A'||op=='O'||op=='X'){
      int X,Y,L;
      scanf("%d %d %d",&X,&Y,&L);
      X--,Y--;
      assert(0<=X&&X+L<=N);
      assert(0<=Y&&Y+L<=N);
      switch(op){
      case 'C': bitrange_assign(buffer.data(),X,buffer.data(),Y,L); break;
      case 'A': bitrange_assign<std::bit_and<uint64_t> >(buffer.data(),X,buffer.data(),Y,L); break;
      case 'O': bitrange_assign<std::bit_or<uint64_t> >(buffer.data(),X,buffer.data(),Y,L); break;
      case 'X': bitrange_assign<std::bit_xor<uint64_t> >(buffer.data(),X,buffer.data(),Y,L); break;
      }
    }else if(op=='N'){
      int X,L;
      scanf("%d %d",&X,&L);
      X--;
      assert(0<=X&&X+L<=N);
      printf("%d\n",(int)bitrange_popcount(buffer.data(),X,L));
    }
  }
  for(int i=0;i<N;i++){
    printf("%c",static_cast<char>('0'+get_bit(buffer.data(),i)));
  }
  printf("\n");
}
