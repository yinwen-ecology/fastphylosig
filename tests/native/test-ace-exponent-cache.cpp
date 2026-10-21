// Standalone lifecycle/allocation audit. Run with libstdc++ debug iterators:
// g++ -std=c++17 -O2 -D_GLIBCXX_DEBUG -Isrc tests/native/test-ace-exponent-cache.cpp
#include "ace_exponent_cache.h"
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <stdexcept>

namespace {
// Guard both ends of every C++ allocation; track requested bytes (not malloc
// metadata). No R process, RNG, BLAS, or package source instrumentation needed.
struct alignas(std::max_align_t) Block {std::size_t size;std::uint64_t guard[2];};
constexpr std::uint64_t guard=0x4a19638fbce205d7ULL;
std::size_t live=0,allocations=0,blocks=0,peak=0;
int fail_after=-1;
void* allocate(std::size_t n) {
  if(fail_after==0)throw std::bad_alloc();
  if(fail_after>0)--fail_after;
  auto* b=static_cast<Block*>(std::malloc(sizeof(Block)+n+sizeof(guard)));
  if(!b)throw std::bad_alloc();
  b->size=n;b->guard[0]=b->guard[1]=guard;
  auto* p=b+1;std::memcpy(static_cast<char*>(static_cast<void*>(p))+n,&guard,sizeof(guard));
  live+=n;peak=std::max(peak,live);++allocations;++blocks;return p;
}
void deallocate(void* p) noexcept {
  if(!p)return;
  auto* b=static_cast<Block*>(p)-1;std::uint64_t end;
  std::memcpy(&end,static_cast<char*>(p)+b->size,sizeof(end));
  assert(b->guard[0]==guard&&b->guard[1]==guard&&end==guard);
  live-=b->size;--blocks;std::memset(p,0xdd,b->size);std::free(b);
}
using Z=std::complex<double>;
Z scalar(const Z& value,double length) {
  const Z z=value*Z(length,0.0);
  return z.imag()==0&&z.real()>=-700&&z.real()<=700 ? Z(std::exp(z.real()),z.imag()):std::exp(z);
}
bool exact(const Z& a,const Z& b) {return std::memcmp(&a,&b,sizeof(a))==0;}
void check(AceExponentCache& cache,const std::array<Z,3>& lambda,const std::vector<double>& lengths) {
  assert(cache.prepare(lambda,lengths,scalar));
  for(std::size_t j=0;j<3;++j)for(std::size_t i=0;i<lengths.size();++i)
    if(!exact(cache.value(j,i),scalar(lambda[j],lengths[i]))) {
      std::fprintf(stderr,"stale/incorrect exponent at column %zu, edge %zu\n",j,i);std::exit(2);
    }
}
}
void* operator new(std::size_t n) {return allocate(n);}
void* operator new[](std::size_t n) {return allocate(n);}
void operator delete(void* p) noexcept {deallocate(p);}
void operator delete[](void* p) noexcept {deallocate(p);}
void operator delete(void* p,std::size_t) noexcept {deallocate(p);}
void operator delete[](void* p,std::size_t) noexcept {deallocate(p);}

int main(int argc,char** argv) {
  const auto before_live=live,before_blocks=blocks;
  const bool reconfigure=argc>1&&std::strcmp(argv[1],"reconfigure")==0;
  const bool no_allocations=argc>1&&std::strcmp(argv[1],"no-allocations")==0;
  const bool accounting=argc>1&&std::strcmp(argv[1],"accounting")==0;
  {
    std::vector<double> lengths(64);
    for(std::size_t i=0;i<lengths.size();++i)lengths[i]=.01*(i+1);
    const auto external=live;
    AceExponentCache cache;cache.configure(lengths.size());
    std::array<Z,3> lambda={Z(-.2,0),Z(-.5,0),Z(0,0)};
    check(cache,lambda,lengths);
    if(reconfigure) {
      lengths.resize(128);
      for(std::size_t i=0;i<lengths.size();++i)lengths[i]=.02*(i+1);
      cache.configure(lengths.size());check(cache,lambda,lengths);
      cache.configure(0);
      assert(!cache.prepare(lambda,lengths,scalar));
    } else {
      // Allocation failures can occur before publication of a new column or
      // its hash entry. Keep the same cache and retry after each failure.
      for(int position=0;position<2;++position) {
        AceExponentCache retry;retry.configure(lengths.size());
        fail_after=position;bool threw=false;
        try {retry.prepare(lambda,lengths,scalar);}catch(const std::bad_alloc&){threw=true;}
        fail_after=-1;assert(threw);check(retry,lambda,lengths);
        auto changed=lambda;changed[0]=Z(-.3,0);check(retry,changed,lengths);
        check(retry,lambda,lengths);
      }
      // Repeated eigenvalues must not evict earlier selections in one call.
      lambda={Z(-.2,0),Z(-.2,0),Z(-.5,0)};check(cache,lambda,lengths);
      // Fail both a partially filled cache and a full-cache replacement. A
      // later fit must never read a partially evaluated column.
      for(int trial=0;trial<100;++trial) {
        lambda={Z(-.1-trial*.001,0),Z(-.8-trial*.002,.2),Z(-.8-trial*.002,-.2)};
        if(trial%7==0) {
          int remaining=trial%63+1;bool threw=false;
          try {cache.prepare(lambda,lengths,[&](const Z& z,double t){
            if(--remaining==0)throw std::runtime_error("injected evaluation failure");
            return scalar(z,t);
          });}catch(const std::runtime_error&){threw=true;}
          assert(threw);
        }
        check(cache,lambda,lengths);
      }
      const auto warm_allocations=allocations;
      for(int trial=100;trial<10100;++trial) {
        lambda={Z(-.1-trial*.001,0),Z(-.8-trial*.002,.2),Z(-.8-trial*.002,-.2)};
        check(cache,lambda,lengths);
      }
      const auto replacements=allocations-warm_allocations;
      std::printf("replacement allocations: %zu; reported cache bytes: %zu; allocated bytes: %zu\n",
        replacements,cache.bytes(),live-external);
      if(no_allocations)assert(replacements==0);
      if(accounting&&cache.bytes()<live-external)return 4;
      assert(cache.bytes()<=24*lengths.size()*sizeof(Z)+8192);
      const int previous=std::fegetround();
      for(int rounding:{FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO}) {
        assert(std::fesetround(rounding)==0);
        lambda={Z(-.0,+.0),Z(+.0,-.0),Z(-.0,-.0)};check(cache,lambda,lengths);
        lambda={Z(-.2,.3),Z(-.2,-.3),Z(-800,0)};check(cache,lambda,lengths);
      }
      assert(std::fesetround(previous)==0);
    }
  }
  assert(live==before_live&&blocks==before_blocks);
  std::printf("guards/lifetime/bitwise values OK; peak requested bytes: %zu\n",peak);
}
