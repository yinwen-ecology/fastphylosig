#ifndef FASTPHYLOSIG_ACE_EXPONENT_CACHE_H
#define FASTPHYLOSIG_ACE_EXPONENT_CACHE_H

#include <algorithm>
#include <array>
#include <cfenv>
#include <complex>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <utility>
#include <vector>

// Exponent values depend on the eigenvalue, immutable branch length and FP
// controls, not eigenvectors or tip states. Cache only exact bitwise keys.
class AceExponentCache {
  static constexpr std::size_t capacity=24;
  struct Key {
    std::array<std::uint64_t,2> bits{};
    explicit Key(const std::complex<double>& z = {}) {
      const double re=z.real(),im=z.imag();
      std::memcpy(&bits[0],&re,sizeof(re));std::memcpy(&bits[1],&im,sizeof(im));
    }
    bool operator==(const Key& b)const{return bits==b.bits;}
  };
  struct Hash {
    std::size_t operator()(const Key& key)const {
      return std::hash<std::uint64_t>{}(key.bits[0]) ^
        (std::hash<std::uint64_t>{}(key.bits[1])+0x9e3779b9u+(key.bits[0]<<6)+(key.bits[0]>>2));
    }
  };
  struct Slot {
    Key key;
    std::uint64_t stamp=0;
    std::vector<std::complex<double>> values;
  };
  std::vector<Slot> slots;
  std::unordered_map<Key,std::size_t,Hash> index;
  // Data pointers remain valid until the next prepare/configure. At most three
  // slots are pinned by the current stamp, so a miss cannot evict a selection.
  std::array<const std::complex<double>*,3> selected{};
  std::array<std::uint64_t,2> controls{};
  std::size_t edges=0,used=0;
  std::uint64_t clock=0;

  static bool fp_controls(std::array<std::uint64_t,2>& out) {
    out[0]=static_cast<std::uint64_t>(std::fegetround());
#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
    unsigned short cw;__asm__ __volatile__("fnstcw %0" : "=m"(cw));
    unsigned int mxcsr=0;
#if defined(__SSE__)
    __asm__ __volatile__("stmxcsr %0" : "=m"(mxcsr));
#endif
    // Sticky exception flags do not affect values; retain all control bits,
    // including x87 precision and SSE flush/denormal modes.
    out[1]=(static_cast<std::uint64_t>(cw)<<32)|(mxcsr&~0x3fu);
    return true;
#elif defined(__GNUC__) && defined(__aarch64__)
    __asm__ __volatile__("mrs %0, fpcr" : "=r"(out[1]));
    return true;
#else
    // An unknown control representation takes the uncached production path.
    return false;
#endif
  }
public:
  AceExponentCache()=default;
  AceExponentCache(const AceExponentCache&)=delete;
  AceExponentCache& operator=(const AceExponentCache&)=delete;
  void configure(const std::size_t n_edge) {
    selected.fill(nullptr);index.clear();slots.clear();
    used=0;clock=0;controls={};edges=n_edge;
    if(n_edge) {slots.resize(capacity);index.reserve(capacity);}
  }
  template<class Eigenvalues,class Lengths,class Evaluate>
  bool prepare(const Eigenvalues& lambda,const Lengths& lengths,Evaluate evaluate) {
    if(slots.empty())return false;
    std::array<std::uint64_t,2> current{};
    if(!fp_controls(current))return false;
    if(current!=controls||clock==UINT64_MAX) {
      index.clear();used=0;clock=0;controls=current;
      for(auto& slot:slots){slot.stamp=0;slot.values.clear();}
    }
    ++clock;selected.fill(nullptr);
    for(std::size_t j=0;j<3;++j) {
      const Key key(lambda[j]);auto found=index.find(key);
      std::size_t slot;
      if(found!=index.end())slot=found->second;
      else {
        slot=used<slots.size() ? used++ : static_cast<std::size_t>(
          std::min_element(slots.begin(),slots.end(),[](const Slot& a,const Slot& b){return a.stamp<b.stamp;})-slots.begin());
        // Reuse the evicted node without allocating in a full-cache miss.
        auto node=slots[slot].stamp ? index.extract(slots[slot].key) : decltype(index)::node_type{};
        // Publish only complete columns. A failed evaluation/allocation leaves
        // an unindexed slot that can safely be filled by a later retry.
        slots[slot].stamp=0;
        slots[slot].key=key;slots[slot].values.resize(edges);
        for(std::size_t row=0;row<edges;++row)slots[slot].values[row]=evaluate(lambda[j],lengths[row]);
        if(node.empty())index.emplace(key,slot);
        else {node.key()=key;node.mapped()=slot;index.insert(std::move(node));}
      }
      slots[slot].stamp=clock;selected[j]=slots[slot].values.data();
    }
    return true;
  }
  const std::complex<double>& value(const std::size_t j,const std::size_t edge)const {
    return selected[j][edge];
  }
  std::size_t bytes()const {
    std::size_t total=slots.capacity()*sizeof(Slot)+index.bucket_count()*sizeof(void*)+
      // Include node links/cached hashes; malloc metadata remains excluded.
      index.size()*(sizeof(std::pair<const Key,std::size_t>)+2*sizeof(void*));
    for(const auto& slot:slots)total+=slot.values.capacity()*sizeof(std::complex<double>);
    return total;
  }
};
#endif
