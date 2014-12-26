#include <TSelector.h>

#include <string>
#include <unordered_map>
#include <utility>
#include <memory>

#include <TH1.h>
#include <TLorentzVector.h>

#include "BHEvent.h"
#include "SJClusterAlg.h"
#include "weight.h"

#include <kiwi/csshists.h>
#include "timed_counter.h"

class TDirectory;

class selector_base: public TSelector {
protected:
  TTree *fTree;

  mutable Long64_t numOK;
  timed_counter *counter;

public:
  selector_base(timed_counter *counter);
  virtual ~selector_base();

  virtual void Init(TTree *tree);
  virtual Bool_t Process(Long64_t entry);
};

// Histogram wrappers ***********************************************
struct hist {
  static std::unique_ptr<const kiwi::csshists> css;
  static const SJClusterAlg* alg_ptr;
  virtual void Fill(Double_t x) noexcept =0;
};

/*
class hist_alg: public hist {
  unordered_map<const SJClusterAlg*,TH1*> h;
public:
  hist_alg(const string& name) {
    TH1* hist = css->mkhist(name);
    for (auto& alg : SJClusterAlg::all)
      h[alg.get()] = static_cast<TH1*>(
        hist->Clone( (alg->name+"__"+name).c_str() )
      );
    delete hist;
  }
  virtual void Fill(Double_t x) noexcept { h[alg_ptr]->Fill(x); }
};
*/

class hist_wt: public hist {
  std::unordered_map<const weight*,TH1*> h;

public:
  hist_wt(const std::string& name);
  virtual void Fill(Double_t x) noexcept;

  static std::unordered_map<const weight*,TDirectory*> dirs;
};

template <typename T>
struct pair_hash {
  size_t operator()(const T &x ) const {
    return std::hash<typename T::first_type >()(x.first ) ^
           std::hash<typename T::second_type>()(x.second);
  }
};

class hist_alg_wt: public hist {
  typedef std::pair<const SJClusterAlg*,const weight*> key;
  std::unordered_map<key,TH1*,pair_hash<key>> h;

public:
  hist_alg_wt(const std::string& name);
  virtual void Fill(Double_t x) noexcept;
  virtual void FillOverflow() noexcept;

  static std::unordered_map<key,TDirectory*,pair_hash<key>> dirs;
};
