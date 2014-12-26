#include "selector_base.h"

#include <TTree.h>
#include <TDirectory.h>

using namespace std;

selector_base::selector_base(timed_counter *counter)
: TSelector(), fTree(nullptr), numOK(0), counter(counter)
{ }
selector_base::~selector_base() { }

void selector_base::Init(TTree *tree) {
  fTree = tree;
}

Bool_t selector_base::Process(Long64_t entry) {
  (*counter)(entry);
  fTree->GetEntry(entry);

  if (BHEvent::event.nparticle>BHMAXNP) {
    cerr << "More particles in the event then MAXNP" << endl
         << "Increase array length to " << BHEvent::event.nparticle << endl;
    exit(1);
  }

  return kTRUE;
}

// Histogram wrappers ***********************************************
unique_ptr<const kiwi::csshists> hist::css;
const SJClusterAlg* hist::alg_ptr;

hist_wt::hist_wt(const string& name) {
  TH1* hist = css->mkhist(name);
  hist->Sumw2(false); // in ROOT6 true seems to be the default
  for (auto& wt : weight::all) {
    const weight *w = wt.get();
    dirs[w]->cd();
    h[w] = static_cast<TH1*>( hist->Clone() );
  }
  delete hist;
}

void hist_wt::Fill(Double_t x) noexcept {
  for (auto& wt : weight::all) h[wt.get()]->Fill(x,wt->w);
}

unordered_map<const weight*,TDirectory*> hist_wt::dirs;

hist_alg_wt::hist_alg_wt(const string& name) {
  TH1* hist = css->mkhist(name);
  hist->Sumw2(false); // in ROOT6 true seems to be the default
  for (auto& alg : SJClusterAlg::all) {
    for (auto& wt : weight::all) {
      const auto k = make_pair(alg.get(),wt.get());
      dirs[k]->cd();
      h.emplace(k, static_cast<TH1*>( hist->Clone() ) );
    }
  }
  delete hist;
}

void hist_alg_wt::Fill(Double_t x) noexcept {
  key k(alg_ptr,nullptr);
  for (auto& wt : weight::all) {
    k.second = wt.get();
    h[k]->Fill(x,wt->w);
  }
}
void hist_alg_wt::FillOverflow() noexcept {
  key k(alg_ptr,nullptr);
  for (auto& wt : weight::all) {
    k.second = wt.get();
    TH1* _h = h[k];
    Int_t obin = _h->GetNbinsX()+1;
    _h->SetBinContent(obin,_h->GetBinContent(obin)+wt->w);
  }
}

unordered_map<hist_alg_wt::key, TDirectory*,
              pair_hash<hist_alg_wt::key>> hist_alg_wt::dirs;
