#include "SJClusterAlg.h"

#include <cmath>

#include "TTree.h"

using namespace std;

#define branch(var) \
  tree->SetBranchAddress((alg+'_'+#var).c_str(), &var);

SJClusterAlg::SJClusterAlg(TTree* tree, const char* algorithm)
: N(0), eta(0), phi(0), e(0), mass(0), pt(0), numC(0), ind(0),
  alg(algorithm)
{
  branch(N)
  branch(eta)
  branch(phi)
  //branch(e)
  branch(mass)
  branch(pt)
  //branch(ind)
  //branch(numC)
}
SJClusterAlg::~SJClusterAlg() { }

class SJClusterAlg::sortByPt { // class used for jet sorting
  const SJClusterAlg * const p;
  public:
  sortByPt(const SJClusterAlg* p): p(p) { }
  bool operator() (Int_t i, Int_t j) {
    return ( p->pt->at(i) > p->pt->at(j) ); // decending order
  }
};

vector<const SJClusterAlg*> SJClusterAlg::AllAlgs;

const SJClusterAlg* SJClusterAlg::AddAlg(TTree* tree, const char* algorithm) {
  AllAlgs.push_back( new SJClusterAlg(tree,algorithm) );
  return AllAlgs.back();
}
void SJClusterAlg::clean() {
  for (SJAlgIter it = begin(); it!=end(); ++it)
    delete *it;
}

const string& SJClusterAlg::name() const { return alg; }

vector<TLorentzVector> SJClusterAlg::jetsByPt(double pt_cut, double eta_cut) const {
  Int_t j[N]; // order of jets
  for (Int_t i=0;i<N;++i) j[i] = i;
  sort(j,j+N,sortByPt(this));

  vector<TLorentzVector> jets;
  for (Int_t i=0;i<N;++i) {
    const double _pt  = pt->at(j[i]);
    const double _eta = eta->at(j[i]);

    // cuts
    if (_pt  < pt_cut ) continue;
    if ( abs(_eta) > eta_cut) continue;

    jets.push_back(TLorentzVector());
    jets.back().SetPtEtaPhiM(_pt,_eta,phi->at(j[i]),mass->at(j[i]));
  }

  return jets;
}

SJAlgIter SJClusterAlg::begin() { return AllAlgs.begin(); }
SJAlgIter SJClusterAlg::  end() { return AllAlgs.end();   }
