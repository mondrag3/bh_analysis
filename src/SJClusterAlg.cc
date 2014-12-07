#include "SJClusterAlg.h"

#include <cmath>
#include <cstdlib>

#include "TTree.h"

using namespace std;

#define branch(var) \
  if ( tree->SetBranchAddress((name+'_'+#var).c_str(), &var) \
       == TTree::kMissingBranch ) exit(1);

SJClusterAlg::SJClusterAlg(TTree* tree, const string& name)
: N(0), eta(0), phi(0), e(0), mass(0), pt(0), numC(0), ind(0),
  name(name)
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

vector<unique_ptr<const SJClusterAlg>> SJClusterAlg::all;

void SJClusterAlg::add(TTree* tree, const string& name) {
  all.emplace_back( new SJClusterAlg(tree,name) );
}

class sortByPt { // class used for jet sorting
  const SJClusterAlg *p;
  public:
  sortByPt(const SJClusterAlg* p): p(p) { }
  bool operator() (Int_t i, Int_t j) const {
    return ( p->pt->at(i) > p->pt->at(j) ); // decending order
  }
};

vector<TLorentzVector> SJClusterAlg::jetsByPt(double pt_cut, double eta_cut) const {
  Int_t j[N]; // order of jets
  for (Int_t i=0;i<N;++i) j[i] = i;
  sort(j,j+N,sortByPt(this));

  vector<TLorentzVector> jets;
  jets.reserve(N);
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
