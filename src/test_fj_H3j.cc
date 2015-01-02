#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <tuple>
#include <array>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "BHEvent.h"
// #include "SJClusterAlg.h"
// #include "finder.h"
#include "timed_counter.h"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
using namespace fastjet;

// template<typename T> inline T sq(const T& x) { return x*x; }

int main(int argc, char** argv)
{
  if (argc!=3) {
    cout << "Usage: " << argv[0] << " bh.root wt.root" << endl;
    exit(0);
  }

  // Setup input files **********************************************
  TChain*    tree = new TChain("t3");
  TChain* wt_tree = new TChain("weights");

  // Add trees from all the files to the TChains
  if (!   tree->AddFile(argv[1],-1) ) exit(1);
  if (!wt_tree->AddFile(argv[2],-1) ) exit(1);

  // List weight branches
  tuple<string,Double_t,Double_t> _weight;
  vector<tuple<string,Float_t,Double_t>> weight;
  {
    get<0>(_weight)   = "weight";
    const TObjArray *wts = wt_tree->GetListOfBranches();
    const Int_t nwts = wts->GetEntries();
    weight.resize(nwts);
    for (Int_t i=0;i<nwts;++i)
      get<0>(weight[i]) = wts->At(i)->GetName();
  }

  // Friend BlackHat tree with SpartyJet and Weight trees
  tree->AddFriend(wt_tree,"weights");

  // BlackHat tree branches
  BHEvent event;
  event.SetTree(tree, BHEvent::kinematics);

  // Set weight branches address
  tree->SetBranchAddress( get<0>(_weight).c_str(), &get<1>(_weight) );
  for (auto& w : weight)
    tree->SetBranchAddress( get<0>(w).c_str(), &get<1>(w) );

  // Choose Clustering Algorithm
  const JetDefinition jet_def(antikt_algorithm, 0.4);
  cout << "\nClustering with " << jet_def.description() << endl << endl;

  // Reading events from the input TChain ***************************
  const Long64_t nent = tree->GetEntries();
  cout << "Reading " << nent << " entries" << endl;
  timed_counter counter;

  Long64_t selected = 0;

  for (Long64_t ent = 0; ent < nent; ++ent) {
    counter(ent);
    tree->GetEntry(ent);

    // test(ent)

    // Collect Final State particles
    vector<PseudoJet> particles;
    particles.reserve(event.nparticle);

    for (Int_t i=0; i<event.nparticle; ++i) {
      if (event.kf[i]==25) continue;
      particles.emplace_back(
        event.px[i],event.py[i],event.pz[i],event.E[i]
      );
    }

    // Perform clustering
    ClusterSequence cs(particles, jet_def);

    // Sort jets by Pt
    const vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    short njets = 0;
    for (auto& j : jets) {
      if (sqrt(j.kt2()) < 30.) continue;
      if (abs(j.eta()) > 4.4) continue;
      ++njets;
    }

    if (njets>=3) {
      get<2>(_weight) += get<1>(_weight);
      for (auto& w : weight) get<2>(w) += get<1>(w);
      ++selected;
    }
  }

  counter.prt(nent);
  cout << endl << endl;

  cout << "Accepted " << selected << " of " << nent << " events" << endl;
  cout << "Ntuple weight" << endl;
  cout << "σ = "
       << showpoint << setprecision(6) << get<2>(_weight)/nent
       << " pb" << endl;
  for (auto& w : weight) {
    cout << get<0>(w) << endl;
    cout << "σ = "
         << showpoint << setprecision(6) << get<2>(w)/nent
         << " pb" << endl;
  }

  delete tree;
  delete wt_tree;

  return 0;
}
