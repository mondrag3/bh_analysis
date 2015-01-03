#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <tuple>
#include <array>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TLorentzVector.h>

#include "SJClusterAlg.h"
#include "timed_counter.h"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

// template<typename T> inline T sq(const T& x) { return x*x; }

int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, sj_files, wt_files;
  // string weight_branch;

  try {
    // General Options ------------------------------------
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "produce help message")
      ("bh", po::value< vector<string> >(&bh_files)->required(),
       "add input BlackHat root file")
      ("sj", po::value< vector<string> >(&sj_files)->required(),
       "add input SpartyJet root file")
      ("wt", po::value< vector<string> >(&wt_files)->required(),
       "add input weights root file")
      // ("weight,b", po::value<string>(&weight_branch)
      //  ->default_value("Fac0.5Ht_Ren0.5Ht_PDFCT10nlo_cent"),
      //  "weight branch")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (argc == 1 || vm.count("help")) {
      cout << desc << endl;
      return 0;
    }
    po::notify(vm);
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Setup input files **********************************************
  TChain*    tree = new TChain("t3");
  TChain* sj_tree = new TChain("SpartyJet_Tree");
  TChain* wt_tree = new TChain("weights");

  // Add trees from all the files to the TChains
  for (auto& f : bh_files) if (!   tree->AddFile(f.c_str(),-1) ) exit(1);
  for (auto& f : sj_files) if (!sj_tree->AddFile(f.c_str(),-1) ) exit(1);
  for (auto& f : wt_files) if (!wt_tree->AddFile(f.c_str(),-1) ) exit(1);

  // Friend BlackHat tree with SpartyJet and Weight trees
  tree->AddFriend(sj_tree,"SJ");
  tree->AddFriend(wt_tree,"weights");

  // BlackHat tree branches
  // BHEvent event;
  // event.SetTree(tree, BHEvent::cross_section);

  // SpartyJet tree branches
  SJClusterAlg::add(tree,"AntiKt4");

/*
  Float_t weight;
  Double_t _weight;
  bool orig_weight = false;
  if ( tree->GetBranch(weight_branch.c_str()) ) {
    if (weight_branch.compare("weight")) {
      tree->SetBranchAddress(weight_branch.c_str(), &weight);
    } else {
      tree->SetBranchAddress(weight_branch.c_str(), &_weight);
      orig_weight = true;
    }
  } else {
    cerr << "Unknown branch: " << weight_branch << endl;
    exit(1);
  }
*/

  tuple<string,Double_t,Double_t> _weight;
  array<tuple<string,Float_t,Double_t>,3> weight;
  get<0>(_weight)   = "weight";
  get<0>(weight[0]) = "Fac0.25Ht_Ren0.25Ht_PDFCT10nlo_cent";
  get<0>(weight[1]) = "Fac0.5Ht_Ren0.5Ht_PDFCT10nlo_cent";
  get<0>(weight[2]) = "Fac1Ht_Ren1Ht_PDFCT10nlo_cent";
  tree->SetBranchAddress( get<0>(_weight).c_str(),
                         &get<1>(_weight) );
  for (short i=0;i<3;++i)
    tree->SetBranchAddress( get<0>(weight[i]).c_str(),
                           &get<1>(weight[i]) );

  // Reading events from the input TChain ***************************
  const Long64_t nent = tree->GetEntries();
  cout << "Reading " << nent << " entries" << endl;
  timed_counter counter;

  Long64_t selected = 0;

  for (Long64_t ent = 0; ent < nent; ++ent) {
    counter(ent);
    tree->GetEntry(ent);

    // number of jets; sort jets by Pt
    const size_t njets = SJClusterAlg::all.front()->jetsByPt(30.,4.4).size();

    if (njets>=3) {
      get<2>(_weight) += get<1>(_weight);
      for (short i=0;i<3;++i)
        get<2>(weight[i]) += get<1>(weight[i]);
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
  delete sj_tree;
  delete wt_tree;

  return 0;
}
