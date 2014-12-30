// #include <cmath>
#include <iostream>
#include <iomanip>
// #include <sstream>
// #include <fstream>
#include <string>
#include <vector>
// #include <unordered_map>
// #include <utility>
// #include <stdexcept>
// #include <memory>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
// #include <TDirectory.h>
// #include <TH1.h>
#include <TLorentzVector.h>

// #include "BHEvent.h"
#include "SJClusterAlg.h"
#include "finder.h"
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
  string weight_branch;

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
      ("weight,b", po::value<string>(&weight_branch)
       ->default_value("Fac0.5Ht_Ren0.5Ht_PDFCT10nlo_cent"),
       "weight branch")
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

  // Reading events from the input TChain ***************************
  const Long64_t nent = tree->GetEntries();
  cout << "Reading " << nent << " entries" << endl;
  timed_counter counter;

  Double_t total_weight = 0.;
  Long64_t selected = 0;

  for (Long64_t ent = 0; ent < nent; ++ent) {
    counter(ent);
    tree->GetEntry(ent);

    // number of jets; sort jets by Pt
    const size_t njets = SJClusterAlg::all.front()->jetsByPt(30.,4.4).size();

    if (njets>=3) {
      if (orig_weight) total_weight += _weight;
      else total_weight += weight;
      ++selected;
    }
  }

  counter.prt(nent);
  cout << endl << endl;

  cout << "σ = "
       << showpoint << setprecision(6) << total_weight/nent
       << " pb" << endl
       << "Selected: " << selected << endl;

  delete tree;
  delete sj_tree;
  delete wt_tree;

  return 0;

}
