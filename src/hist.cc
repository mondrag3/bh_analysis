// #include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
// #include <fstream>
#include <string>
#include <vector>
// #include <unordered_map>
#include <utility>
// #include <stdexcept>
// #include <memory>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
// #include <TDirectory.h>

#include "selector.h"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

// istream range operator *******************************************
template<typename T>
struct range: public pair<T,T> {
  range(): pair<T,T>(0,0) { }
};
template<typename T>
istream& operator>> (istream& is, range<T>& r) {
  string str;
  is >> str;
  const size_t sep = str.find(':');
  if (sep==string::npos) {
    r.first = 0;
    stringstream(str) >> r.second;
  } else {
    stringstream(str.substr(0,sep)) >> r.first;
    stringstream(str.substr(sep+1)) >> r.second;
  }
  return is;
}

// NOTE: Cannot make istream operator for std::pair work with boost::po

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, sj_files, wt_files,
                 jet_algs, weights;
  string output_file, css_file;
  double pt_cut, eta_cut;
  range<Long64_t> num_events;
  bool quiet;

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
      ("output,o", po::value<string>(&output_file)->required(),
       "output root file with histograms")
      ("jet-alg", po::value<vector<string>>(&jet_algs)
       ->default_value({"AntiKt4"},"AntiKt4"),
       "jet algorithms from SJ file")
      ("weight,w", po::value<vector<string>>(&weights),
       "weight branch from weights file, e.g. Fac0.5Ht_Ren0.5Ht_PDFCT10_cent; "
       "if skipped, all weights from wt files are used")
      ("pt-cut", po::value<double>(&pt_cut)->default_value(30.),
       "jet pT cut")
      ("eta-cut", po::value<double>(&eta_cut)->default_value(4.4,"4.4"),
       "jet eta cut")
      ("style,s", po::value<string>(&css_file)
       ->default_value(CONFDIR"/Hj.css","Hj.css"),
       "CSS style file for histogram binning and formating")
      ("num-events,n", po::value< range<Long64_t> >(&num_events),
       "process only this many events, num or first:num")
      ("quiet,q", po::bool_switch(&quiet),
       "Do not print exception messages")
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

  // Find number of events to process
  if (num_events.second>0) {
    const Long64_t need_events = num_events.first + num_events.second;
    if (need_events>tree->GetEntries()) {
      cerr << "Fewer entries in BH chain (" << tree->GetEntries()
         << ") then requested (" << need_events << ')' << endl;
      exit(1);
    }
    if (need_events>sj_tree->GetEntries()) {
      cerr << "Fewer entries in SJ chain (" << sj_tree->GetEntries()
         << ") then requested (" << need_events << ')' << endl;
      exit(1);
    }
    if (need_events>wt_tree->GetEntries()) {
      cerr << "Fewer entries in weights chain (" << wt_tree->GetEntries()
         << ") then requested (" << need_events << ')' << endl;
      exit(1);
    }
  } else {
    num_events.second = tree->GetEntries();
    if (num_events.second!=sj_tree->GetEntries()) {
      cerr << num_events.second << " entries in BH chain, but "
           << sj_tree->GetEntries() << " entries in SJ chain" << endl;
      exit(1);
    }
    if (num_events.second!=wt_tree->GetEntries()) {
      cerr << num_events.second << " entries in BH chain, but "
           << wt_tree->GetEntries() << " entries in weights chain" << endl;
      exit(1);
    }
  }

  // Friend BlackHat tree with SpartyJet and Weight trees
  tree->AddFriend(sj_tree,"SJ");
  tree->AddFriend(wt_tree,"weights");

  // BlackHat tree branches
  BHEvent::event.SetTree(tree, BHEvent::kinematics);

  // SpartyJet tree branches
  if (jet_algs.size()) {
    cout << "Selected jet clustering algorithms:" << endl;
    for (auto& j : jet_algs) {
      cout << j << endl;
      SJClusterAlg::add(tree,j);
    }
  } else {
    cout << "Using all jet clustering algorithms:" << endl;
    const TObjArray *br = sj_tree->GetListOfBranches();
    for (Int_t i=0,n=br->GetEntries();i<n;++i) {
      auto j = br->At(i)->GetName();
      cout << j << endl;
      SJClusterAlg::add(tree,j);
    }
  }
  cout << endl;

  // Weights tree branches
  if (weights.size()) {
    cout << "Selected weights:" << endl;
    for (auto& w : weights) {
      cout << w << endl;
      weight::add(tree,w);
    }
  } else {
    cout << "Using all weights:" << endl;
    const TObjArray *br = wt_tree->GetListOfBranches();
    for (Int_t i=0,n=br->GetEntries();i<n;++i) {
      auto w = br->At(i)->GetName();
      cout << w << endl;
      weight::add(tree,w);
    }
  }
  cout << endl;

  // Read CSS file with histogram properties
  cout << "Histogram CSS file: " << css_file << endl;
  hist::css.reset( new kiwi::csshists(css_file) );
  cout << endl;

  // TODO: Check if CSS file exists

  // Open output file with histograms *******************************
  TFile* fout = new TFile(output_file.c_str(),"recreate");
  if (fout->IsZombie()) exit(1);
  else cout << "Output file: " << fout->GetName() << endl << endl;

  // Make directories ***********************************************
  for (auto& w : weight::all) {
    hist_wt::dirs[w.get()] = fout->mkdir(w->name.c_str());
  }

  for (auto& j : SJClusterAlg::all) {
    const auto dir = fout->mkdir(j->name.c_str());
    for (auto& w : weight::all) {
      hist_alg_wt::dirs[make_pair(j.get(),w.get())]
        = dir->mkdir(w->name.c_str());
    }
  }

  // Reading events from the input TChain ***************************
  cout << "Reading " << num_events.second << " entries";
  if (num_events.first>0) cout << " starting at " << num_events.first << endl;
  else cout << endl;
  timed_counter counter;

  selector *analysis = new selector(&counter);
  tree->Process(analysis,"",num_events.first,num_events.second);
  delete analysis;

  counter.prt(num_events.second);
  cout << endl;
  //cout << "Successfully processed events: " << numOK << endl;

  // Close files
  fout->Write();
  fout->Close();
  delete fout;
  delete tree;
  delete sj_tree;
  delete wt_tree;

  return 0;
}
