#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <stdexcept>
#include <memory>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TLorentzVector.h>

#include <kiwi/csshists.h>

#include "BHEvent.h"
#include "SJClusterAlg.h"
#include "finder.h"
#include "timed_counter.h"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

template<typename T> inline T sq(const T& x) { return x*x; }

template <typename T>
struct pair_hash {
  size_t operator()(const T &x ) const {
    return hash<typename T::first_type >()(x.first ) ^
           hash<typename T::second_type>()(x.second);
  }
};

// Weights collector ************************************************
struct weight {
  string name;
  Float_t w;
  weight(TTree *tree, const string& name): name(name) {
    if ( tree->SetBranchAddress(name.c_str(), &w) == TTree::kMissingBranch )
      exit(1);
  }

  static vector<unique_ptr<const weight>> all;

  static void add(TTree* tree, const string& name) {
    all.emplace_back( new weight(tree,name) );
  }
};
vector<unique_ptr<const weight>> weight::all;

// TODO: place histograms in directories

// Histogram wrappers ***********************************************
struct hist {
  static unique_ptr<const kiwi::csshists> css;
  static const SJClusterAlg* alg_ptr;
  virtual void Fill(Double_t x) noexcept =0;
};
unique_ptr<const kiwi::csshists> hist::css;
const SJClusterAlg* hist::alg_ptr;

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
  unordered_map<const weight*,TH1*> h;
public:
  hist_wt(const string& name) {
    TH1* hist = css->mkhist(name);
    hist->Sumw2(false); // in ROOT6 true seems to be the default
    for (auto& wt : weight::all) {
      const weight *w = wt.get();
      dirs[w]->cd();
      h[w] = static_cast<TH1*>( hist->Clone() );
    }
    delete hist;
  }

  virtual void Fill(Double_t x) noexcept {
    for (auto& wt : weight::all) h[wt.get()]->Fill(x,wt->w);
  }

  static unordered_map<const weight*,TDirectory*> dirs;
};
unordered_map<const weight*,TDirectory*> hist_wt::dirs;

class hist_alg_wt: public hist {
  typedef pair<const SJClusterAlg*,const weight*> key;
  unordered_map<key,TH1*,pair_hash<key>> h;
public:
  hist_alg_wt(const string& name) {
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

  virtual void Fill(Double_t x) noexcept {
    key k(alg_ptr,nullptr);
    for (auto& wt : weight::all) {
      k.second = wt.get();
      h[k]->Fill(x,wt->w);
    }
  }
  virtual void FillOverflow() noexcept {
    key k(alg_ptr,nullptr);
    for (auto& wt : weight::all) {
      k.second = wt.get();
      TH1* _h = h[k];
      Int_t obin = _h->GetNbinsX()+1;
      _h->SetBinContent(obin,_h->GetBinContent(obin)+wt->w);
    }
  }

  static unordered_map<key,TDirectory*,pair_hash<key>> dirs;
};
unordered_map<hist_alg_wt::key, TDirectory*,
              pair_hash<hist_alg_wt::key>> hist_alg_wt::dirs;

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

inline Double_t tau(Double_t j_pt, Double_t j_mass, Double_t j_eta, Double_t H_eta) {
  return sqrt( sq(j_pt) + sq(j_mass) )/( 2.*cosh(j_eta - H_eta) );
}

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
      ("jet-alg,j", po::value<vector<string>>(&jet_algs)
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
  BHEvent event;
  event.SetTree(tree, BHEvent::kinematics);

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

  fout->cd();

  // Book histograms ************************************************
  TH1* h_N   = hist::css->mkhist("N");
  TH1* h_pid = hist::css->mkhist("pid");

  hist_wt h_xs("xs"),
          h_H_mass("H_mass"), h_H_pt("H_pT"), h_H_y("H_y")
  ;

  #define h_(name) h_##name(#name)

  /* NOTE:
   * excl = exactly the indicated number of jets, zero if no j in name
   * incl = that many or more jets
   * VBF = vector boson fusion cut
   */

  // Book Histograms
  hist_alg_wt
    h_(NJet_incl), h_(NJet_excl), h_(NJet_incl_50), h_(NJet_excl_50),

    h_(H_pT_excl),

    h_(jet1_pT), h_(jet1_pT_excl), h_(jet1_y), h_(jet1_tau),
    h_(jet2_pT), h_(jet2_y), h_(jet2_tau),
    h_(jet3_pT), h_(jet3_y), h_(jet3_tau),

    h_(jj_mass),
    h_(j_j_deltaphi), h_(j_j_deltaphi_excl), h_(j_j_deltaphi_VBF),
    h_(j_j_deltay),

    h_(H_j_pT), h_(H_j_pT_excl),
    h_(Hj_pT), h_(Hj_pT_excl),

    h_(H_jj_pT), h_(H_jj_pT_excl), h_(Hjj_mass),
    h_(H_jj_deltaphi), h_(H_jj_deltaphi_excl),
    h_(H_jj_deltay),
    h_(Hjj_pT), h_(Hjj_pT_excl),

    h_(H_jjj_pT),

    h_(loose), h_(tight),
    h_(jets_HT), h_(jets_tau_max), h_(jets_tau_sum)
  ;

  // Reading events from the input TChain ***************************
  Long64_t numOK = 0;
  cout << "Reading " << num_events.second << " entries";
  if (num_events.first>0) cout << " starting at " << num_events.first << endl;
  else cout << endl;
  num_events.second += num_events.first;
  timed_counter counter;

  for (Long64_t ent = num_events.first; ent < num_events.second; ++ent) {
    counter(ent);
    tree->GetEntry(ent);

    if (event.nparticle>BHMAXNP) {
      cerr << "More particles in the event then MAXNP" << endl
           << "Increase array length to " << event.nparticle << endl;
      exit(1);
    }

    // map particle pdg id to index number in the array
    finder<Int_t> pdg(event.kf,event.nparticle);
    size_t hi; // Higgs index

    try {
      hi = pdg(25,1); // find Higgs

      numOK++; // count number of good events
    } catch (exception& e) {
      if (!quiet) {
        cerr << "In event " << ent << ": ";
        cerr << e.what() << endl;

        for (Int_t i=0;i<event.nparticle;i++) cerr << event.kf[i] << ' ';
        cerr << endl;

      }
      continue; // skip to next event
    }

    h_N->Fill(0.5);

    // Higgs 4-vector
    const TLorentzVector higgs(event.px[hi],event.py[hi],event.pz[hi],event.E[hi]);

    const Double_t H_mass = higgs.M();        // Higgs Mass
    const Double_t H_pt   = higgs.Pt();       // Higgs Pt
    const Double_t H_eta  = higgs.Rapidity(); // Higgs Rapidity

    // Fill histograms ***********************************
    for (Int_t i=0;i<event.nparticle;i++) h_pid->Fill(event.kf[i]);

    h_xs.Fill(0.5);

    h_H_mass.Fill(H_mass);
    h_H_pt  .Fill(H_pt);
    h_H_y   .Fill(H_eta);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Loop over SpartyJet clustering algorithms
    for (auto& alg : SJClusterAlg::all) {
      hist::alg_ptr = alg.get(); // set static algorithm pointer for histograms

      // sort jets by Pt
      const vector<TLorentzVector> jets = alg->jetsByPt(pt_cut,eta_cut);
      const size_t njets = jets.size(); // number of jets

      int njets50 = 0;
      for (auto& j : jets) {
        if (j.Pt()>=50.) ++njets50;
        else break;
      }

      // Number of jets hists
      h_NJet_excl.Fill(njets);
      h_NJet_excl_50.Fill(njets50);
      for (unsigned char i=0;i<4;i++) {
        if (njets  >=i) h_NJet_incl   .Fill(i);
        if (njets50>=i) h_NJet_incl_50.Fill(i);
      }

      if (njets==0) { // njets == 0;

        h_H_pT_excl.Fill(H_pt);

      }
      else { // njets > 0;

        Double_t jets_HT = 0;
        static const Double_t jet_tau_cut=8;
        Double_t max_tj=0;
        Double_t sum_tj=0;

        for (auto& jet : jets) {
          const Double_t jet_pt  = jet.Pt();
          const Double_t jet_tau = tau(jet_pt,jet.M(),jet.Rapidity(),H_eta);

          jets_HT += jet_pt;

          if ( jet_tau > jet_tau_cut ) {
            sum_tj += jet_tau;
            if (jet_tau > max_tj) max_tj = jet_tau;
          }
        }
        h_jets_HT.Fill(jets_HT);
        h_jets_tau_max.Fill(max_tj);
        h_jets_tau_sum.Fill(sum_tj);

        const TLorentzVector& j1 = jets[0]; // First jet

        const Double_t j1_mass = j1.M();
        const Double_t j1_pt   = j1.Pt();
        const Double_t j1_eta  = j1.Rapidity();

        const Double_t Hj_pt = (higgs + j1).Pt();

        h_jet1_pT.Fill(j1_pt);
        h_jet1_y .Fill(j1_eta);
        h_Hj_pT  .Fill(Hj_pt);
        h_H_j_pT .Fill(H_pt);

        h_jet1_tau.Fill( tau(j1_pt,j1_mass,j1_eta,H_eta) );

        if (njets==1) { // njets == 1;

          h_Hj_pT_excl  .Fill(Hj_pt);
          h_H_j_pT_excl .Fill(H_pt);
          h_jet1_pT_excl.Fill(j1_pt);

        }
        else { // njets > 1;

          const TLorentzVector& j2 = jets[1]; // Second jet

          const Double_t j2_mass = j2.M();
          const Double_t j2_pt   = j2.Pt();
          const Double_t j2_eta  = j2.Rapidity();

          const TLorentzVector jj(j1+j2);
          const TLorentzVector Hjj(higgs + jj);

          const Double_t jj_mass        = jj.M();
          const Double_t deltaPhi_j_j   = j1.Phi() - j2.Phi();
          const Double_t deltaPhi_H_jj  = higgs.Phi() - jj.Phi();
          const Double_t Hjj_pt         = Hjj.Pt();
          const Double_t Hjj_mass       = Hjj.M();
          const Double_t deltay_j_j     = abs(j1_eta - j2_eta);
          const Double_t H_jj_deltay    = abs(H_eta-jj.Rapidity());

          h_j_j_deltaphi .Fill(deltaPhi_j_j);
          h_H_jj_deltaphi.Fill(deltaPhi_H_jj);
          h_Hjj_pT       .Fill(Hjj_pt);
          h_H_jj_pT      .Fill(H_pt);
          h_jet2_pT      .Fill(j2_pt);
          h_jet2_y       .Fill(j2_eta);
          h_jj_mass      .Fill(jj_mass);
          h_Hjj_mass     .Fill(Hjj_mass);
          h_j_j_deltay   .Fill(deltay_j_j);
          h_H_jj_deltay  .Fill(H_jj_deltay);

          if (deltay_j_j>2.8) {
            if (jj_mass>400) {
              h_j_j_deltaphi_VBF.Fill(deltaPhi_H_jj);
              h_loose.Fill(1);
              if (deltaPhi_H_jj>2.6) h_tight.Fill(1);
            }
          }

          h_jet2_tau.Fill( tau(j2_pt,j2_mass,j2_eta,H_eta) );

          if (njets==2) { // njets == 2;

            h_j_j_deltaphi_excl  .Fill(deltaPhi_j_j);
            h_H_jj_deltaphi_excl .Fill(deltaPhi_H_jj);
            h_Hjj_pT_excl        .Fill(Hjj_pt);
            h_H_jj_pT_excl       .Fill(H_pt);

          }
          else { // njets > 2;
            const TLorentzVector& j3 = jets[2]; // Second jet

            h_H_jjj_pT.Fill(H_pt);

            const Double_t j3_mass = j3.M();
            const Double_t j3_pt   = j3.Pt();
            const Double_t j3_eta  = j3.Rapidity();

            h_jet3_pT.Fill(j3_pt);
            h_jet3_y .Fill(j3_eta);

            h_jet3_tau.Fill( tau(j3_pt,j3_mass,j3_eta,H_eta) );

          } // END njets > 2;

        } // END njets > 1;

      } // END njets > 0;

    } // END Loop over SpartyJet clustering algorithms

  } // END of event loop

  counter.prt(num_events.second);
  cout << endl;
  cout << "Successfully processed events: " << numOK << endl;

  // Close files
  fout->Write();
  fout->Close();
  delete fout;
  delete tree;
  delete sj_tree;
  delete wt_tree;

  return 0;
}
