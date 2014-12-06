#include <cmath>
#include <ctime>
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
#include <TH1.h>
#include <TLorentzVector.h>

#include <kiwi/csshists.h>

#include "BHEvent.h"
#include "SJClusterAlg.h"
#include "finder.h"

#define debug_var(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

template<typename T> inline T sq(const T& x) { return x*x; }

// Weights collector ************************************************
struct weight {
  string name;
  Float_t w;
  weight(TTree *tree, const string& name): name(name)
  { tree->SetBranchAddress(name.c_str(), &w); }

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
    for (auto& wt : weight::all)
      h[wt.get()] = static_cast<TH1*>(
        hist->Clone( (wt->name+"__"+name).c_str() )
      );
    delete hist;
  }

  virtual void Fill(Double_t x) noexcept {
    for (auto& wt : weight::all) h[wt.get()]->Fill(x,wt->w);
  }
};

class hist_alg_wt: public hist {
  typedef pair<const SJClusterAlg*,const weight*> key;
  struct key_hash {
    size_t operator()(const key& k) const {
      return hash<const SJClusterAlg*>()(k.first) ^ hash<const weight*>()(k.second);
    }
  };
  unordered_map<key,TH1*,key_hash> h;
public:
  hist_alg_wt(const string& name) {
    TH1* hist = css->mkhist(name);
    for (auto& alg : SJClusterAlg::all) {
      for (auto& wt : weight::all) {
        h.emplace( make_pair(alg.get(),wt.get()), static_cast<TH1*>(
          hist->Clone( (alg->name+'_'+wt->name+"__"+name).c_str() )
        ) );
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
};

namespace std {
  template <>
  struct hash<pair<const SJClusterAlg*,const weight*>>{
    size_t operator()(const pair<const SJClusterAlg*,const weight*> &x ) const {
      return hash<const SJClusterAlg*>()(x.first) ^ hash<const weight*>()(x.second);
    }
  };
}

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

// INFO: Cannot make istream operator for std::pair work with boost::po

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, sj_files, wt_files,
                 jet_algs, weights;
  string output_file, css_file;
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
      ("jet-alg,j", po::value<vector<string>>(&jet_algs)->default_value(vector<string>(),"AntiKt4"),
       "jet algorithms from SJ file") // TODO: read all if not provided
      ("weight,w", po::value<vector<string>>(&weights)->required(),
       "weight from weights file, e.g. Fac0.5Ht_Ren0.5Ht_PDFCT10_cent") // TODO: read all if not provided
      ("style,s", po::value<string>(&css_file)->required(),
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
    cerr << "Error: " <<  e.what() << endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Setup input files **********************************************
  TChain*    tree = new TChain("t3");
  TChain* sj_tree = new TChain("SpartyJet_Tree");
  TChain* wt_tree = new TChain("weights");

  // Add trees from all the files to the TChains
  for (auto& f : bh_files)    tree->Add(f.c_str());
  for (auto& f : sj_files) sj_tree->Add(f.c_str());
  for (auto& f : wt_files) wt_tree->Add(f.c_str());

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

  // Friend SpartyJet Tree with the parton level Tree
  tree->AddFriend(sj_tree,"SJ");
  tree->AddFriend(wt_tree,"weights");

  // BlackHat tree branches
  BHEvent event;
  event.SetTree(tree, BHEvent::kinematics);

  // SpartyJet tree branches
  cout << "Selected jet clustering algorithms:" << endl;
  for (auto& j : jet_algs) {
    cout << j << endl;
    SJClusterAlg::add(tree,j);
  }
  cout << endl;

  // Weights tree branches
  cout << "Selected weights:" << endl;
  for (auto& w : weights) {
    cout << w << endl;
    weight::add(tree,w);
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

  // Book histograms ************************************************
  TH1* h_pid = hist::css->mkhist("pid");

  hist_wt h_xs("xs"),
          h_H_mass("H_mass"), h_H_pt("H_pt"), h_H_y("H_y")
  ;

  #define h_(name) h_##name(#name)

  // Histogram maps
  hist_alg_wt
    h_(NJet_incl), h_(NJet_excl), h_(NJet_incl_50), h_(NJet_excl_50),
    h_(H_pT_excl), h_(H_pT_fine_excl),
    h_(jet1_pT), h_(jet1_pT_fine), h_(jet1_pT_excl), h_(jet1_pT_excl_fine),
    h_(jet1_y), h_(jet1_y_fine), h_(tau_jet1),
    h_(jet2_pT), h_(jet2_pT_fine), h_(jet2_y), h_(jet2_y_fine), h_(tau_jet2),
    h_(jet3_pT), h_(jet3_pT_fine), h_(jet3_y), h_(jet3_y_fine), h_(tau_jet3),
    h_(H_j_pT), h_(H_j_pT_fine), h_(H_j_pT_excl), h_(H_j_pT_fine_excl),
    h_(Hj_pT), h_(Hj_pT_fine), h_(Hj_pT_excl), h_(Hj_pT_fine_excl),
    h_(H_jj_pT), h_(H_jj_pT_fine), h_(H_jj_pT_excl), h_(H_jj_pT_fine_excl),
    h_(Hjj_pT), h_(Hjj_pT_fine), h_(Hjj_pT_excl), h_(Hjj_pT_fine_excl),
    h_(dijet_mass), h_(dijet_mass_fine), h_(H_dijet_mass),
    h_(deltaphi_jj), h_(deltaphi_jj_fine), h_(deltaphi_jj_excl),
    h_(deltaphi_Hjj), h_(deltaphi_Hjj_excl), h_(deltaphi_jj_VBF),
    h_(deltay_jj), h_(deltay_jj_fine), h_(deltay_H_jj), h_(deltay_H_jj_fine),
    h_(loose), h_(tight),
    h_(HT_jets_hist),
    h_(tau_jet_max), h_(sum_tau_jet)
  ;

  // Reading events from the input TChain ***************************
  Long64_t numOK = 0;
  cout << "Preparing to read " << num_events.second
       << " entries starting at " << num_events.first << endl;
  time_t last_time = time(0), cur_time;
  unsigned seconds = 0;

  for (Long64_t ent = num_events.first; ent < num_events.second; ++ent) {

    // timed counter
    if ( difftime( (cur_time=time(0)), last_time ) > 1 ) {
      cout << setw(10) << ent << setw(9) << seconds << 's';
      cout.flush();
      for (char i=0;i<20;i++) cout << '\b';
      last_time = cur_time;
      ++seconds;
    }

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
    h_H_y   .Fill(abs(H_eta));

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Loop over SpartyJet clustering algorithms
    for (auto& alg : SJClusterAlg::all) {
      hist::alg_ptr = alg.get(); // set static algorithm pointer for histograms

      // sort jets by Pt
      const vector<TLorentzVector> jets = alg->jetsByPt(30.,5.);
      const size_t njets = jets.size(); // number of jets

      h_NJet_excl.Fill(njets);
      for (unsigned char i=0;i<4;i++)
        if (njets>=i) h_NJet_incl.Fill(i);

      if (njets==0) { // njets == 0;

        h_H_pT_excl     .Fill(H_pt);
        h_H_pT_fine_excl.Fill(H_pt);
        // 6/2 added fill for jet1_pT for 0-30 GeV bin, i.e. no jets
        h_jet1_pT       .Fill(10);
        // 6/2 added fill for overflow bin for 0 jets in event
        h_deltay_jj     .FillOverflow();
        // 6/2 added fill for overflow bin for 0 jets in event
        h_Hjj_pT        .FillOverflow();
        //6/2 added fill for overflow bin for 0 jets in event
        h_jet2_y        .FillOverflow();
        //6/2 added fill for overflow bin for 0 jets in event
        h_jet2_pT       .FillOverflow();

        h_NJet_incl_50.Fill(0);
        h_NJet_excl_50.Fill(0);

      }
      else { // njets > 0;

        const TLorentzVector& j1 = jets[0]; // First jet

        const Double_t j1_mass = j1.M();
        const Double_t j1_pt   = j1.Pt();
        const Double_t j1_eta  = j1.Rapidity();

        const Double_t Hj_pt = (higgs + j1).Pt();

        h_jet1_pT     .Fill(j1_pt);
        h_jet1_pT_fine.Fill(j1_pt);
        h_jet1_y      .Fill(j1_eta);
        h_jet1_y_fine .Fill(j1_eta);
        h_Hj_pT       .Fill(Hj_pt);
        h_Hj_pT_fine  .Fill(Hj_pt);
        h_H_j_pT      .Fill(H_pt);
        h_H_j_pT_fine .Fill(H_pt);

        h_tau_jet1.Fill(
          sqrt( sq(j1_pt) + sq(j1_mass) )/( 2.0*cosh(j1_eta - H_eta) )
        );

        if (j1_pt>50.) {
          h_NJet_incl_50.Fill(1);
        } else {
          h_NJet_incl_50.Fill(0);
          h_NJet_excl_50.Fill(0);
        }

        if (njets==1) { // njets == 1;

          h_Hj_pT_excl        .Fill(Hj_pt);
          h_Hj_pT_fine_excl   .Fill(Hj_pt);
          h_H_j_pT_excl       .Fill(H_pt);
          h_H_j_pT_fine_excl  .Fill(H_pt);
          h_jet1_pT_excl      .Fill(j1_pt);
          h_jet1_pT_excl_fine .Fill(j1_pt);
          // 6/2 added fill for j2_pT for 0-30 GeV bins, i.e. no 2nd jet
          h_jet2_pT           .Fill(10);
          // 6/2 added fill for overflow bin for 1 jet in event
          h_deltay_jj         .FillOverflow();
          // 6/2 added fill for overflow bin for 1 jet in event
          h_Hjj_pT            .FillOverflow();
          //6/2 added fill for overflow bin for 1 jet in event
          h_jet2_y            .FillOverflow();

          if(j1_pt>50.) h_NJet_excl_50.Fill(1);

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
          const Double_t deltay_H_jj    = abs(H_eta-jj.Rapidity());

          h_deltaphi_jj     .Fill(deltaPhi_j_j);
          h_deltaphi_jj_fine.Fill(deltaPhi_j_j);
          h_deltaphi_Hjj    .Fill(deltaPhi_H_jj);
          h_Hjj_pT          .Fill(Hjj_pt);
          h_Hjj_pT_fine     .Fill(Hjj_pt);
          h_H_jj_pT         .Fill(H_pt);
          h_H_jj_pT_fine    .Fill(H_pt);
          h_jet2_pT         .Fill(j2_pt);
          h_jet2_pT_fine    .Fill(j2_pt);
          h_jet2_y          .Fill(abs(j2_eta));
          h_jet2_y_fine     .Fill(abs(j2_eta));
          h_dijet_mass      .Fill(jj_mass);
          h_dijet_mass_fine .Fill(jj_mass);
          h_H_dijet_mass    .Fill(Hjj_mass);
          h_deltay_jj       .Fill(deltay_j_j);
          h_deltay_jj_fine  .Fill(deltay_j_j);
          h_deltay_H_jj     .Fill(deltay_H_jj);
          h_deltay_H_jj_fine.Fill(deltay_H_jj);

          if (deltay_j_j>2.8) {
            if (jj_mass>400) {
              // 6/13 realized that delta-phi cut should be between Higgs and dijet system
              h_deltaphi_jj_VBF.Fill(deltaPhi_H_jj);
              // 6/2 added loose and tight histograms to tally cross section as cross-checks
              h_loose.Fill(1);
              if (deltaPhi_H_jj>2.6) h_tight.Fill(1);
            }
          }

          h_tau_jet2.Fill(
            sqrt( sq(j2_pt) + sq(j2_mass) )/( 2.0*cosh(j2_eta - H_eta) )
          );

          if (j2_pt>50.) h_NJet_incl_50.Fill(2);

          if (njets==2) { // njets == 2;

            h_deltaphi_jj_excl  .Fill(deltaPhi_j_j);
            h_deltaphi_Hjj_excl .Fill(deltaPhi_H_jj);
            h_Hjj_pT_excl       .Fill(Hjj_pt);
            h_Hjj_pT_fine_excl  .Fill(Hjj_pt);
            h_H_jj_pT_excl      .Fill(H_pt);
            h_H_jj_pT_fine_excl .Fill(H_pt);
            // 6/2 added fill for j3_pT 0-30 GeV, i.e. no jet3
            h_jet3_pT           .Fill(10);

            if (j2_pt>50.) h_NJet_excl_50.Fill(2);

          }
          else { // njets > 2;
            const TLorentzVector& j3 = jets[2]; // Second jet

            const Double_t j3_mass = j3.M();
            const Double_t j3_pt   = j3.Pt();
            const Double_t j3_eta  = j3.Rapidity();

            h_jet3_pT     .Fill(j3_pt);
            h_jet3_pT_fine.Fill(j3_pt);
            h_jet3_y      .Fill(j3_eta);
            h_jet3_y_fine .Fill(j3_eta);

            h_tau_jet3.Fill(
              sqrt( sq(j3_pt) + sq(j3_mass) )/( 2.*cosh(j3_eta - H_eta) )
            );

            if (j3_pt>50.) h_NJet_incl_50.Fill(3);

            if (njets==3) {

              if (j3_pt>50.) h_NJet_excl_50.Fill(3);

            }

          } // END njets > 2;

        } // END njets > 1;

      } // END njets > 0;

      Double_t HT_jets = 0;

      static const Double_t tau_jet_cut=8;
      Double_t max_tj=0;
      Double_t sum_tj=0;

      for(size_t i=0;i<njets;i++) {
        const TLorentzVector& jet = jets[i];

        const Double_t jet_pt   = jet.Pt();
        const Double_t jet_mass = jet.M();
        const Double_t jet_eta  = jet.Rapidity();

        HT_jets += jet_pt;

        Double_t tauJet =
          sqrt( sq(jet_pt) + sq(jet_mass) )/( 2.*cosh(jet_eta - H_eta) );

        if ( tauJet > tau_jet_cut ) {
          sum_tj += tauJet;
          if (tauJet > max_tj) max_tj = tauJet;
        }

      }
      h_HT_jets_hist.Fill(HT_jets);
      h_tau_jet_max.Fill(max_tj);
      h_sum_tau_jet.Fill(sum_tj);

    } // END Loop over SpartyJet clustering algorithms

  } // END of event loop

  cout << setw(10) << num_events.second-1 << setw(9) << seconds <<'s' << endl;

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
