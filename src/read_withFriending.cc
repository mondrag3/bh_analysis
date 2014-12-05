#include <cmath>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <utility>
#include <stdexcept>
#include <ctime>
#include <tuple>
#include <memory>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TLorentzVector.h>

#include <flock/csshists.h>

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
    all.emplace_back( new weight(name,tree) );
  }
};

// TODO: place histograms in directories

// Histogram wrappers ***********************************************
struct hist {
  static unique_ptr<const flock::csshists> css;
  static const SJClusterAlg* alg_ptr;
  virtual Fill(Double_t x) =0;
};
unique_ptr<const flock::csshists> hist::css;
const SJClusterAlg* hist::alg_ptr;

class hist_alg: public hist {
  unordered_map<const SJClusterAlg*,TH1*> h;
public:
  hist_alg(const string& name) {
    TH1* hist = css->mkhist(name);
    for (auto alg : SJClusterAlg::all)
      h[alg.get()] = hist->Clone((alg->name+"__"+name).c_str());
    delete hist;
  }

  virtual void Fill(Double_t x) { h[alg_ptr]->Fill(x); }
};

class hist_wt: public hist {
  unordered_map<const weight*,TH1*> h;
public:
  hist_alg_wt(const string& name) {
    TH1* hist = css->mkhist(name);
    for (auto wt : weight::all)
      h[wt.get()] = hist->Clone((wt->name+"__"+name).c_str());
    delete hist;
  }

  virtual void Fill(Double_t x) {
    for (auto wt : weight::all) h[wt.get()]->Fill(x,wt->w);
  }
};

class hist_alg_wt: public hist {
  unordered_map<tuple<const SJClusterAlg*,const weight*>,TH1*> h;
public:
  hist_alg_wt(const string& name) {
    TH1* hist = css->mkhist(name);
    for (auto alg : SJClusterAlg::all) {
      for (auto wt : weight::all) {
        h[make_tuple(alg.get(),wt.get())]
          = hist->Clone((alg->name+'_'+wt->name+"__"+name).c_str());
      }
    }
    delete hist;
  }

  virtual void Fill(Double_t x) {
    for (auto wt : weight::all)
      h[tie(alg_ptr,wt.get())]->Fill(x,wt->w);
  }
};

// istream pair operator ********************************************
template<typename T>
istream& operator >> (istream& is, pair<T,T>& p)
{
  string str;
  is >> str;
  size_t sep = str.find(':');
  if (sep==string::npos) {
    p.first = 0;
    stringstream ss(str);
    ss >> p.second;
  } else {
    stringstream ss(str.substr(0,sep));
    ss >> p.first;
    ss.str(str.substr(sep+1));
    ss >> p.second;
  }
  return is;
}

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, sj_files, wt_files,
                 jet_algs, weights;
  string output_file, css_file;
  pair<Long64_t,Long64_t> num_events {0,0};
  bool quiet;

  try {
    // General Options ------------------------------------
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "produce help message")
      ("bh", po::value< vector<string> >(&event_files)->required(),
       "add input BlackHat root file")
      ("sj", po::value< vector<string> >(&sj_files)->required(),
       "add input SpartyJet root file")
      ("wt", po::value< vector<string> >(&weight_files)->required(),
       "add input weights root file")
      ("output,o", po::value<string>(&output_file)->required(),
       "output root file with histograms")
      ("jet-alg,j", po::value<vector<string>>(&jet_algs)->required(),
       "jet algorithms from SJ file, e.g. AntiKt4") // TODO: read all if not provided
      ("weight,w", po::value<vector<string>>(&weights)->required(),
       "weight from weights file, e.g. Fac0.5Ht_Ren0.5Ht_PDFCT10_cent") // TODO: read all if not provided
      ("style,s", po::value<string>(&css_file)->required(),
       "CSS style file for histogram binning and formating")
      ("num-events,n", po::value<pair<Long64_t,Long64_t>>(&num_events),
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
    return 1;
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
           << treeSJ->GetEntries() << " entries in SJ chain" << endl;
      exit(1);
    }
    if (num_events.second!=wt_tree->GetEntries()) {
      cerr << num_events.second << " entries in BH chain, but "
           << treeW->GetEntries() << " entries in weights chain" << endl;
      exit(1);
    }
  }

  // Friend SpartyJet Tree with the parton level Tree
  tree->AddFriend(sj_tree,"SJ");
  tree->AddFriend(wt_tree,"weights");

  // BlackHat tree branches
  BHEvent event(tree, BHEvent::kinematics);

  // SpartyJet tree branches
  for (auto& j : jet_algs) SJClusterAlg::add(tree,j);

  // Weights tree branches
  for (auto& w : weights) weight::add(tree,w);

  // Read CSS file with histogram properties
  hist::css.reset( new csshists(css_file) );

  // Open output file with histograms *******************************
  TFile* fout = new TFile(output_file.c_str(),"recreate");

  // Book histograms ************************************************
  hist h_pid("pid",40,-10,30,"Particle id");

  hist h_xs("xs",1,0.,1.,"XS");
  hist h_H_mass("H_mass",100,0.,200.,"Higgs mass");
  hist h_H_pt("H_pt",H_pT_bins,"Higgs pT");
  hist h_H_pt_fine("H_pt_fine",H_pT_jj_fine_bins,"Higgs pT");
  hist h_H_y("H_y",H_y_bins,"Higgs rapidity");
  hist h_H_y_fine("H_y_fine",H_y_fine_bins,"Higgs rapidity");

  // Histogram maps
  hist_map
    h_NJet_incl, h_NJet_excl, h_NJet_incl_50, h_NJet_excl_50,
    h_H_pT_excl, h_H_pT_fine_excl,
    h_jet1_pT, h_jet1_pT_fine, h_jet1_pT_excl, h_jet1_pT_excl_fine,
    h_jet1_y, h_jet1_y_fine, h_tau_jet1,
    h_jet2_pT, h_jet2_pT_fine, h_jet2_y, h_jet2_y_fine, h_tau_jet2,
    h_jet3_pT, h_jet3_pT_fine, h_jet3_y, h_jet3_y_fine, h_tau_jet3,
    h_H_j_pT, h_H_j_pT_fine, h_H_j_pT_excl, h_H_j_pT_fine_excl,
    h_Hj_pT, h_Hj_pT_fine, h_Hj_pT_excl, h_Hj_pT_fine_excl,
    h_H_jj_pT, h_H_jj_pT_fine, h_H_jj_pT_excl, h_H_jj_pT_fine_excl,
    h_Hjj_pT, h_Hjj_pT_fine, h_Hjj_pT_excl, h_Hjj_pT_fine_excl,
    h_dijet_mass, h_dijet_mass_fine, h_H_dijet_mass,
    h_deltaphi_jj, h_deltaphi_jj_fine, h_deltaphi_jj_excl,
    h_deltaphi_Hjj, h_deltaphi_Hjj_excl, h_deltaphi_jj_VBF,
    h_deltay_jj, h_deltay_jj_fine, h_deltay_H_jj, h_deltay_H_jj_fine,
    h_loose, h_tight,
    h_HT_jets_hist,
    h_tau_jet_max, h_sum_tau_jet
  ;

  // Loop over SpartyJet clustering algorithms
  for ( hist_map::it=SJClusterAlg::begin();
        hist_map::it!=SJClusterAlg::end(); ++hist_map::it) {

    // Create directory in the file
    TDirectory* const dir = fout->mkdir((*hist_map::it)->name().c_str());
    dir->cd();

    h_NJet_incl.book("NJet_incl",5,-0.5,4.5,"Inclusive number of jets");
    h_NJet_excl.book("NJet_excl",5,-0.5,4.5,"Exclusive number of jets");
    h_NJet_incl_50.book("NJet_incl_50",5,-0.5,4.5);
    h_NJet_excl_50.book("NJet_excl_50",5,-0.5,4.5);

    h_H_pT_excl.book("H_pT_excl",H_pT_0j_bins);
    h_H_pT_fine_excl.book("H_pT_fine_excl",H_pT_jj_fine_bins);

    h_jet1_pT.book("jet1_pT",jet1_pT_bins);
    h_jet1_pT_fine.book("jet1_pT_fine",H_pT_jj_fine_bins);
    h_jet1_pT_excl.book("jet1_pT_excl",jet1_pT_bins);
    h_jet1_pT_excl_fine.book("jet1_pT_excl_fine",H_pT_jj_fine_bins);
    h_jet1_y.book("jet1_y",jet1_y_bins);
    h_jet1_y_fine.book("jet1_y_fine",jet1_y_fine_bins);

    h_jet2_pT.book("jet2_pT",jet2_pT_bins);
    h_jet2_pT_fine.book("jet2_pT_fine",H_pT_jj_fine_bins);
    h_jet2_y.book("jet2_y",jet1_y_bins);
    h_jet2_y_fine.book("jet2_y_fine",jet1_y_fine_bins);

    h_jet3_pT.book("jet3_pT",jet3_pT_bins);
    h_jet3_pT_fine.book("jet3_pT_fine",H_pT_jj_fine_bins);
    h_jet3_y.book("jet3_y",jet1_y_bins);
    h_jet3_y_fine.book("jet3_y_fine",jet1_y_fine_bins);

    h_H_j_pT.book("H_j_pT",H_pT_j_bins);
    h_H_j_pT_fine.book("H_j_pT_fine",H_pT_jj_fine_bins);
    h_H_j_pT_excl.book("H_j_pT_excl",H_pT_1j_bins);
    h_H_j_pT_fine_excl.book("H_j_pT_fine_excl",H_pT_jj_fine_bins);

    h_Hj_pT.book("Hj_pT",Hj_pT_bins);
    h_Hj_pT_fine.book("Hj_pT_fine",H_pT_jj_fine_bins);
    h_Hj_pT_excl.book("Hj_pT_excl",H_pT_1j_bins);
    h_Hj_pT_fine_excl.book("Hj_pT_fine_excl",H_pT_jj_fine_bins);

    h_H_jj_pT.book("H_jj_pT",H_pT_jj_bins);
    h_H_jj_pT_fine.book("H_jj_pT_fine",H_pT_jj_fine_bins);
    h_H_jj_pT_excl.book("H_jj_pT_excl",H_pT_2j_bins);
    h_H_jj_pT_fine_excl.book("H_jj_pT_fine_excl",H_pT_jj_fine_bins);

    h_Hjj_pT.book("Hjj_pT",Hjj_pT_bins);
    h_Hjj_pT_fine.book("Hjj_pT_fine",H_pT_jj_fine_bins);
    h_Hjj_pT_excl.book("Hjj_pT_excl",H_pT_2j_bins);
    h_Hjj_pT_fine_excl.book("Hjj_pT_fine_excl",H_pT_jj_fine_bins);

    h_dijet_mass.book("dijet_mass",dijet_mass_bins);
    h_dijet_mass_fine.book("dijet_mass_fine",dijet_mass_fine_bins);
    h_H_dijet_mass.book("H_dijet_mass",H_dijet_mass_bins);

    h_deltaphi_jj.book("deltaphi_jj",deltaphi_jj_bins);
    h_deltaphi_jj_fine.book("deltaphi_jj_fine",deltaphi_jj_fine_bins);
    h_deltaphi_jj_excl.book("deltaphi_jj_excl",deltaphi_jj_bins);
    h_deltaphi_Hjj.book("deltaphi_Hjj",deltaphi_Hjj_bins);
    h_deltaphi_Hjj_excl.book("deltaphi_Hjj_excl",deltaphi_Hjj_bins);
    h_deltaphi_jj_VBF.book("deltaphi_jj_VBF",deltaphi_jj_fine_bins);

    h_deltay_jj.book("deltay_jj",deltay_jj_bins);
    h_deltay_jj_fine.book("deltay_jj_fine",deltay_jj_fine_bins);
    h_deltay_H_jj.book("deltay_H_jj",deltay_H_jj_bins);
    h_deltay_H_jj_fine.book("deltay_H_jj_fine",deltay_jj_fine_bins);

    h_loose.book("loose",3,0.,3.);
    h_tight.book("tight",3,0.,3.);

    h_HT_jets_hist.book("HT_jets_hist",HT_jets_bins);

    h_tau_jet1.book("tau_jet1",tau_jet_bins);
    h_tau_jet2.book("tau_jet2",tau_jet_bins);
    h_tau_jet3.book("tau_jet3",tau_jet_bins);
    h_tau_jet_max.book("tau_jet_max",tau_jet_bins);
    h_sum_tau_jet.book("sum_tau_jet",tau_jet_bins);

  }

  // Reading events from the input TChain ***************************
  Long64_t numOK = 0;
  cout << "Preparing to read " << nEntries << " events" << endl;
  time_t last_time = time(0);
  unsigned seconds = 0;

  for (Long64_t ent = num_events.first; ent < num_events.second; ++ent) {

    // counter
    time_t cur_time = time(0);
    if ( timediff(cur_time-last_time) > 1 ) {
      cout << setw(10) << ent << setw(9) << seconds << 's';
      cout.flush();
      for (char i=0;i<20;i++) cout << '\b';
      last_time = cur_time;
      ++seconds;
    }

    tree->GetEntry(ent);

    if (nparticle>MAXNP) {
      cerr << "More particles in the event then MAXNP" << endl
           << "Increase array length to " << nparticle << endl;
      exit(1);
    }

    // map particle pdg id to index number in the array
    finder<Int_t> pdg(kf,nparticle);
    size_t hi; // Higgs index

    try {
      hi = pdg(25,1); // find Higgs

      numOK++; // count number of good events
    } catch (exception& e) {
      if (!quiet) {
        cerr << "In event " << ent << ": ";
        cerr << e.what() << endl;

        for (Int_t i=0;i<nparticle;i++) cerr << kf[i] << ' ';
        cerr << endl;

      }
      continue; // skip to next event
    }

    // Higgs 4-vector
    const TLorentzVector higgs(px[hi],py[hi],pz[hi],E[hi]);

    const Double_t H_mass = higgs.M();        // Higgs Mass
    const Double_t H_pt   = higgs.Pt();       // Higgs Pt
    const Double_t H_eta  = higgs.Rapidity(); // Higgs Rapidity

    // Fill histograms ***********************************
    for (Int_t i=0;i<nparticle;i++) h_pid.Fill(kf[i]);

    h_xs.Fill(0.5,weight);

    h_H_mass   .Fill(H_mass     ,weight);
    h_H_pt     .Fill(H_pt       ,weight);
    h_H_pt_fine.Fill(H_pt       ,weight);
    h_H_y      .Fill(abs(H_eta) ,weight);
    h_H_y_fine .Fill(abs(H_eta) ,weight);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Loop over SpartyJet clustering algorithms
    for (hist_map::it = SJClusterAlg::begin();
         hist_map::it!= SJClusterAlg::end(); ++hist_map::it)
    {
      // sort jets by Pt
      const vector<TLorentzVector> jets = (*hist_map::it)->jetsByPt(30.,5.);
      const size_t njets = jets.size(); // number of jets

      h_NJet_excl.Fill(njets,weight);
      for (unsigned char i=0;i<4;i++)
        if (njets>=i) h_NJet_incl.Fill(i,weight);

      if (njets==0) { // njets == 0;

        h_H_pT_excl     .Fill(H_pt,weight);
        h_H_pT_fine_excl.Fill(H_pt,weight);
        // 6/2 added fill for jet1_pT for 0-30 GeV bin, i.e. no jets
        h_jet1_pT       .Fill(10,weight);
        // 6/2 added fill for overflow bin for 0 jets in event
        h_deltay_jj     .FillOverflow(weight);
        // 6/2 added fill for overflow bin for 0 jets in event
        h_Hjj_pT        .FillOverflow(weight);
        //6/2 added fill for overflow bin for 0 jets in event
        h_jet2_y        .FillOverflow(weight);
        //6/2 added fill for overflow bin for 0 jets in event
        h_jet2_pT       .FillOverflow(weight);

        h_NJet_incl_50.Fill(0,weight);
        h_NJet_excl_50.Fill(0,weight);

      }
      else { // njets > 0;

        const TLorentzVector& j1 = jets[0]; // First jet

        const Double_t j1_mass = j1.M();
        const Double_t j1_pt   = j1.Pt();
        const Double_t j1_eta  = j1.Rapidity();

        const Double_t Hj_pt = (higgs + j1).Pt();

        h_jet1_pT     .Fill(j1_pt   ,weight);
        h_jet1_pT_fine.Fill(j1_pt   ,weight);
        h_jet1_y      .Fill(j1_eta  ,weight);
        h_jet1_y_fine .Fill(j1_eta  ,weight);
        h_Hj_pT       .Fill(Hj_pt   ,weight);
        h_Hj_pT_fine  .Fill(Hj_pt   ,weight);
        h_H_j_pT      .Fill(H_pt    ,weight);
        h_H_j_pT_fine .Fill(H_pt    ,weight);

        h_tau_jet1.Fill(
          sqrt( sq(j1_pt) + sq(j1_mass) )/( 2.0*cosh(j1_eta - H_eta) ),
          weight
        );

        if (j1_pt>50.) {
          h_NJet_incl_50.Fill(1,weight);
        } else {
          h_NJet_incl_50.Fill(0,weight);
          h_NJet_excl_50.Fill(0,weight);
        }

        if (njets==1) { // njets == 1;

          h_Hj_pT_excl        .Fill(Hj_pt ,weight);
          h_Hj_pT_fine_excl   .Fill(Hj_pt ,weight);
          h_H_j_pT_excl       .Fill(H_pt  ,weight);
          h_H_j_pT_fine_excl  .Fill(H_pt  ,weight);
          h_jet1_pT_excl      .Fill(j1_pt ,weight);
          h_jet1_pT_excl_fine .Fill(j1_pt ,weight);
          // 6/2 added fill for j2_pT for 0-30 GeV bins, i.e. no 2nd jet
          h_jet2_pT           .Fill(10,weight);
          // 6/2 added fill for overflow bin for 1 jet in event
          h_deltay_jj         .FillOverflow(weight);
          // 6/2 added fill for overflow bin for 1 jet in event
          h_Hjj_pT            .FillOverflow(weight);
          //6/2 added fill for overflow bin for 1 jet in event
          h_jet2_y            .FillOverflow(weight);

          if(j1_pt>50.) h_NJet_excl_50.Fill(1,weight);

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

          h_deltaphi_jj     .Fill(deltaPhi_j_j  ,weight);
          h_deltaphi_jj_fine.Fill(deltaPhi_j_j  ,weight);
          h_deltaphi_Hjj    .Fill(deltaPhi_H_jj ,weight);
          h_Hjj_pT          .Fill(Hjj_pt        ,weight);
          h_Hjj_pT_fine     .Fill(Hjj_pt        ,weight);
          h_H_jj_pT         .Fill(H_pt          ,weight);
          h_H_jj_pT_fine    .Fill(H_pt          ,weight);
          h_jet2_pT         .Fill(j2_pt         ,weight);
          h_jet2_pT_fine    .Fill(j2_pt         ,weight);
          h_jet2_y          .Fill(abs(j2_eta)   ,weight);
          h_jet2_y_fine     .Fill(abs(j2_eta)   ,weight);
          h_dijet_mass      .Fill(jj_mass       ,weight);
          h_dijet_mass_fine .Fill(jj_mass       ,weight);
          h_H_dijet_mass    .Fill(Hjj_mass      ,weight);
          h_deltay_jj       .Fill(deltay_j_j    ,weight);
          h_deltay_jj_fine  .Fill(deltay_j_j    ,weight);
          h_deltay_H_jj     .Fill(deltay_H_jj   ,weight);
          h_deltay_H_jj_fine.Fill(deltay_H_jj   ,weight);

          if (deltay_j_j>2.8) {
            if (jj_mass>400) {
              // 6/13 realized that delta-phi cut should be between Higgs and dijet system
              h_deltaphi_jj_VBF.Fill(deltaPhi_H_jj,weight);
              // 6/2 added loose and tight histograms to tally cross section as cross-checks
              h_loose.Fill(1,weight);
              if (deltaPhi_H_jj>2.6) h_tight.Fill(1,weight);
            }
          }

          h_tau_jet2.Fill(
            sqrt( sq(j2_pt) + sq(j2_mass) )/( 2.0*cosh(j2_eta - H_eta) ),
            weight
          );

          if (j2_pt>50.) h_NJet_incl_50.Fill(2,weight);

          if (njets==2) { // njets == 2;

            h_deltaphi_jj_excl  .Fill(deltaPhi_j_j  ,weight);
            h_deltaphi_Hjj_excl .Fill(deltaPhi_H_jj ,weight);
            h_Hjj_pT_excl       .Fill(Hjj_pt        ,weight);
            h_Hjj_pT_fine_excl  .Fill(Hjj_pt        ,weight);
            h_H_jj_pT_excl      .Fill(H_pt          ,weight);
            h_H_jj_pT_fine_excl .Fill(H_pt          ,weight);
            // 6/2 added fill for j3_pT 0-30 GeV, i.e. no jet3
            h_jet3_pT           .Fill(10            ,weight);

            if (j2_pt>50.) h_NJet_excl_50.Fill(2,weight);

          }
          else { // njets > 2;
            const TLorentzVector& j3 = jets[2]; // Second jet

            const Double_t j3_mass = j3.M();
            const Double_t j3_pt   = j3.Pt();
            const Double_t j3_eta  = j3.Rapidity();

            h_jet3_pT     .Fill(j3_pt   ,weight);
            h_jet3_pT_fine.Fill(j3_pt   ,weight);
            h_jet3_y      .Fill(j3_eta  ,weight);
            h_jet3_y_fine .Fill(j3_eta  ,weight);

            h_tau_jet3.Fill(
              sqrt( sq(j3_pt) + sq(j3_mass) )/( 2.0*cosh(j3_eta - H_eta) ),
              weight
            );

            if (j3_pt>50.) h_NJet_incl_50.Fill(3,weight);

            if (njets==3) {

              if (j3_pt>50.) h_NJet_excl_50.Fill(3,weight);

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
          sqrt( sq(jet_pt) + sq(jet_mass) )/( 2.0*cosh(jet_eta - H_eta) );

        if ( tauJet > tau_jet_cut ) {
          sum_tj += tauJet;
          if (tauJet > max_tj)
            max_tj = tauJet;
        }
        
      }
      h_HT_jets_hist.Fill(HT_jets,weight);
      h_tau_jet_max.Fill(max_tj,weight);
      h_sum_tau_jet.Fill(sum_tj,weight);

    } // END Loop over SpartyJet clustering algorithms

  } // END of event loop

  cout << setw(10) << nEntries << setw(9) << seconds <<'s' << endl;

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
