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
#include <algorithm>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "BHEvent.hh"
#include "SJClusterAlg.hh"
#include "weight.hh"
#include "timed_counter.hh"
#include "xml_analysis.hh"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

template<typename T> inline T sq(const T& x) { return x*x; }

// istream operators ************************************************
namespace std {
  template<class A, class B>
  istream& operator>> (istream& is, pair<A,B>& r) {
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
}

// Parse FastJet jet definition *************************************
fastjet::JetDefinition* JetDef(string& str) {
  string::iterator it = --str.end();
  while (isdigit(*it)) --it;
  ++it;
  string name;
  transform(str.begin(), it, back_inserter(name), ::tolower);
  fastjet::JetAlgorithm alg;
  if (!name.compare("antikt")) alg = fastjet::antikt_algorithm;
  else if (!name.compare("kt")) alg = fastjet::kt_algorithm;
  else if (!name.compare("cambridge")) alg = fastjet::cambridge_algorithm;
  else throw runtime_error("Undefined jet clustering algorithm: "+name);
  return new fastjet::JetDefinition(
    alg,
    atof( string(it,str.end()).c_str() )/10.
  );
}

// ******************************************************************

inline Double_t tau(Double_t j_pt, Double_t j_mass, Double_t j_eta, Double_t H_eta) {
  return sqrt( sq(j_pt) + sq(j_mass) )/( 2.*cosh(j_eta - H_eta) );
}

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, sj_files, wt_files, weights;
  string xml_file, output_file, css_file, jet_alg;
  pair<Long64_t,Long64_t> num_events;

  bool sj_given = false, wt_given = false;

  try {
    // General Options ------------------------------------
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "produce help message")
      ("analysis,a", po::value<string>(&xml_file)->required(),
       "*analysis XML file")
      ("bh", po::value< vector<string> >(&bh_files)->required(),
       "*add input BlackHat root file")
      ("sj", po::value< vector<string> >(&sj_files),
       "add input SpartyJet root file")
      ("wt", po::value< vector<string> >(&wt_files),
       "add input weights root file")
      ("output,o", po::value<string>(&output_file)->required(),
       "*output root file with histograms")
      ("cluster,c", po::value<string>(&jet_alg)->default_value("AntiKt4"),
       "jet clustering algorithm: e.g. antikt4, kt6\n"
       "without --sj: select FastJet algorithm\n"
       "with --sj: read jets from SpartyJet ntuple")
      ("weight,w", po::value<vector<string>>(&weights),
       "weight branchs; if skipped:\n"
       "  without --wt: ntuple weight is used\n"
       "  with --wt: all weights from wt files")
      ("style,s", po::value<string>(&css_file)
       ->default_value(CONFDIR"/Hj.css","Hj.css"),
       "CSS style file for histogram binning and formating")
      ("num-events,n", po::value<pair<Long64_t,Long64_t>>(&num_events),
       "process only this many events,\nnum or first:num")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (argc == 1 || vm.count("help")) {
      cout << desc << endl;
      return 0;
    }
    po::notify(vm);
    if (sj_files.size()) sj_given = true;
    if (wt_files.size()) wt_given = true;
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Setup input files **********************************************
  TChain*    tree = new TChain("t3");
  TChain* sj_tree = (sj_given ? new TChain("SpartyJet_Tree") : nullptr);
  TChain* wt_tree = (wt_given ? new TChain("weights") : nullptr);

  // Add trees from all the files to the TChains
  cout << "BH files:" << endl;
  for (auto& f : bh_files) {
    cout << "  " << f << endl;
    if (!tree->AddFile(f.c_str(),-1) ) exit(1);
  }
  if (sj_given) {
    cout << "SJ files:" << endl;
    for (auto& f : sj_files) {
      cout << "  " << f << endl;
      if (!sj_tree->AddFile(f.c_str(),-1) ) exit(1);
    }
  }
  if (wt_given) {
    cout << "Weight files:" << endl;
    for (auto& f : wt_files) {
      cout << "  " << f << endl;
      if (!wt_tree->AddFile(f.c_str(),-1) ) exit(1);
    }
  }
  cout << endl;

  // Find number of events to process
  if (num_events.second>0) {
    const Long64_t need_events = num_events.first + num_events.second;
    if (need_events>tree->GetEntries()) {
      cerr << "Fewer entries in BH chain (" << tree->GetEntries()
         << ") then requested (" << need_events << ')' << endl;
      exit(1);
    }
    if (sj_given) if (need_events>sj_tree->GetEntries()) {
      cerr << "Fewer entries in SJ chain (" << sj_tree->GetEntries()
         << ") then requested (" << need_events << ')' << endl;
      exit(1);
    }
    if (wt_given) if (need_events>wt_tree->GetEntries()) {
      cerr << "Fewer entries in weights chain (" << wt_tree->GetEntries()
         << ") then requested (" << need_events << ')' << endl;
      exit(1);
    }
  } else {
    num_events.second = tree->GetEntries();
    if (sj_given) if (num_events.second!=sj_tree->GetEntries()) {
      cerr << num_events.second << " entries in BH chain, but "
           << sj_tree->GetEntries() << " entries in SJ chain" << endl;
      exit(1);
    }
    if (wt_given) if (num_events.second!=wt_tree->GetEntries()) {
      cerr << num_events.second << " entries in BH chain, but "
           << wt_tree->GetEntries() << " entries in weights chain" << endl;
      exit(1);
    }
  }

  // Friend BlackHat tree with SpartyJet and Weight trees
  if (sj_given) tree->AddFriend(sj_tree,"SJ");
  if (wt_given) tree->AddFriend(wt_tree,"weights");

  // BlackHat tree branches
  BHEvent event;
  event.SetTree(tree, BHEvent::kinematics);

  // Jet Clustering Algorithm
  unique_ptr<fastjet::JetDefinition> jet_def;
  unique_ptr<SJClusterAlg> sj_alg;

  if (sj_given) {
    sj_alg.reset( new SJClusterAlg(tree,jet_alg) );
  } else {
    jet_def.reset( JetDef(jet_alg) );
    cout << "Clustering with " << jet_def->description() << endl << endl;
  }

  // Weights tree branches
  if (wt_given) {
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
  } else weight::add(tree,"weight",false); // Use default ntuple weight
  cout << endl;

  // Read CSS file with histogram properties
  cout << "Histogram CSS file: " << css_file << endl;
  hist::css.reset( new csshists(css_file) );
  cout << endl;

  // Open output file with histograms *******************************
  TFile* fout = new TFile(output_file.c_str(),"recreate");
  if (fout->IsZombie()) exit(1);
  else cout << "Output file: " << fout->GetName() << endl << endl;

  // Make directories ***********************************************
  for (auto& w : weight::all) {
    hist_wt::dirs[w.get()] = fout->mkdir(w->name.c_str());
  }

  fout->cd();

  // Book histograms ************************************************
  TH1* h_N   = hist::css->mkhist("N");
  TH1* h_pid = hist::css->mkhist("pid");

  #define h_(name) h_##name(#name)

  /* NOTE:
   * excl = exactly the indicated number of jets, zero if no j in name
   * incl = that many or more jets
   * VBF = vector boson fusion cut
   */

  // Book Histograms
  hist_wt
    h_(xs),

    h_(NJet_incl), h_(NJet_excl), h_(NJet_incl_50), h_(NJet_excl_50),

    h_(H_mass), h_(H_pT), h_(H_y),

    h_(H_0j_pT), h_(H_0+j_pT),

    h_(H_1j_pT), h_(H_1+j_pT), h_(jet1_pt),

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
  Int_t prev_id = -1;
  cout << "Reading " << num_events.second << " entries";
  if (num_events.first>0) cout << " starting at " << num_events.first << endl;
  else cout << endl;
  num_events.second += num_events.first;
  timed_counter counter;

  for (Long64_t ent = num_events.first; ent < num_events.second; ++ent) {
    counter(ent);
    tree->GetEntry(ent);

    if (event.nparticle>BHMAXNP) {
      cerr << "More particles in the event then BHMAXNP" << endl
           << "Increase array length to " << event.nparticle << endl;
      exit(1);
    }

    // Find Higgs
    Int_t hi = 0; // Higgs index
    while (hi<event.nparticle) {
      if (event.kf[hi]==25) break;
      else ++hi;
    }
    if (hi==event.nparticle) {
      cerr << "No Higgs in event " << ent << endl;
      continue;
    }
    ++numOK;

    // Count number of events (not entries)
    if (prev_id!=event.eid) h_N->Fill(0.5);
    prev_id = event.eid;

    // Higgs 4-vector
    const TLorentzVector higgs(event.px[hi],event.py[hi],event.pz[hi],event.E[hi]);

    const Double_t H_mass = higgs.M();   // Higgs Mass
    const Double_t H_pT   = higgs.Pt();  // Higgs Pt
    const Double_t H_eta  = higgs.Eta(); // Higgs Pseudo-rapidity

    // Fill histograms ***********************************
    for (Int_t i=0;i<event.nparticle;i++) h_pid->Fill(event.kf[i]);

    h_xs.Fill(0.5);

    h_H_mass.Fill(H_mass);
    h_H_pT  .Fill(H_pT);
    h_H_y   .Fill(H_eta);

    // Jet clustering *************************************
    vector<TLorentzVector> jets;
    if (sj_given) { // Read jets from SpartyJet ntuple
      jets = sj_alg->jetsByPt(pt_cut,eta_cut);

    } else { // Clusted with FastJet on the fly
      vector<fastjet::PseudoJet> particles;
      particles.reserve(event.nparticle-1);

      for (Int_t i=0; i<event.nparticle; ++i) {
        if (i==hi) continue;
        particles.emplace_back(
          event.px[i],event.py[i],event.pz[i],event.E[i]
        );
      }

      // Cluster, sort jets by pT, and apply pT cut
      const vector<fastjet::PseudoJet> fj_jets = sorted_by_pt(
        fastjet::ClusterSequence(particles, *jet_def).inclusive_jets(pt_cut)
      );

      // Apply eta cut
      jets.reserve(fj_jets.size());
      for (auto& j : fj_jets) {
        if (abs(j.eta()) < eta_cut)
          jets.emplace_back(j.px(),j.py(),j.pz(),j.E());
      }
    }
    const size_t njets = jets.size(); // number of jets

    // ****************************************************

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

      h_H_pT_excl.Fill(H_pT);

    }
    else { // njets > 0;

      Double_t jets_HT = 0;
      static const Double_t jet_tau_cut=8;
      Double_t max_tj=0;
      Double_t sum_tj=0;

      for (auto& jet : jets) {
        const Double_t jet_pt  = jet.Pt();
        const Double_t jet_tau = tau(jet_pt,jet.M(),jet.Eta(),H_eta);

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
      const Double_t j1_eta  = j1.Eta();

      const Double_t Hj_pt = (higgs + j1).Pt();

      h_jet1_pT.Fill(j1_pt);
      h_jet1_y .Fill(j1_eta);
      h_Hj_pT  .Fill(Hj_pt);
      h_H_j_pT .Fill(H_pT);

      h_jet1_tau.Fill( tau(j1_pt,j1_mass,j1_eta,H_eta) );

      if (njets==1) { // njets == 1;

        h_Hj_pT_excl  .Fill(Hj_pt);
        h_H_j_pT_excl .Fill(H_pT);
        h_jet1_pT_excl.Fill(j1_pt);

      }
      else { // njets > 1;

        const TLorentzVector& j2 = jets[1]; // Second jet

        const Double_t j2_mass = j2.M();
        const Double_t j2_pt   = j2.Pt();
        const Double_t j2_eta  = j2.Eta();

        const TLorentzVector jj(j1+j2);
        const TLorentzVector Hjj(higgs + jj);

        const Double_t jj_mass        = jj.M();
        const Double_t deltaPhi_j_j   = j1.Phi() - j2.Phi();
        const Double_t deltaPhi_H_jj  = higgs.Phi() - jj.Phi();
        const Double_t Hjj_pt         = Hjj.Pt();
        const Double_t Hjj_mass       = Hjj.M();
        const Double_t deltay_j_j     = abs(j1_eta - j2_eta);
        const Double_t H_jj_deltay    = abs(H_eta-jj.Eta());

        h_j_j_deltaphi .Fill(deltaPhi_j_j);
        h_H_jj_deltaphi.Fill(deltaPhi_H_jj);
        h_Hjj_pT       .Fill(Hjj_pt);
        h_H_jj_pT      .Fill(H_pT);
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
          h_H_jj_pT_excl       .Fill(H_pT);

        }
        else { // njets > 2;
          const TLorentzVector& j3 = jets[2]; // Second jet

          h_H_jjj_pT.Fill(H_pT);

          const Double_t j3_mass = j3.M();
          const Double_t j3_pt   = j3.Pt();
          const Double_t j3_eta  = j3.Eta();

          h_jet3_pT.Fill(j3_pt);
          h_jet3_y .Fill(j3_eta);

          h_jet3_tau.Fill( tau(j3_pt,j3_mass,j3_eta,H_eta) );

        } // END njets > 2;

      } // END njets > 1;

    } // END njets > 0;

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
