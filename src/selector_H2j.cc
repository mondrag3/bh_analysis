#include "selector.h"

#include <iostream>

#include "finder.h"

using namespace std;

struct hist_defs {
  TH1 *N, *pid;

  hist_wt xs, H_mass, H_pT, H_y;

  hist_alg_wt
    NJet_incl, NJet_excl, NJet_incl_50, NJet_excl_50,

    H_pT_excl,

    jet1_pT, jet1_pT_excl, jet1_y, jet1_tau,
    jet2_pT, jet2_y, jet2_tau,
    jet3_pT, jet3_y, jet3_tau,

    jj_mass,
    j_j_deltaphi, j_j_deltaphi_excl, j_j_deltaphi_VBF,
    j_j_deltay,

    H_j_pT, H_j_pT_excl,
    Hj_pT, Hj_pT_excl,

    H_jj_pT, H_jj_pT_excl, Hjj_mass,
    H_jj_deltaphi, H_jj_deltaphi_excl,
    H_jj_deltay,
    Hjj_pT, Hjj_pT_excl,

    loose, tight,
    jets_HT, jets_tau_max, jets_tau_sum
  ;

  #define h_(name) name(#name)

  hist_defs()
  : N(hist::css->mkhist("N")), pid(hist::css->mkhist("pid")),
    h_(xs), h_(H_mass), h_(H_pT), h_(H_y),
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
    h_(loose), h_(tight),
    h_(jets_HT), h_(jets_tau_max), h_(jets_tau_sum)
    { }
};

selector::selector(timed_counter *counter)
: selector_base(counter), h(nullptr)
{ }
selector::~selector() {
  delete h;
}

void selector::SlaveBegin() {
  h = new hist_defs;
}

template<typename T> inline T sq(const T& x) { return x*x; }

inline Double_t tau(Double_t j_pt, Double_t j_mass, Double_t j_eta, Double_t H_eta) {
  return sqrt( sq(j_pt) + sq(j_mass) )/( 2.*cosh(j_eta - H_eta) );
}

Bool_t selector::Process(Long64_t entry) {
  selector_base::Process(entry);
  static const BHEvent& event = BHEvent::event;

  // map particle pdg id to index number in the array
  finder<Int_t> pdg(event.kf,event.nparticle);
  size_t hi; // Higgs index

  try {
    hi = pdg(25,1); // find Higgs

    numOK++; // count number of good events
  } catch (exception& e) {
    // if (!quiet) {
      cerr << "In event " << entry << ": ";
      cerr << e.what() << entry;

      for (Int_t i=0;i<event.nparticle;i++) cerr << event.kf[i] << ' ';
      cerr << endl;

    // }
    // continue; // skip to next event
    return kTRUE;
  }

  h->N->Fill(0.5);

  // Higgs 4-vector
  const TLorentzVector higgs(event.px[hi],event.py[hi],event.pz[hi],event.E[hi]);

  const Double_t H_mass = higgs.M();        // Higgs Mass
  const Double_t H_pt   = higgs.Pt();       // Higgs Pt
  const Double_t H_eta  = higgs.Rapidity(); // Higgs Rapidity

  // Fill histograms ***********************************
  for (Int_t i=0;i<event.nparticle;i++) h->pid->Fill(event.kf[i]);

  h->xs.Fill(0.5);

  h->H_mass.Fill(H_mass);
  h->H_pT  .Fill(H_pt);
  h->H_y   .Fill(H_eta);

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Loop over SpartyJet clustering algorithms
  for (auto& alg : SJClusterAlg::all) {
    hist::alg_ptr = alg.get(); // set static algorithm pointer for histograms

    // sort jets by Pt
    // const vector<TLorentzVector> jets = alg->jetsByPt(pt_cut,eta_cut);
    const vector<TLorentzVector> jets = alg->jetsByPt(30.,4.4);
    const size_t njets = jets.size(); // number of jets

    int njets50 = 0;
    for (auto& j : jets) {
      if (j.Pt()>=50.) ++njets50;
      else break;
    }

    // Number of jets hists
    h->NJet_excl.Fill(njets);
    h->NJet_excl_50.Fill(njets50);
    for (unsigned char i=0;i<4;i++) {
      if (njets  >=i) h->NJet_incl   .Fill(i);
      if (njets50>=i) h->NJet_incl_50.Fill(i);
    }

    if (njets==0) { // njets == 0;

      h->H_pT_excl.Fill(H_pt);

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
      h->jets_HT.Fill(jets_HT);
      h->jets_tau_max.Fill(max_tj);
      h->jets_tau_sum.Fill(sum_tj);

      const TLorentzVector& j1 = jets[0]; // First jet

      const Double_t j1_mass = j1.M();
      const Double_t j1_pt   = j1.Pt();
      const Double_t j1_eta  = j1.Rapidity();

      const Double_t Hj_pt = (higgs + j1).Pt();

      h->jet1_pT.Fill(j1_pt);
      h->jet1_y .Fill(j1_eta);
      h->Hj_pT  .Fill(Hj_pt);
      h->H_j_pT .Fill(H_pt);

      h->jet1_tau.Fill( tau(j1_pt,j1_mass,j1_eta,H_eta) );

      if (njets==1) { // njets == 1;

        h->Hj_pT_excl  .Fill(Hj_pt);
        h->H_j_pT_excl .Fill(H_pt);
        h->jet1_pT_excl.Fill(j1_pt);

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

        h->j_j_deltaphi .Fill(deltaPhi_j_j);
        h->H_jj_deltaphi.Fill(deltaPhi_H_jj);
        h->Hjj_pT       .Fill(Hjj_pt);
        h->H_jj_pT      .Fill(H_pt);
        h->jet2_pT      .Fill(j2_pt);
        h->jet2_y       .Fill(j2_eta);
        h->jj_mass      .Fill(jj_mass);
        h->Hjj_mass     .Fill(Hjj_mass);
        h->j_j_deltay   .Fill(deltay_j_j);
        h->H_jj_deltay  .Fill(H_jj_deltay);

        if (deltay_j_j>2.8) {
          if (jj_mass>400) {
            h->j_j_deltaphi_VBF.Fill(deltaPhi_H_jj);
            h->loose.Fill(1);
            if (deltaPhi_H_jj>2.6) h->tight.Fill(1);
          }
        }

        h->jet2_tau.Fill( tau(j2_pt,j2_mass,j2_eta,H_eta) );

        if (njets==2) { // njets == 2;

          h->j_j_deltaphi_excl  .Fill(deltaPhi_j_j);
          h->H_jj_deltaphi_excl .Fill(deltaPhi_H_jj);
          h->Hjj_pT_excl        .Fill(Hjj_pt);
          h->H_jj_pT_excl       .Fill(H_pt);

        }
        else { // njets > 2;
          const TLorentzVector& j3 = jets[2]; // Second jet

          const Double_t j3_mass = j3.M();
          const Double_t j3_pt   = j3.Pt();
          const Double_t j3_eta  = j3.Rapidity();

          h->jet3_pT.Fill(j3_pt);
          h->jet3_y .Fill(j3_eta);

          h->jet3_tau.Fill( tau(j3_pt,j3_mass,j3_eta,H_eta) );

        } // END njets > 2;

      } // END njets > 1;

    } // END njets > 0;

  } // END Loop over SpartyJet clustering algorithms

  return kTRUE;
}
