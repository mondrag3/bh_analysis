#include <iostream>
#include <iomanip>
#include <string>

#include <boost/program_options.hpp>

#include "TFile.h"
#include "TTree.h"

#include "rew_calc.h"
#include "timed_counter.h"

using namespace std;
namespace po = boost::program_options;

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

BHEvent event; // extern

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  string BH_file, weights_file, pdf_set;
  bool old_bh, counter_newline;
  Long64_t num_events = 0;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
      ("help,h", "produce help message")
      ("bh", po::value<string>(&BH_file)->required(),
       "input event root file (Blackhat ntuple)")
      ("weights,o", po::value<string>(&weights_file)->required(),
       "output root file with new event weights")
      ("pdf", po::value<string>(&pdf_set)->default_value("CT10nlo"),
       "LHAPDF set name")
      ("num-events,n", po::value<Long64_t>(&num_events),
       "process only this many events. Zero means all.")
      ("old-bh", po::bool_switch(&old_bh),
       "read an old BH tree (no part & alphas_power branches)")
      ("counter-newline", po::bool_switch(&counter_newline),
       "do not overwrite previous counter message")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, all_opt), vm);
    if (argc == 1 || vm.count("help")) {
      cout << all_opt << endl;
      exit(0);
    }
    po::notify(vm);
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Open input event file
  TFile *fin = new TFile(BH_file.c_str(),"READ");
  if (fin->IsZombie()) exit(1);

  cout << "Input BH event file: " << fin->GetName() << endl;

  TTree *tin = (TTree*)fin->Get("t3");

  // Find number of events to process
  if (num_events>0) num_events = min( num_events, tin->GetEntries() );
  else num_events = tin->GetEntries();

  // Set up BlackHat event
  event.SetTree(tin, BHEvent::reweighting, old_bh);

  tin->SetBranchAddress("weight", &event.weight);

  if (old_bh) {

    if      (BH_file.find("born")!=string::npos) event.SetPart('B');
    else if (BH_file.find("real")!=string::npos) event.SetPart('R');
    else if (BH_file.find("loop")!=string::npos) event.SetPart('V');
    else if (BH_file.find("vsub")!=string::npos) event.SetPart('I');
    else {
      cerr << "\033[31mCannot determine part type from file name\033[0m" << endl;
      exit(1);
    }

    const size_t jpos = BH_file.find_first_of('j');
    if (jpos!=string::npos) {
      const char njc = BH_file[jpos-1];
      if (isdigit(njc)) {

        Char_t nj = njc - '0';
        if (event.part[0]=='V' || event.part[0]=='I') ++nj;
        event.SetAlphasPower(nj);

      } else {
        cerr << "\033[31mCannot determine number of jets from file name\033[0m" << endl;
        exit(1);
      }
    } else {
      cerr << "\033[31mCannot determine number of jets from file name\033[0m" << endl;
      exit(1);
    }
  }

  // Load PDF
  cout << endl;
  usePDFset(pdf_set);
  cout << endl;

  // Open output weights file
  TFile *fout = new TFile(weights_file.c_str(),"recreate");
  if (fout->IsZombie()) exit(1);

  cout << "Output weights file: " << fout->GetName() << endl;

  TTree *tree = new TTree("weights","");

  // pointers to calc (automatically collected)
  // (as well as argument pointers)
  auto FacHt1 = mk_fac_calc(new mu_fHt_Higgs(1.));
  auto FacHt2 = mk_fac_calc(new mu_fHt_Higgs(0.5)/*,true*/);
  auto FacHt4 = mk_fac_calc(new mu_fHt_Higgs(0.25));

  auto RenHt1 = mk_ren_calc(new mu_fHt_Higgs(1.),alphas_fcn::two_mH);
  auto RenHt2 = mk_ren_calc(new mu_fHt_Higgs(0.5),alphas_fcn::two_mH);
  auto RenHt4 = mk_ren_calc(new mu_fHt_Higgs(0.25),alphas_fcn::two_mH);

  auto Fac2MH = mk_fac_calc(new mu_const(125.*2.));
  auto FacMH  = mk_fac_calc(new mu_const(125.));
  auto FacMH2 = mk_fac_calc(new mu_const(125./2.));

  auto Ren2MH = mk_ren_calc(new mu_const(125.*2.),alphas_fcn::two_mH);
  auto RenMH  = mk_ren_calc(new mu_const(125.),alphas_fcn::two_mH);
  auto RenMH2 = mk_ren_calc(new mu_const(125./2.),alphas_fcn::two_mH);

  // auto FacDef = mk_fac_calc(new mu_fac_default());
  // auto RenDef = mk_ren_calc(new mu_ren_default(),alphas_fcn::two_mH);

  // define reweighting scales combinatios
  // and add branches to tree
  vector<reweighter*> rew {
    // new reweighter(FacHt2,RenHt2,tree/*,true*/),
    // new reweighter(FacDef,RenDef,tree/*,true*/),
    // new reweighter(FacHt2,RenHt1,tree),
    // new reweighter(FacHt2,RenHt4,tree),
    // new reweighter(FacHt4,RenHt2,tree),
    // new reweighter(FacHt4,RenHt4,tree),
    // new reweighter(FacHt1,RenHt1,tree),
    // new reweighter(FacHt1,RenHt2,tree)
    // new reweighter(FacMH,RenMH,tree)

    new reweighter(FacHt4,RenHt4,tree),
    new reweighter(FacHt2,RenHt2,tree),
    new reweighter(FacHt1,RenHt1,tree),

    new reweighter(Fac2MH,Ren2MH,tree),
    new reweighter(FacMH, RenMH, tree),
    new reweighter(FacMH2,RenMH2,tree)
  };

  // Reading events from the input ntuple ***************************
  cout << "\nReading " << num_events << " events" << endl;
  timed_counter counter(counter_newline);

  cout << scientific;
  cout.precision(10);

  for (Long64_t ent=0; ent<num_events; ++ent) {
    counter(ent);
    tin->GetEntry(ent);

    // test(ent)
    // test(event.weight)

    // use event id for event number
    event.eid = ent;

    // REWEIGHTING
    calc_all_scales();
    for (auto r : rew) r->stitch();
    tree->Fill();
  }
  counter.prt(num_events);
  cout << endl;

  fout->Write();
  cout << "\n\033[32mWrote\033[0m: " << fout->GetName() << endl;
  fout->Close();
  delete fout;

  fin->Close();
  delete fin;

  for (auto r : rew) delete r;

  return 0;
}
