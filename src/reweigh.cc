#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <ctime>

#include <boost/program_options.hpp>

#include "TFile.h"
#include "TTree.h"

#include "rew_calc.h"

using namespace std;
namespace po = boost::program_options;

BHEvent event; // extern

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  string BH_file, weights_file, pdf_set, conf_file_name;
  bool old_bh, continue_on_err, counter_newline;
  Long64_t num_events = 0;
  vector<Float_t> _fHt_fac, _fHt_ren;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
      ("help,h", "produce help message")
      ("bh", po::value<string>(&BH_file),
       "input event root file (Blackhat ntuple)")
      ("weights,o", po::value<string>(&weights_file),
       "output root file with new event weights")
      ("pdf", po::value<string>(&pdf_set),
       "LHAPDF set name")
      ("num-events,n", po::value<Long64_t>(&num_events),
       "process only this many events. Zero means all.")
    /*("fac", po::value< vector<Float_t> >(&_fHt_fac),
       "add a factorization scale")
      ("ren", po::value< vector<Float_t> >(&_fHt_ren),
       "add a renormalization scale")*/
      ("old-bh", po::bool_switch(&old_bh),
       "read an old BH tree (no part & alphas_power branches)")
      ("continue-on-err", po::bool_switch(&continue_on_err),
       "continue after faulty event; default is break")
      ("counter-newline", po::bool_switch(&counter_newline),
       "do not overwrite previous counter message")
      ("conf,c", po::value<string>(&conf_file_name),
       "read a configuration file")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, all_opt), vm);
    po::notify(vm);
    if (vm.count("conf")) {
      ifstream conf_file(conf_file_name.c_str());
      po::store(po::parse_config_file<char>(conf_file, all_opt), vm);
    }
    po::notify(vm);

    // Options Properties ---------------------------------
    if (argc == 1 || vm.count("help")) {
      cout << all_opt << endl;
      exit(0);
    }

    // Necessary options ----------------------------------
    vector<string> rec_opt;
    rec_opt.push_back("bh");
    rec_opt.push_back("weights");
    rec_opt.push_back("pdf");
    //rec_opt.push_back("fac");
    //rec_opt.push_back("ren");

    for (size_t i=0, size=rec_opt.size(); i<size; ++i) {
      if (!vm.count(rec_opt[i]))
      { cerr << "Missing command --" << rec_opt[i] << endl; exit(1); }
    }
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Open input event file
  TFile *fin = new TFile(BH_file.c_str(),"READ");
  if (fin->IsZombie()) exit(1);

  TTree *tin = (TTree*)fin->Get("t3");

  // Find number of events to process
  if (num_events>0) num_events = min( num_events, tin->GetEntries() );
  else num_events = tin->GetEntries();

  // Set up BlackHat event
  event.SetTree(tin, BHEvent::reweighting);
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

  // Open output weights file
  TFile *fout = new TFile(weights_file.c_str(),"recreate");
  if (fout->IsZombie()) exit(1);

  TTree *tree = new TTree("weights","");

  // pointers to calc (automatically collected)
  // (as well as argument pointers)
  auto FacHt4 = mk_fac_calc(new mu_fHt(0.25));
  auto RenHt4 = mk_ren_calc(new mu_fHt(0.5));
  auto FacHt2 = mk_fac_calc(new mu_fHt(0.25));
  auto RenHt2 = mk_ren_calc(new mu_fHt(0.5));

  // define reweighting scales combinatios
  // and add branches to tree
  vector<reweighter> rew {
    reweighter(FacHt2,RenHt2,tree),
    reweighter(FacHt2,RenHt4,tree),
    reweighter(FacHt4,RenHt2,tree)
  };

  // Reading events from the input ntuple ***************************
  cout << "Prepared to read " << num_events << " events" << endl;
  time_t time1=time(0), time2;
  unsigned seconds = 0;
  bool complete = true;

  for (Long64_t ent=0; ent<num_events; ++ent) {
    tin->GetEntry(ent);

    // REWEIGHTING
    calc_all_scales();
    for (auto& r : rew) r.stitch();
    tree->Fill();

    // timed counter
    if ( difftime(time2=time(0),time1) > 1 ) {
      ++seconds;
      cout << setw(10) << ent
           << setw( 7) << seconds << 's';
      cout.flush();
      if (counter_newline) cout << endl;
      else for (char i=0;i<18;++i) cout << '\b';
      time1 = time2;
    }
  }
  cout << setw(10) << num_events
       << setw( 7) << seconds << 's' << endl << endl;

  fout->Write();
  fout->Close();
  delete fout;

  fin->Close();
  delete fin;

  // TODO: Implement error reporting in rew_calc
  if (complete) cout << "\033[32mComplete!\033[0m" << endl;
  else cout << "\033[31mIncomplete!\033[0m" << endl;

  return 0;
}
