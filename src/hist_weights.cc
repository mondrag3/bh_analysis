#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TBranch.h>
#include <TH1.h>

#include "hist_wrap.h"

using namespace std;

int main(int argc, char** argv)
{
  if (argc!=4) {
    cout << "Usage: " << argv[0]
         << " weights.root binning.bins hist.root" << endl;
    return 0;
  }

  TFile *fin = new TFile(argv[1],"read");
  TTree *tree = (TTree*)fin->Get("weights");
  if (!tree) {
    cerr << "No tree weights in file " << fin->GetName() << endl;
    fin->ls();
    exit(1);
  }
  const TObjArray *brarr = tree->GetListOfBranches();
  const size_t numbr = brarr->GetEntries();

  hist::read_binnings(argv[2]);

  TFile *fout = new TFile(argv[3],"recreate");

  vector<Float_t> x(numbr);
  vector<hist*> h;

  for (size_t i=0;i<numbr;++i) {
    TBranch *br = dynamic_cast<TBranch*>(brarr->At(i));
    string name = br->GetName();
    cout << "Branch: " << name << endl;
    br->SetAddress(&x[i]);
    h.push_back(new hist(name));
  }

  const Long64_t nent = tree->GetEntries();
  time_t time1=time(0), time2;
  unsigned seconds = 0;

  cout << "Prepared to read " << nent << " entries" << endl;
  for (Long64_t ent=0; ent<nent; ++ent) {

    tree->GetEntry(ent);
    for (size_t i=0;i<numbr;++i) h[i]->Fill(x[i]);

    // timed counter
    if ( difftime(time2=time(0),time1) > 1 ) {
      ++seconds;
      cout << setw(10) << ent
           << setw( 7) << seconds << 's';
      cout.flush();
      for (char i=0;i<18;++i) cout << '\b';
      time1 = time2;
    }
  }
  cout << setw(10) << nent
       << setw( 7) << seconds << 's' << endl << endl;

  fout->Write();
  fout->Close();
  delete fout;

  fin->Close();
  delete fin;

  return 0;
}
