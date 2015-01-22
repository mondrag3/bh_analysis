#include <iostream>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>

#include "BHEvent.hh"

using namespace std;

int main(int argc, char** argv)
{
  if (argc!=2) {
    cout << "Usage: " << argv[0] << " bh_ntuple.root" << endl;
    exit(0);
  }

  TFile *fin  = new TFile(argv[1],"read");
  if (fin->IsZombie()) exit(1);
  TTree *tree = (TTree*)fin->Get("t3");

  // BlackHat tree branches
  BHEvent event;
  event.SetTree(tree, BHEvent::cross_section);

  Double_t sigma = 0.;

  const Long64_t nent = tree->GetEntries();
  cout << "Entries: " << nent << endl;
  for (Long64_t ent = 0; ent < nent; ++ent) {
    tree->GetEntry(ent);
    sigma += event.weight;
  }
  sigma /= nent;

  cout << "Cross section: "
       << showpoint << setprecision(6) << sigma
       << " pb" << endl;

  delete fin;

  return 0;
}
