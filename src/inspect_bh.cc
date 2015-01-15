#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include <TFile.h>
#include <TTree.h>

#include "BHEvent.h"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;

int main(int argc, char** argv)
{
  if (argc!=4) {
    cout << "Usage: " << argv[0] << " bh.root first num" << endl;
    exit(0);
  }
  const Long64_t first = atol(argv[2]);
  const Long64_t num   = atol(argv[3]);

  TFile *bh = new TFile(argv[1],"read");
  if (bh->IsZombie()) exit(1);
  TTree *tree = (TTree*)bh->Get("t3");

  BHEvent event;
  event.SetTree(tree);

  cout << setw(4) << "id"
       << setw(4) << "kf"
       << endl;

  const Long64_t end = min(tree->GetEntries(),first+num);
  for (Long64_t ent = first; ent < end; ++ent) {
    tree->GetEntry(ent);

    cout << setw(4) << event.eid
         << setw(4) << event.kf
         << endl;

  }

  return 0;
}
