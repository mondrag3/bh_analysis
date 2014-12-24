#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <TFile.h>
#include <TH1.h>

using namespace std;

template<class T>
T* get_obj(TFile *fin, const char* name) noexcept {
  T *obj = dynamic_cast<T*>(fin->Get(name));
  if (obj) return obj;
  else {
    cout << "No " << name << " in file " << fin->GetName() << endl;
    exit(1);
  }
}

int main(int argc, char** argv)
{
  if (argc!=3) {
    cout << "Usage: " << argv[0] << " n_jets hists.root" << endl;
    exit(0);
  }

  TFile *fin  = new TFile(argv[2],"read");
  if (fin->IsZombie()) exit(1);

  const Double_t N = get_obj<TH1>(fin,"N")->GetBinContent(1);

  for (short j=0, nj=atoi(argv[1]); j<=nj; ++j) {
    stringstream ss;
    ss << "H_";
    if (j) {
      ss << '_';
      for (short i=0;i<j;++i) ss << 'j';
    }
    ss << "_pT";

    TH1 *H_pT = get_obj<TH1>(fin,ss.str().c_str());

    cout << "H+" << j << "jets σ = "
         << showpoint << setprecision(6)
         << H_pT->Integral(0,H_pT->GetNbinsX()+1)/N
         << " pb" << endl;
  }

  return 0;
}
