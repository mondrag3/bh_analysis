#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>

using namespace std;

Double_t N;
short nj;

template<class T>
T* get_obj(TDirectory *d, const char* name) noexcept {
  T *obj = dynamic_cast<T*>(d->Get(name));
  if (obj) return obj;
  else {
    cout << "No " << name << " in " << d->GetName() << endl;
    exit(1);
  }
}

void sigma(TDirectory *d) {
  for (short j=0; j<4; ++j) {
    stringstream ss;
    ss << "H";
    if (j) {
      ss << '_';
      for (short i=0;i<j;++i) ss << 'j';
    }
    ss << "_pT";

    TH1 *h = (TH1*)d->Get(ss.str().c_str());
    if (!h) continue;

    cout << "H+" << j << "jets σ = "
         << showpoint << setprecision(6)
         << h->Integral(0,h->GetNbinsX()+1)/N
         << " pb" << endl;
  }
}

void read(TDirectory *d) noexcept {
  cout << d->GetName() << endl;
  static TKey *key;
  TIter nextkey(d->GetListOfKeys());
  int nsub = 0;
  while ((key = (TKey*)nextkey())) {
    static TObject *obj;
    obj = key->ReadObj();
    if (obj->InheritsFrom(TDirectory::Class())) {
      read(static_cast<TDirectory*>(obj));
      ++nsub;
    }
  }
  if (nsub==0) sigma(d);
}

int main(int argc, char** argv)
{
  if (argc!=2) {
    cout << "Usage: " << argv[0] << " hists.root" << endl;
    exit(0);
  }

  TFile *fin = new TFile(argv[1],"read");
  if (fin->IsZombie()) exit(1);
  cout << "File: " << fin->GetName() << endl;

  N = get_obj<TH1>(fin,"N")->GetBinContent(1);
  cout << "Events: " << N << endl;
  cout << endl;

  read(fin);

  return 0;
}
