#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TLatex.h>
#include <TAxis.h>

using namespace std;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

int main(int argc, char **argv)
{
  // Arguments ******************************************************
  if (argc!=2 && argc!=3) {
    cout << "Usage: " << argv[0] << " input.root [output.pdf]" << endl;
    exit(1);
  }

  string fin(argv[1]);
  string fout;
  if (argc==3) {
    fout = argv[2];
  } else {
    size_t first = fin.rfind('/')+1;
    if (first==string::npos) first = 0;
    fout = fin_name.substr(first,fin_name.rfind('.')-first)+".pdf";
  }

  // Collect ********************************************************
  vector<pair<string,vector<pair<string,TH1*>>>> hists;

  TFile *f = new TFile(fin.c_str(),"read");

  TIter nextkey_d(f->GetListOfKeys());
  TKey *key_d;
  while ((key_d = static_cast<TKey*>(nextkey()))) {
    TObject *obj = key_d->ReadObj();
    if (obj->InheritsFrom(TDirectory::Class())) {
      TDirectory *dir = static_cast<TDirectory*>(obj);

        TIter nextkey_h(dir->GetListOfKeys());
        TKey *key_h;
        while ((key_h = static_cast<TKey*>(nextkey()))) {
          TObject *obj = key_h->ReadObj();
          if (obj->InheritsFrom(TH1::Class())) {

            

          }
        }

    }
  }


  return 0;
}
