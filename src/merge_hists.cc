#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>

using namespace std;
namespace po = boost::program_options;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

template<class T>
inline T* get_obj(TDirectory *d, const char* name) noexcept {
  T *obj = (T*)d->Get(name);
  if (obj) return obj;
  else {
    cerr << "\033[31mError: no object \"" << name
         << " in \"" << d->GetName() << "\"\033[0m" << endl;
    exit(1);
  }
}

struct hist_tree {
  TDirectory *d;
  vector<TH1*> h;
  vector<hist_tree> t;
} htree;

bool same_part = false;
Double_t N = 0.;

void merge_first(TDirectory *from, TDirectory *to, hist_tree& t=htree) {
  to->cd();
  t.d = to;
  static TKey *key;
  TIter nextkey(from->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    static TObject *obj;
    obj = key->ReadObj();
    if (obj->InheritsFrom(TDirectory::Class())) {
      t.t.emplace_back();
      merge_first(static_cast<TDirectory*>(obj),
                  to->mkdir(obj->GetName()), t.t.back());
      to->cd();
    } else if (obj->InheritsFrom(TH1::Class())) {
      t.h.push_back(static_cast<TH1*>( obj->Clone() ));
      if (!same_part) t.h.back()->Scale(1/N);
    } else {
      cerr << "Unexpected object class " << obj->ClassName() << endl;
      exit(1);
    }
  }
}
void merge(TDirectory *from, const hist_tree& t=htree) {
  for (auto&  h : t.h) {
    if (same_part) h->Add( get_obj<TH1>(from,h->GetName()) );
    else {
      TH1* _h = get_obj<TH1>(from,h->GetName());
      _h->Scale(1/N);
      h->Add(_h);
    }
  }
  for (auto& _t : t.t) merge( get_obj<TDirectory>(from,_t.d->GetName()), _t );
}

void scale_all(const hist_tree& t=htree) {
  for (auto&  h : t.h) h->Scale(1/N);
  for (auto& _t : t.t) scale_all(_t);
}

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> fin_name;
  string fout_name;

  try {
    // General Options ------------------------------------
    po::options_description desc("Options");
    desc.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<vector<string>>(&fin_name)->required()->multitoken(),
     "input file with histograms")
    ("output,o", po::value<string>(&fout_name)->required(),
     "output pdf plots")
    ("same-part,s", po::bool_switch(&same_part),
     "true : add, then scale\nfalse: scale, then add")
    ;
    po::positional_options_description desc_pos;
    desc_pos.add("input", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc)
              .positional(desc_pos).run(), vm);

    if (argc == 1 || vm.count("help")) {
      cout << desc << endl;
      return 0;
    }
    po::notify(vm);
  }
  catch(exception& e) {
    cerr << "\033[31mError: " << e.what() <<"\033[0m"<< endl;
    exit(1);
  }

  if ( find(fin_name.begin(),fin_name.end(),fout_name) != fin_name.end() ) {
    cerr << "\033[31mError: output file name is in the input files list\033[0m"
         << endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  TFile *fout = new TFile(fout_name.c_str(),"read");
  if (fout->IsZombie()) exit(1);
  cout << "Output file: " << fout->GetName() << endl;

  for (auto& fname : fin_name) {
    static bool first = true;

    TFile *fin = new TFile(fname.c_str(),"read");
    if (fin->IsZombie()) exit(1);
    cout << "Input file: " << fin->GetName() << endl;

    if (same_part) N += get_obj<TH1>(fin,"N")->GetBinContent(1);
    else N = get_obj<TH1>(fin,"N")->GetBinContent(1);

    if (first) merge_first(fin,fout);
    else merge(fin);

    first = false;
  }

  if (same_part) scale_all();

  fout->Write();
  cout << "Wrote " << fout->GetName() << endl;
  fout->Close();
  delete fout;

  return 0;
}
