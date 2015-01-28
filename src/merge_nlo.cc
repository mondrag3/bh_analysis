#include <iostream>
#include <string>
#include <unordered_map>

#include <TClass.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>

using namespace std;

// Safely get an object from TDirectory
template<class T>
inline T* get(TDirectory* dir, const char* name) {
  T *obj = nullptr;
  dir->GetObject(name,obj);
  if (obj) return obj;
  else {
    cerr << "No " << T::Class()->ClassName() << ' ' << name
         << " in " << dir->GetName() << endl;
    exit(1);
  }
}

int main(int argc, char** argv)
{
  if (argc<3) {
    cout << "Usage: " << argv[0]
         << " output_nlo.root input_born.root input_real.root ..." << endl;
    exit(0);
  }
  
  // Track histograms
  unordered_map<string,int> hist_names;

  // Output file
  TFile *nlo = new TFile(argv[1],"recreate");
  if (nlo->IsZombie()) exit(1);
  cout << "Output file: " << nlo->GetName() << endl;

  for (int fi=2; fi<argc; ++fi) {
    static bool first_file = true;

    // Current input file
    TFile *f = new TFile(argv[fi],"read");
    if (f->IsZombie()) exit(1);
    cout << "\nInput file: " << f->GetName() << endl;

    // Number of events
    const Double_t N = ((TH1*)f->Get("N"))->GetBinContent(1);
    const Double_t N_scale = 1./N;
    cout << "Events: " << N << endl;

    static TKey *key1;
    TIter nextkey1(f->GetListOfKeys());
    while ((key1 = (TKey*)nextkey1())) { // loop over dirs
      static TObject *obj1;
      obj1 = key1->ReadObj();
      if (obj1->InheritsFrom(TDirectory::Class())) {

        static bool first_dir = true;

        TDirectory *dir = static_cast<TDirectory*>(obj1);
        TDirectory *nlo_dir =
          ( first_file ? nlo->mkdir(dir->GetName())
                       : get<TDirectory>(nlo,dir->GetName()) );
        nlo_dir->cd();

        static TKey *key2;
        TIter nextkey2(dir->GetListOfKeys());
        while ((key2 = (TKey*)nextkey2())) { // loop over hists
          static TObject *obj2;
          obj2 = key2->ReadObj();
          if (obj2->InheritsFrom(TH1::Class())) {

            TH1 *h = static_cast<TH1*>(obj2);
            
            if (first_dir) {
              hist_names[h->GetName()] = 1;
            } else {
              const auto it = hist_names.find(h->GetName());
              if (it==hist_names.end()) {
                cerr << "Unexpected histogram " << h->GetName()
                     << " in file " << f->GetName()
                     << " in dir " << dir->GetName();
                exit(1);
              } else {
                ++it->second;
              }
            }
            
            // Scale and add histograms
            if (first_file) {
              h->Scale(N_scale);
              h->Clone();
            } else {
              h->Scale(N_scale);
              get<TH1>(dir,h->GetName())->Add(h,N_scale);
            }

          }
        }
        
        // Check if all histograms from the first directory were present
        auto it = hist_names.begin();
        const int first_counts = it->second;
        for (auto end=hist_names.end(); it!=end; ++it) {
          if (it->second != first_counts) {
            cerr << "Expected hist "
                 << (first_counts < it->second ? it->first
                                               : hist_names.begin()->first)
                 << " is not in file " << f->GetName()
                 << " in dir " << dir->GetName();
            exit(1);
          }
        }
        
        // no longer first directory
        if (first_dir) first_dir = false;
      }
    }

    delete f;
    // no longer first file
    if (first_file) first_file = false;
  }
  
  // Write and close output root file
  nlo->Write();
  cout << "\nWrote file: " << nlo->GetName() << endl;
  delete nlo;

  return 0;
}
