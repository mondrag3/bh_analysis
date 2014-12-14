#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>

#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>

#include <kiwi/propmap.h>

using namespace std;
namespace po = boost::program_options;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

// ******************************************************************
void hist_key(vector<prop_ptr>& hkey, TH1* h) noexcept {
  hkey[0] = new prop<string>(h->GetName());
}

bool dir_key(vector<prop_ptr>& hkey, TDirectory* d) noexcept {
  static const boost::regex regex("Fac(.*)_Ren(.*)_PDF(.*)_(.*)");
  static boost::smatch result;
  if ( boost::regex_search(string(d->GetName()), result, regex) ) {
    hkey[1] = new prop<string>( string(result[1].first, result[1].second) );
    hkey[2] = new prop<string>( string(result[2].first, result[2].second) );
    hkey[3] = new prop<string>( string(result[3].first, result[3].second) );
    hkey[4] = new prop<string>( string(result[4].first, result[4].second) );
    return true;
  } else return false;
}

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  string ifile, ofile, jet_alg;
  vector<string> scales;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<string>(&ifile)->required(),
     "input file with histograms")
    ("output,o", po::value<string>(&ofile)->required(),
     "output pdf plots")
    ("jet-alg,j", po::value<string>(&jet_alg)->default_value("AntiKt4"),
     "jet algorithm")
    ("scales", po::value<vector<string>>(&scales)
     ->default_value({"0.25Ht","0.5Ht","1Ht"},"Ht/4, Ht/2, Ht"),
     "fac and ren scales")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, all_opt), vm);
    if (argc == 1 || vm.count("help")) {
      cout << all_opt << endl;
      return 0;
    }
    po::notify(vm);
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  TFile *f = new TFile(ifile.c_str(),"read");
  if (f->IsZombie()) exit(1);

  // property map of histograms
  propmap<TH1*> hists(5);
  vector<prop_ptr> hkey(5);

  // Get histograms *************************************************
/*  {
    static TKey *key1;
    TIter nextkey1(f->GetListOfKeys());
    while ((key1 = (TKey*)nextkey1())) {
      static TObject *obj1;
      obj1 = key1->ReadObj();
      if (obj1->InheritsFrom(TDirectory::Class())) {

        TDirectory *d1 = static_cast<TDirectory*>(obj1);
        static TKey *key2;
        TIter nextkey2(d1->GetListOfKeys());
        while ((key2 = (TKey*)nextkey2())) {
          static TObject *obj2;
          obj2 = key2->ReadObj();
          if (obj2->InheritsFrom(TDirectory::Class())) {

            TDirectory *d2 = static_cast<TDirectory*>(obj2);
            static TKey *key3;
            TIter nextkey3(d2->GetListOfKeys());
            while ((key3 = (TKey*)nextkey3())) {
              static TObject *obj3;
              obj3 = key3->ReadObj();

              if (obj3->InheritsFrom(TH1::Class())) {
                TH1* h = static_cast<TH1*>(obj3);
                hist_key(hkey,h,d2,d1);
                hists.insert(hkey,h);
              }
            }

          } else if (obj2->InheritsFrom(TH1::Class())) {
            TH1* h = static_cast<TH1*>(obj2);
            hist_key(hkey,h,d1);
            hists.insert(hkey,h);
          }
        }

      } else if (obj1->InheritsFrom(TH1::Class())) {
        TH1* h = static_cast<TH1*>(obj1);
        hist_key(hkey,h);
        hists.insert(hkey,h);
      }
    }
  }*/

  {
    static TKey *key1;
    TIter nextkey1(f->GetListOfKeys());
    while ((key1 = (TKey*)nextkey1())) {
      static TObject *obj1;
      obj1 = key1->ReadObj();
      if (obj1->InheritsFrom(TDirectory::Class())) {

        TDirectory *d1 = static_cast<TDirectory*>(obj1);

        if (dir_key(hkey,d1)) {
          static TKey *key2;
          TIter nextkey2(d1->GetListOfKeys());
          while ((key2 = (TKey*)nextkey2())) {
            static TObject *obj2;
            obj2 = key2->ReadObj();
            if (obj2->InheritsFrom(TH1::Class())) {
              TH1* h = static_cast<TH1*>(obj2);
              hist_key(hkey,h);
              hists.insert(hkey,h);
            }
          }
        } else if (!jet_alg.compare(d1->GetName())) {
          static TKey *key2;
          TIter nextkey2(d1->GetListOfKeys());
          while ((key2 = (TKey*)nextkey2())) {
            static TObject *obj2;
            obj2 = key2->ReadObj();
            if (obj2->InheritsFrom(TDirectory::Class())) {

              TDirectory *d2 = static_cast<TDirectory*>(obj2);

              if (dir_key(hkey,d2)) {
                static TKey *key3;
                TIter nextkey3(d2->GetListOfKeys());
                while ((key3 = (TKey*)nextkey3())) {
                  static TObject *obj3;
                  obj3 = key3->ReadObj();
                  if (obj3->InheritsFrom(TH1::Class())) {
                    TH1* h = static_cast<TH1*>(obj3);
                    hist_key(hkey,h);
                    hists.insert(hkey,h);
                  }
                }
              }

            }
          }
        }
      }
    }
  }

  // Get histograms *************************************************


  return 0;
}
