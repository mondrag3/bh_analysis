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
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>

#include <kiwi/propmap11.h>

using namespace std;
namespace po = boost::program_options;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

// ******************************************************************
typedef propmap<TH1*,5> hmap_t;
typedef hmap_t::Key hkey_t;

void hist_key(hkey_t& hkey, TH1* h) noexcept {
  hkey[0] = new prop<string>(h->GetName());
}

bool dir_key(hkey_t& hkey, TDirectory* d) noexcept {
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
  string fin, fout, jet_alg;
  vector<string> scales;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<string>(&fin)->required(),
     "input file with histograms")
    ("output,o", po::value<string>(&fout)->required(),
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

  TFile *f = new TFile(fin.c_str(),"read");
  if (f->IsZombie()) exit(1);

  // property map of histograms
  hmap_t hmap;
  hkey_t hkey;

  // Get histograms *************************************************
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
              hmap.insert(hkey,h);
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
                    hmap.insert(hkey,h);
                  }
                }
              }

            }
          }
        }
      }
    }
  }

  // Make plots *****************************************************
  cout << "\033[36mPDF:\033[0m " << hkey[3]->str() << endl << endl;

  const prop_ptr pdf_cent = new prop<string>("cent");
  const prop_ptr pdf_down = new prop<string>("down");
  const prop_ptr pdf_up   = new prop<string>("up");

  TCanvas canv;
  canv.SaveAs((fout+'[').c_str());

  TH1* h;
  hmap.loop<0>(hkey,[&](hkey_t key) noexcept {
    cout << "\033[36mHistogram:\033[0m " << key[0]->str() << endl;
    key[4] = pdf_cent;
    int i=0;
    hmap.loop<1>(hkey,[&](hkey_t key) noexcept {
      hmap.loop<2>(hkey,[&](hkey_t key) noexcept {
        hmap.loop<3>(hkey,[&](hkey_t key) noexcept {

          if (hmap.get(key,h)) {
            h->SetLineColor(i+2);
            h->SetTitle(key[0]->str().c_str());

            if (i==0) h->Draw();
            else h->Draw("same");
            ++i;
          }

        });
      });
    });

    // leg.Draw();
    canv.SaveAs(fout.c_str());

  });

  canv.SaveAs((fout+']').c_str());

  delete f;

  return 0;
}
