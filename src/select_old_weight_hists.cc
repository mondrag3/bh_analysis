#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>

#include "propmap.h"

using namespace std;
namespace po = boost::program_options;

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

class scale_p: public prop<float> {
  float init(const string& s) {
    const size_t sep = s.find("Ht");
    float x=1.;
    if (sep==string::npos) return -1;
    if (sep>0) x *= atof(s.substr(0,sep).c_str());
    if (s.size()-sep>2) x /= atof(s.substr(sep+2).c_str());
    return x;
  }
public:
  scale_p(const string& s): prop<float>(init(s)) { }
};

bool parse_h_weight(const string& name, vector<prop_ptr>& key) {
  static const boost::regex regex(
    "(.*)_(.*)_PDF(.*)_SCALEFac(.*)Ren(.*)_CUTS(.*)");
  static boost::smatch result;
  if ( boost::regex_search(name, result, regex) ) {
    key[0] = new prop<string>( string(result[1].first, result[1].second) );
    key[1] = new prop<string>( string(result[2].first, result[2].second) );
    key[2] = new prop<string>( string(result[3].first, result[3].second) );
    key[3] = new scale_p     ( string(result[4].first, result[4].second) );
    key[4] = new scale_p     ( string(result[5].first, result[5].second) );
    key[5] = new prop<string>( string(result[6].first, result[6].second) );

    return true;
  } else return false;
}

bool parse_h_weight_default(const string& name, vector<prop_ptr>& key) {
  static const boost::regex regex(
    "(.*)_(.*)_PDF(.*)_SCALE.*_CUTS(.*)");
  static boost::smatch result;
  if ( boost::regex_search(name, result, regex) ) {
    key[0] = new prop<string>( string(result[1].first, result[1].second) );
    key[1] = new prop<string>( string(result[2].first, result[2].second) );
    key[2] = new prop<string>( string(result[3].first, result[3].second) );
    key[3] = new prop<float> (0);
    key[4] = new prop<float> (0);
    key[5] = new prop<string>( string(result[4].first, result[4].second) );

    return true;
  } else return false;
}

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  string ifile, ofile, pdf_name;
  bool uncut;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
      ("help,h", "produce help message")
      ("input,i", po::value<string>(&ifile),
       "input file with all histograms")
      ("output,o", po::value<string>(&ofile),
       "output file with only weight histograms")
      ("pdf", po::value<string>(&pdf_name)->default_value("MSTW2008nlo68cl"),
       "select pdf name")
      ("uncut", po::bool_switch(&uncut),
       "select full weights histograms instead of those in \"Jets\" directory")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, all_opt), vm);
    po::notify(vm);

    // Options Properties ---------------------------------
    if (argc == 1 || vm.count("help")) {
      cout << all_opt << endl;
      return 0;
    }

    // Necessary options ----------------------------------
    vector<string> rec_opt;
    rec_opt.push_back("input");
    rec_opt.push_back("output");

    for (size_t i=0, size=rec_opt.size(); i<size; ++i) {
      if (!vm.count(rec_opt[i]))
      { cerr << "Missing command --" << rec_opt[i] << endl; return 1; }
    }
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    return 1;
  }
  // END OPTIONS ****************************************************

  propmap<TH1*> _h_old(6);
  vector<prop_ptr> hkey(6);

  // Open input file
  TFile *fin = new TFile(ifile.c_str(),"read");
  if (fin->IsZombie()) return 1;

  // Get histograms from the "old" file
  if (uncut) {

    TKey *key1;
    TIter nextkey1(fin->GetListOfKeys());
    while ((key1 = (TKey*)nextkey1())) {
      TObject *obj = key1->ReadObj();

      if (obj->IsFolder()) {
        TDirectory* dir = dynamic_cast<TDirectory*>(obj);
        if ( !strcmp(dir->GetName(),"Jets") ) continue;

        TKey *key2;
        TIter nextkey2(dir->GetListOfKeys());
        while ((key2 = (TKey*)nextkey2())) {
          obj = key2->ReadObj();

          string histname(obj->GetName());
          const size_t sep = histname.find("__");
          if ( !histname.substr(0,sep).compare("h_weight_all") ) {
            histname = "all_"+histname.substr(sep+2);
            if ( parse_h_weight(histname,hkey) ||
                 parse_h_weight_default(histname,hkey) )
            {
              _h_old.insert(hkey,dynamic_cast<TH1*>(obj));
            }
          }
        }

      } // end if folder
    } // end while key

  } else {

    TKey *key1;
    TIter nextkey1(dynamic_cast<TDirectory*>(fin->Get("Jets"))->GetListOfKeys());
    while ((key1 = (TKey*)nextkey1())) {
      TObject *obj = key1->ReadObj();

      if (obj->IsFolder()) {
        TDirectory* dir = dynamic_cast<TDirectory*>(obj);

        TKey *key2;
        TIter nextkey2(dir->GetListOfKeys());
        while ((key2 = (TKey*)nextkey2())) {
          obj = key2->ReadObj();

          string histname(obj->GetName());
          const size_t sep = histname.find("__");
          if ( !histname.substr(0,sep).compare("h_weight") ) {
            histname = histname.substr(sep+2);
            if ( parse_h_weight(histname,hkey) ||
                 parse_h_weight_default(histname,hkey) )
            {
              _h_old.insert(hkey,dynamic_cast<TH1*>(obj));
            }
          }
        }

      } // end if folder
    } // end while key

  }

  // Open output file
  TFile *fout = new TFile(ofile.c_str(),"recreate");
  if (fout->IsZombie()) return 1;

  hkey[2] = new prop<string>(pdf_name);
  hkey[1] = new prop<string>("Elec");
  pmloop(_h_old,fac,3) {
    hkey[3] = *fac;
    pmloop(_h_old,ren,4) {
      hkey[4] = *ren;

      TH1* h;
      if ( _h_old.get(hkey,h) ) {
        stringstream ss;
        /*ss << hkey[1]->str() << '_'
           << hkey[2]->str() << "_Fac"
           << hkey[3]->str() << "_Ren"
           << hkey[4]->str();*/
        ss << "Fac" << hkey[3]->str() << "Ht_"
           << "Ren" << hkey[4]->str() << "Ht_"
           << "PDF" << hkey[2]->str() << "_cent";
        h->Clone(ss.str().c_str());
        cout << h->GetName() << endl;
      }

    }
  }

  fout->Write();
  fout->Close();
  fin->Close();
  delete fin;
  delete fout;

  return 0;
}
