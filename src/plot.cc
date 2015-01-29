#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>

using namespace std;
namespace po = boost::program_options;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

template<class T>
inline T* get(TDirectory* dir, const char* name) {
  T *obj = nullptr;
  dir->GetObject(name,obj);
  if (obj) return obj;
  else {
    cerr << "No object " << name << " in " << dir->GetName() << endl;
    exit(1);
  }
}

template<class T>
void get_all(TDirectory* dir, vector<T*>& hists) noexcept {
  TIter nextkey(dir->GetListOfKeys());
  TKey *key;
  while ((key = static_cast<TKey*>(nextkey()))) {
    TObject *obj = key->ReadObj();
    if (obj->InheritsFrom(T::Class()))
      hists.push_back( static_cast<T*>(obj) );
  }
}

enum class var : char { other, N, pT, mass, y, phi, HT, tau };

int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  string fin_name, title;
  string fout_name;
  float of_lim;
  bool of_lim_set = false;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<string>(&fin_name)->required(),
     "*input files with histograms")
    ("output,o", po::value<string>(&fout_name),
     "output pdf plots")
    ("title,t", po::value<string>(&title),
     "string appended to each title")
    ("overflow", po::value<float>(&of_lim),
     "print under- and overflow messages on histograms above the limit")
    ;
    
    po::positional_options_description pos;
    pos.add("input",1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(all_opt).positional(pos).run(), vm);
    if (argc == 1 || vm.count("help")) {
      cout << all_opt << endl;
      return 0;
    }
    po::notify(vm);
    
    if (!fout_name.size()) {
      size_t first = fin_name.rfind('/')+1;
      if (first==string::npos) first = 0;
      fout_name = fin_name.substr(first,fin_name.rfind('.')-first)+".pdf";
    }
    if (vm.count("overflow")) of_lim_set = true;
      
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************
  
  // Open files
  TFile *fin = new TFile(fin_name.c_str(),"read");
  if (fin->IsZombie()) exit(1);
  cout << "Input file: " << fin->GetName() << endl;
  cout << "Output file: " << fout_name << endl;
  cout << endl;
  
  // Number of events
  const Double_t N = get<TH1>(f,"N")->GetBinContent(1);
  const Double_t N_scale = 1./N;
  cout << "Events: " << N << endl;

  // Collect histograms
  vector<TDirectory*> dirs;
  dirs.reserve(16);
  get_all(fin,dirs);
  for (auto d : dirs) cout << d->GetName() << endl;
  cout << endl;
  
  vector<TH1*> central;
  central.reserve(32);
  get_all(dirs[0],central);
  
  // Draw histograms
  TCanvas canv;
  canv.SaveAs((fout_name+'[').c_str());
  
  for (auto h : central) {
    const string h_name(h->GetName());
    cout << h_name << endl;
    
    var h_var = var::other;
    
    if (_hname.find("_pT") != string::npos) {
      h_var = pT;
    } else if (_hname.find("_mass") != string::npos){
    
    } else if (_hname.find("_deltay") != string::npos){
    
    } else if (_hname.find("_deltaphi") != string::npos){
    
    } else if (_hname.find("_HT") != string::npos){
    
    } else if (_hname.find("_tau") != string::npos){
    
    } else if (_hname.find("_pT") != string::npos){
    
    }
    
    const size_t nbins = h->GetNbinsX();
    vector<double_t> bins_edge(nbins,0.),
                     bins_wdth(nbins,0.),
                     cent     (nbins,0.),
                     scales_lo(nbins,0.),
                     scales_hi(nbins,0.),
                     pdf_lo   (nbins,0.),
                     pdf_hi   (nbins,0.);

    const Double_t sigma   = ( h_name.find("jets_N_incl") == string::npos
                           ? h->Integral(0,nbins+1)
                           : h->GetBinContent(1) );
    const Double_t sigma_u = h->GetBinContent(0);
    const Double_t sigma_o = h->GetBinContent(nbins+1);
    
    h_cent->Sumw2(false);
    h_cent->SetLineWidth(2);
    h_cent->SetLineColor(1);
    h_cent->Draw("same");
    
    /*
    for (size_t di=1,dn=dirs.size();di<dn;++di) {
      string dir_name(dirs[di]->GetName());
      if (dir_name.substr(dir_name.rfind('-')+1)=="up") { // pdf up
        
      } else if (dir_name.substr(dir_name.rfind('-')+1)=="down") { // pdf up
        
      } else { // scale variation
        
      }
    }
    */
    
    canv.SaveAs(fout_name.c_str());
  }

  canv.SaveAs((fout_name+']').c_str());
  
  delete fin;

  return 0;
}
