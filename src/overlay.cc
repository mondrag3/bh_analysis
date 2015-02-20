#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <tuple>
#include <cstring>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
// #include <TText.h>
// #include <TLatex.h>
// #include <TAxis.h>

#include "propmap.hh"
#include "hist_range.hh"

using namespace std;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

unordered_map<prop_ptr,string> rm_eq_fields(const vector<prop_ptr>& ss, const char* sep) {
  const size_t n = ss.size();
  vector<list<string>> tokenized(ss.size());
  multiset<string> tracker;
  for (size_t i=0; i<n; ++i) {
    string str = ss[i]->str();
    char cstr[str.size()];
    strcpy(cstr,str.c_str());
    char * pch;
    pch = strtok (cstr,sep);
    while (pch) {
      tokenized[i].emplace_back(pch);
      tracker.emplace(pch);
      //cout << pch << endl;
      pch = strtok(NULL,sep);
    }
  }
  for (auto& tok : tracker) {
    for (auto& tl : tokenized) {
      if (tracker.count(tok)==n) tl.remove(tok);
    }
  }
  //cout << endl;
  
  unordered_map<prop_ptr,string> _map;
  for (size_t i=0; i<n; ++i) {
    string& _s = _map[ss[i]];
    bool first = true;
    for (auto& s : tokenized[i]) {
      _s += ( first ? s : '_'+s);
      if (first) first = false;
    }
    //cout << _s << endl;
  }
  
  return _map;
};

// ******************************************************************
typedef propmap<TH1*,2> hmap_t;
typedef hmap_t::key_t hkey_t;

constexpr Color_t color[] = { 2, 3, 4, 6, 7, 8, 9, 28, 30, 46 };

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
    fout = fin.substr(first,fin.rfind('.')-first)+".pdf";
  }

  // Collect ********************************************************
  hmap_t hmap;
  hkey_t hkey;

  TFile *f = new TFile(fin.c_str(),"read");

  { // Collect histograms
    TIter nextkey_d(f->GetListOfKeys());
    TKey *key_d;
    while ((key_d = static_cast<TKey*>(nextkey_d()))) {
      TObject *obj_d = key_d->ReadObj();
      if (obj_d->InheritsFrom(TDirectory::Class())) {
        TDirectory *dir = static_cast<TDirectory*>(obj_d);

          TIter nextkey_h(dir->GetListOfKeys());
          TKey *key_h;
          while ((key_h = static_cast<TKey*>(nextkey_h()))) {
            TObject *obj_h = key_h->ReadObj();
            if (obj_h->InheritsFrom(TH1::Class())) {
              TH1 *h = static_cast<TH1*>(obj_h);

              hkey[0] = new prop<string>(dir->GetName());
              hkey[1] = new prop<string>(h->GetName());
              hmap.insert(hkey,h);

            }
          }

      }
    }
  }
  
  // remove fields identical for all
  unordered_map<prop_ptr,string> lmap = rm_eq_fields(hmap.pvec<0>(),"_");
  
  // Draw ***********************************************************
  gStyle->SetOptStat(0);
  
  TCanvas canv;
  canv.SaveAs((fout+'[').c_str());

  for (const auto& hname : hmap.pvec<1>()) {
    hkey[1] = hname;
      
    hist_range y_range;
    
    TLegend leg(0.72,0.92-0.05*hmap.pvec<0>().size(),0.95,0.92);
    leg.SetFillColor(0);

    int i = 0;
    TH1 *first = nullptr;
    for (const auto& hdir : hmap.pvec<0>()) {
      hkey[0] = hdir;
      
      TH1 *h;
      if (!hmap.get(hkey,h)) {
        cerr << "Warning: No histogram " << hname->str()
             << " in dir " << hdir->str() << endl;
        continue;
      }
      
      y_range(h->GetMinimum(),h->GetMaximum());
      
      if (i==0) {
        first = h;
        h->Draw();
      } else h->Draw("same");

      h->SetLineWidth(2);
      h->SetLineColor(color[i]);
      h->SetMarkerColor(color[i]);
      h->Sumw2(false);
      leg.AddEntry(h,lmap[hdir].c_str());
      
      ++i;
    }

    y_range(first)->SetTitle(hname->str().c_str());

    leg.Draw();
    canv.SaveAs(fout.c_str());
  }
  
  canv.SaveAs((fout+']').c_str());
  
  delete f;

  return 0;
}
