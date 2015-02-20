#include <iostream>
#include <iomanip>
#include <sstream>
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
#include <TLine.h>
// #include <TText.h>
// #include <TLatex.h>
#include <TAxis.h>

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

string sigma_prt(Double_t sigma, unsigned prec) {
  stringstream ss;
  ss << "#sigma = " << showpoint << setprecision(prec) << sigma
     << noshowpoint << " pb";
  return ss.str();
}

enum class var : char { other, N, pT, mass, y, phi, HT, tau };

// ******************************************************************
typedef propmap<TH1D*,2> hmap_t;
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

  // Number of events
  TH1D *h_N = get<TH1D>(f,"N");
  const Double_t N_events = h_N->GetAt(1);
  const Double_t N_scale  = 1./N_events;
  cout << "Events: " << lround(h_N->GetEntries()) << endl;

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
            if (obj_h->InheritsFrom(TH1D::Class())) {
              TH1D *h = static_cast<TH1D*>(obj_h);

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

  TLine zero;
  zero.SetLineStyle(7);

  TCanvas canv;
  canv.SaveAs((fout+'[').c_str());

  for (const auto& hname : hmap.pvec<1>()) {
    hkey[1] = hname;

    const string h_name = hname->str();

    TLegend leg(0.72,0.92-0.0325*hmap.pvec<0>().size(),0.95,0.92);
    TLegend sig(leg.GetX1(),2*leg.GetY1()-leg.GetY2(),
                leg.GetX2(),leg.GetY1());
    leg.SetFillColorAlpha(0,0.65);
    sig.SetFillColorAlpha(0,0.65);

    TLine top_corner_cover(leg.GetX1(),0.9,0.91,0.9);
    TLine right_corner_cover(0.9,sig.GetY1(),0.9,0.91);
    top_corner_cover.SetNDC();
    right_corner_cover.SetNDC();
    top_corner_cover.SetLineColor(0);
    right_corner_cover.SetLineColor(0);
    top_corner_cover.SetLineWidth(3);
    right_corner_cover.SetLineWidth(3);

    hist_range y_range;

    int i = 0;
    TH1D *first = nullptr;
    for (const auto& hdir : hmap.pvec<0>()) {
      hkey[0] = hdir;

      TH1D *h;
      if (!hmap.get(hkey,h)) {
        cerr << "Warning: No histogram " << h_name
             << " in dir " << hdir->str() << endl;
        continue;
      }

      const Int_t nbins = h->GetNbinsX()+1;

      h->Scale(N_scale);

      // Cross section
      const Double_t sigma = ( h_name.find("_N_incl") != string::npos
                           ? h->GetAt(1)
                           : h->Integral(0,nbins) );

      for (Int_t i=0; i<nbins; ++i)
        h->SetAt(h->GetAt(i)/(h->GetBinWidth(i)),i);

      y_range(h->GetMinimum(),h->GetMaximum());

      if (i==0) {
        first = h;
        h->Draw();
      } else h->Draw("same");

      h->SetLineWidth(2);
      h->SetLineColor(color[i]);
      h->SetMarkerColor(color[i]);
      // h->SetFillColorAlpha(0,0);
      h->Sumw2(false);
      leg.AddEntry(h,lmap[hdir].c_str());
      sig.AddEntry(h,sigma_prt(sigma,3).c_str());

      ++i;
    }

    y_range(first)->SetTitle(hname->str().c_str());

    if (first->GetMinimum()<0.)
      zero.DrawLine(first->GetBinLowEdge(1),0,
                    first->GetBinLowEdge(first->GetNbinsX()+1),0);

    // Variable type
    var h_var = var::other;

    if (h_name.find("_N") != string::npos) {
      h_var = var::N;
    } else if (h_name.find("_pT") != string::npos) {
      h_var = var::pT;
    } else if (h_name.find("_mass") != string::npos) {
      h_var = var::mass;
    } else if ( (h_name.find("_deltay") != string::npos) ||
                (h_name.find("_y") != string::npos) ) {
      h_var = var::y;
    } else if ( (h_name.find("_deltaphi") != string::npos) ||
                (h_name.find("_phi") != string::npos) ) {
      h_var = var::phi;
    } else if (h_name.find("_HT") != string::npos) {
      h_var = var::HT;
    } else if (h_name.find("_tau") != string::npos) {
      h_var = var::tau;
    }

    TAxis *xa = first->GetXaxis();
    TAxis *ya = first->GetYaxis();
    ya->SetTitleOffset(1.3);
    switch (h_var) {
      case var::pT:
        xa->SetTitle("pT, GeV");
        ya->SetTitle("d#sigma/dp_{T}, pb/GeV");
        break;
      case var::HT:
        xa->SetTitle("HT, GeV");
        ya->SetTitle("d#sigma/dH_{T}, pb/GeV");
        break;
      case var::mass:
        xa->SetTitle("m, GeV");
        ya->SetTitle("d#sigma/dm, pb/GeV");
        break;
      case var::y:
        xa->SetTitle("y");
        ya->SetTitle("d#sigma/dy, pb");
        break;
      case var::phi:
        xa->SetTitle("#phi, rad");
        ya->SetTitle("d#sigma/d#phi, pb/rad");
        break;
      case var::tau:
        xa->SetTitle("#tau, GeV");
        ya->SetTitle("d#sigma/d#tau, pb/GeV");
        break;
      default:
        ya->SetTitle("#sigma, pb");
    }

    top_corner_cover.Draw();
    right_corner_cover.Draw();
    leg.Draw();
    sig.Draw();
    canv.SaveAs(fout.c_str());
  }

  canv.SaveAs((fout+']').c_str());

  delete f;

  return 0;
}
