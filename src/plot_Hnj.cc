#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdexcept>

#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>

#include "propmap11.h"

using namespace std;
namespace po = boost::program_options;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

// ******************************************************************
typedef propmap<TH1*,5> hmap_t;
typedef hmap_t::Key hkey_t;

vector<string> scales;

void hist_key(hkey_t& hkey, TH1* h) noexcept {
  hkey[0] = new prop<string>(h->GetName());
}

bool dir_key(hkey_t& hkey, TDirectory* d) noexcept {
  static const boost::regex regex("Fac(.*)_Ren(.*)_PDF(.*)_(.*)");
  static boost::smatch result;
  if ( boost::regex_search(string(d->GetName()), result, regex) ) {
    hkey[1] = new prop<string>( string(result[1].first, result[1].second) );
    if ( find(scales.begin(),scales.end(),hkey[1]->str())==scales.end() )
      return false;
    hkey[2] = new prop<string>( string(result[2].first, result[2].second) );
    if ( find(scales.begin(),scales.end(),hkey[2]->str())==scales.end() )
      return false;
    hkey[3] = new prop<string>( string(result[3].first, result[3].second) );
    hkey[4] = new prop<string>( string(result[4].first, result[4].second) );
    return true;
  } else return false;
}

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  string fin, fout, jet_alg, pdf, part;

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
    ("pdf", po::value<string>(&pdf),
     "select PDF set of there are multiple in the file")
    ("scales", po::value<vector<string>>(&scales)
     ->default_value({"0.25Ht","0.5Ht","1Ht"},"Ht/4, Ht/2, Ht"),
     "fac and ren scales")
    ("part", po::value<string>(&part),
     "calculation part; string appended to each title")
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
  if (pdf.size()) {
    hkey[3] = new prop<string>(pdf);
    if (!hmap.pset<3>().count(hkey[3])) {
      cerr << "No PDF property: " << pdf << endl;
      exit(1);
    }
  } else if (hmap.pset<3>().size()>1) {
    cerr << "More then one PDF property:" << endl;
    for (auto& p : hmap.pset<3>()) {
      cerr <<"  "<< p << endl;
    }
    exit(1);
  }
  cout << "\033[36mPDF:\033[0m " << hkey[3] << endl << endl;

  const prop_ptr pdf_cent = new prop<string>("cent");
  const prop_ptr pdf_down = new prop<string>("down");
  const prop_ptr pdf_up   = new prop<string>("up");

  const prop_ptr mu_cent = new prop<string>(scales[1]);

  TCanvas canv;
  canv.SaveAs((fout+'[').c_str());

  for (auto& hname : hmap.pset<0>()) {
    static TH1 *h(nullptr);
    TH1 *h_cent(nullptr), *h_pdf_lo(nullptr), *h_pdf_hi(nullptr);
    vector<TH1*> h_scales;

    cout << "\033[36mHistogram:\033[0m " << hname << endl;
    hkey[0] = hname;

    for (auto& fac : hmap.pset<1>()) {
      hkey[1] = fac;
      for (auto& ren : hmap.pset<2>()) {
        hkey[2] = ren;
        for (auto& dev : hmap.pset<4>()) {
          hkey[4] = dev;

            if (hmap.get(hkey,h)) {
              if (fac==mu_cent && ren==mu_cent) {
                if (dev==pdf_cent)      h_cent   = h;
                else if (dev==pdf_down) h_pdf_lo = h;
                else if (dev==pdf_up)   h_pdf_hi = h;
              } else h_scales.push_back(h);
            }

        }
      }
    }

    const size_t nbins = h_cent->GetNbinsX();
    vector<Float_t> bins_edge(nbins,0.),
                    bins_wdth(nbins,0.),
                    cent     (nbins,0.),
                    scales_lo(nbins,0.),
                    scales_hi(nbins,0.),
                    pdf_lo   (nbins,0.),
                    pdf_hi   (nbins,0.);

    for (size_t i=0;i<nbins;++i) {
      bins_edge[i] = h_cent->GetBinLowEdge(i+1);
      bins_wdth[i] = h_cent->GetBinLowEdge(i+2) - bins_edge[i];
      cent  [i]    = h_cent->GetBinContent(i+1);
      pdf_lo[i]    = cent[i] - h_pdf_lo->GetBinContent(i+1);
      pdf_hi[i]    = h_pdf_hi->GetBinContent(i+1) - cent[i];
      for (TH1 *hs : h_scales) {
        Double_t x = hs->GetBinContent(i+1) - cent[i];
        if (x>0.) {
          if (scales_hi[i]<x) scales_hi[i] = x;
        } else {
          x = -x;
          if (scales_lo[i]<x) scales_lo[i] = x;
        }
      }
    }

    TGraphAsymmErrors g_scales (nbins,bins_edge.data(),cent.data(),
                                      0,bins_wdth.data(),
                                      scales_lo.data(),scales_hi.data()),
                      g_pdf_unc(nbins,bins_edge.data(),cent.data(),
                                      0,bins_wdth.data(),
                                      pdf_lo.data(),pdf_hi.data());

    g_scales .GetXaxis()
      ->SetRangeUser(bins_edge[0],bins_edge.back()+bins_wdth.back());
    g_scales .SetTitle(
      ( part.size() ? (h_cent->GetName()+(' '+part)).c_str()
                    :  h_cent->GetName() )
    );
    g_scales .SetFillColorAlpha(2,0.5);
    // g_scales .SetLineColor(10);
    // g_scales .SetFillStyle(3004);
    g_scales .SetLineWidth(0);
    // g_scales .SetMarkerColor(4);
    // g_scales .SetMarkerStyle(21);
    g_scales .Draw("a2");
    g_pdf_unc.SetFillColorAlpha(4,0.5);
    // g_pdf_unc.SetLineColor(10);
    // g_pdf_unc.SetFillStyle(3005);
    g_pdf_unc.SetLineWidth(0);
    g_pdf_unc.Draw("2");
    h_cent  ->SetLineWidth(2);
    h_cent  ->SetLineColor(1);
    h_cent  ->Draw("same");
    // h_pdf_lo->Draw("same");
    // h_pdf_hi->Draw("same");

    TLegend leg(0.72,0.75,0.89,0.89);
    leg.SetBorderSize(0);
    leg.AddEntry(&g_scales, "Scale unc");
    leg.AddEntry(&g_pdf_unc,"PDF unc");
    leg.Draw();

    TText jet_alg_lbl(0.73,0.74,("Jet alg: "+jet_alg).c_str());
    jet_alg_lbl.SetNDC();
    jet_alg_lbl.SetTextAlign(13);
    jet_alg_lbl.SetTextFont(42);
    jet_alg_lbl.SetTextSize(0.035);
    jet_alg_lbl.Draw();

    // leg.Draw();
    canv.SaveAs(fout.c_str());

  }

  canv.SaveAs((fout+']').c_str());

  delete f;

  return 0;
}
