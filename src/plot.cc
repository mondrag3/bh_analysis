#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
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
#include <TLatex.h>
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

hmap_t hmap;
hkey_t hkey;

vector<string> scales;

Double_t N_scale = 0.;
Double_t N_total = 0.;

void hist_key(TH1* h) noexcept {
  hkey[0] = new prop<string>(h->GetName());
}

bool dir_key(const TDirectory* d) noexcept {
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

void get_hists(const TDirectory *d, bool first) noexcept {
  static TKey *key;
  TIter nextkey(d->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    static TObject *obj;
    obj = key->ReadObj();
    if (obj->InheritsFrom(TH1::Class())) {
      static TH1* h;
      h = static_cast<TH1*>(obj);
      h->Scale(N_scale);
      hist_key(h);
      if (first) hmap.insert(hkey,h);
      else hmap.get(hkey)->Add(h);
    }
  }
}

string sigma_prt(Double_t sigma, unsigned prec) {
  stringstream ss;
  ss << "#sigma = " << showpoint << setprecision(prec) << sigma
     << noshowpoint << " pb";
  return ss.str();
}

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> fin;
  string fout, jet_alg, pdf, title;
  float of_lim;
  bool  of_draw = false;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<vector<string>>(&fin)->required()->multitoken(),
     "input file with histograms")
    ("output,o", po::value<string>(&fout)->required(),
     "output pdf plots")
    ("jet-alg,j", po::value<string>(&jet_alg)->default_value("AntiKt4"),
     "jet algorithm")
    ("pdf", po::value<string>(&pdf),
     "select PDF set if there are multiple in the file")
    ("scales", po::value<vector<string>>(&scales)
     ->default_value({"0.25Ht","0.5Ht","1Ht"},"Ht/4, Ht/2, Ht"),
     "fac and ren scales")
    ("title", po::value<string>(&title),
     "string appended to each title")
    ("overflow", po::value<float>(&of_lim),
     "print under- and overflow messages on histograms above the limit")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, all_opt), vm);
    if (argc == 1 || vm.count("help")) {
      cout << all_opt << endl;
      return 0;
    }
    po::notify(vm);
    if (vm.count("overflow")) of_draw = true;
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Get histograms *************************************************
  TFile *first_file = nullptr;
  for (const auto& fname : fin) {
    static bool first = true;

    TFile *f = new TFile(fname.c_str(),"read");
    if (f->IsZombie()) exit(1);

    const Double_t N = ((TH1*)f->Get("N"))->GetBinContent(1);
    N_scale = 1./N;
    N_total += N;

    cout << "\033[36mFile  :\033[0m " << f->GetName() << endl;
    cout << "\033[36mEvents:\033[0m " << N << endl << endl;

    static TKey *key1;
    TIter nextkey1(f->GetListOfKeys());
    while ((key1 = (TKey*)nextkey1())) {
      static TObject *obj1;
      obj1 = key1->ReadObj();
      if (obj1->InheritsFrom(TDirectory::Class())) {

        const TDirectory *d1 = static_cast<TDirectory*>(obj1);

        if (dir_key(d1)) get_hists(d1,first);

        else if (!jet_alg.compare(d1->GetName())) {
          static TKey *key2;
          TIter nextkey2(d1->GetListOfKeys());
          while ((key2 = (TKey*)nextkey2())) {
            static TObject *obj2;
            obj2 = key2->ReadObj();
            if (obj2->InheritsFrom(TDirectory::Class())) {

              const TDirectory *d2 = static_cast<TDirectory*>(obj2);

              if (dir_key(d2)) get_hists(d2,first);

            }
          }
        }
      }
    }
    if (first) first_file = f;
    else delete f;
    first = false;
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

  cout << "\033[36mPDF    :\033[0m " << hkey[3] << endl << endl;
  cout << "\033[36mJet Alg:\033[0m " << jet_alg << endl << endl;

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
      pdf_lo[i]    = (cent[i] - h_pdf_lo->GetBinContent(i+1));
      pdf_hi[i]    = (h_pdf_hi->GetBinContent(i+1) - cent[i]);
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

    const string _hname( hname->str() );

    const Double_t sigma   = ( _hname.find("NJet_incl") == string::npos
                               ? h_cent->Integral(0,nbins+1)
                               : h_cent->GetBinContent(1) );
    const Double_t sigma_u = h_cent->GetBinContent(0);
    const Double_t sigma_o = h_cent->GetBinContent(nbins+1);

    const bool is_pT   = (_hname.find("_pT")       != string::npos);
    const bool is_HT   = (_hname.find("_HT")       != string::npos);
    const bool is_mass = (_hname.find("_mass")     != string::npos);
    const bool is_eta  = (_hname.find("_deltay")   != string::npos) ||
                         (_hname.find("_y")        != string::npos);
    const bool is_phi  = (_hname.find("_deltaphi") != string::npos) ||
                         (_hname.find("_phi")      != string::npos);
    const bool is_tau  = (_hname.find("_tau")      != string::npos);

    if (is_pT||is_HT||is_mass||is_eta||is_phi||is_tau) {
      Float_t width;
      for (size_t i=0;i<nbins;++i) {
        width = bins_wdth[i];
        cent     [i] /= width;
        h_cent->SetBinContent(i+1,cent[i]);
        pdf_lo   [i] /= width;
        pdf_hi   [i] /= width;
        scales_hi[i] /= width;
        scales_lo[i] /= width;
      }
    }

    TGraphAsymmErrors g_scales (nbins,bins_edge.data(),cent.data(),
                                      0,bins_wdth.data(),
                                      scales_lo.data(),scales_hi.data()),
                      g_pdf_unc(nbins,bins_edge.data(),cent.data(),
                                      0,bins_wdth.data(),
                                      pdf_lo.data(),pdf_hi.data());

    g_scales.GetXaxis()
      ->SetRangeUser(bins_edge[0],bins_edge.back()+bins_wdth.back());
    g_scales.SetTitle(
      ( title.size() ? (h_cent->GetName()+(' '+title)).c_str()
                     :  h_cent->GetName() )
    );
    TAxis *xa = g_scales.GetXaxis();
    TAxis *ya = g_scales.GetYaxis();
    ya->SetTitleOffset(1.3);
    if (is_pT) {
      xa->SetTitle("pT, GeV");
      ya->SetTitle("d#sigma/dp_{T}, pb/GeV");
    } else if (is_HT) {
      xa->SetTitle("HT, GeV");
      ya->SetTitle("d#sigma/dHT, pb/GeV");
    } else if (is_mass) {
      xa->SetTitle("m, GeV");
      ya->SetTitle("d#sigma/dm, pb/GeV");
    } else if (is_eta) {
      xa->SetTitle("#eta");
      ya->SetTitle("d#sigma/d#eta, pb");
    } else if (is_phi) {
      xa->SetTitle("#phi, rad");
      ya->SetTitle("d#sigma/d#phi, pb/rad");
    } else if (is_tau) {
      xa->SetTitle("#tau, GeV");
      ya->SetTitle("d#sigma/d#tau, pb/GeV");
    } else {
      ya->SetTitle("#sigma, pb");
    }

    g_scales.SetFillColorAlpha(2,0.5);
    // g_scales .SetLineColor(10);
    // g_scales .SetFillStyle(3004);
    g_scales.SetLineWidth(0);
    // g_scales .SetMarkerColor(4);
    // g_scales .SetMarkerStyle(21);
    g_scales.Draw("a2");

    g_pdf_unc.SetFillColorAlpha(4,0.5);
    // g_pdf_unc.SetLineColor(10);
    // g_pdf_unc.SetFillStyle(3005);
    g_pdf_unc.SetLineWidth(0);
    g_pdf_unc.Draw("2");

    h_cent->Sumw2(false);
    h_cent->SetLineWidth(2);
    h_cent->SetLineColor(1);
    h_cent->Draw("same");
    // h_pdf_lo->Draw("same");
    // h_pdf_hi->Draw("same");

    TLegend leg(0.72,0.75,0.89,0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(&g_scales, "Scale unc");
    leg.AddEntry(&g_pdf_unc,"PDF unc");
    leg.Draw();

    TLatex cs_lbl(0.73,0.74, sigma_prt(sigma,3).c_str());
    cs_lbl.SetNDC();
    cs_lbl.SetTextAlign(13);
    cs_lbl.SetTextFont(42);
    cs_lbl.SetTextSize(0.035);
    cs_lbl.Draw();

    TLatex N_lbl(0.73,0.70, Form("Events: %.2e",N_total));
    N_lbl.SetNDC();
    N_lbl.SetTextAlign(13);
    N_lbl.SetTextFont(42);
    N_lbl.SetTextSize(0.035);
    N_lbl.Draw();

    TText jet_alg_lbl(0.73,0.66,("Jet alg: "+jet_alg).c_str());
    jet_alg_lbl.SetNDC();
    jet_alg_lbl.SetTextAlign(13);
    jet_alg_lbl.SetTextFont(42);
    jet_alg_lbl.SetTextSize(0.035);
    jet_alg_lbl.Draw();

    TLatex pdf_lbl(0.73,0.62, ("PDF: "+hkey[3]->str()).c_str());
    pdf_lbl.SetNDC();
    pdf_lbl.SetTextAlign(13);
    pdf_lbl.SetTextFont(42);
    pdf_lbl.SetTextSize(0.035);
    pdf_lbl.Draw();

    if (of_draw) {
      if (sigma_u/sigma > of_lim) {
        stringstream ss;
        ss << setprecision(2) << 100.*sigma_u/sigma << "% underflow";
        cout << "\033[33m" << ss.str() << "\033[0m" << endl;
        TText *u_lbl = new TText(0.3,0.45,ss.str().c_str());
        u_lbl->SetNDC();
        u_lbl->SetTextAlign(13);
        u_lbl->SetTextColor(2);
        u_lbl->SetTextSize(0.1);
        u_lbl->Draw();
      }
      if (sigma_o/sigma > of_lim) {
        stringstream ss;
        ss << setprecision(2) << 100.*sigma_o/sigma << "% overflow";
        cout << "\033[33m" << ss.str() << "\033[0m" << endl;
        TText *o_lbl = new TText(0.3,0.55,ss.str().c_str());
        o_lbl->SetNDC();
        o_lbl->SetTextAlign(13);
        o_lbl->SetTextColor(2);
        o_lbl->SetTextSize(0.1);
        o_lbl->Draw();
      }
    }

    canv.SaveAs(fout.c_str());

  }

  canv.SaveAs((fout+']').c_str());

  delete first_file;

  return 0;
}
