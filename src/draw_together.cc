#include <iostream>
#include <vector>
#include <unordered_map>
#include <memory>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TH1.h>

using namespace std;
namespace po = boost::program_options;

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

pair<Double_t,Double_t> TH1Range(const TH1* h, int minbin=1) noexcept {
  Double_t min = 0.;
  int i = minbin,
      n = h->GetNbinsX()+1;
  for (;i<n;++i) {
    Double_t min = h->GetBinContent(i);
    if (min>0.) break;
  }
  Double_t max = min;
  for (;i<n;++i) {
    Double_t x = h->GetBinContent(i);
    if (x<min) min = x;
    if (x>max) max = x;
  }
  return make_pair(min,max);
}
pair<Double_t,Double_t> TH1PositiveRange(const TH1* h, int minbin=1) noexcept {
  Double_t min = 0.;
  int i = minbin,
      n = h->GetNbinsX()+1;
  for (;i<n;++i) {
    min = h->GetBinContent(i);
    if (min>0.) break;
  }
  Double_t max = min;
  for (;i<n;++i) {
    Double_t x = h->GetBinContent(i);
    if (x>0.) {
      if (x<min) min = x;
      if (x>max) max = x;
    }
  }
  return make_pair(min,max);
}

/*
inline void TH1SetLogXRange(TH1* h) noexcept {
  h->SetAxisRange(
    h->GetBinLowEdge(h->FindFixBin(0.)+1)/10.,
    h->GetBinLowEdge(h->GetNbinsX()+1),
    "X"
  );
}
*/

int main(int argc, char *argv[])
{
  // START OPTIONS **************************************************
  vector<string> fin;
  string fout;
  bool norm, logx, logy, name_title;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
      ("help,h", "produce help message")
      ("input,i", po::value< vector<string> >(&fin)->required(),
       "label:file, input root files with histograms")
      ("output,o", po::value<string>(&fout)->required(),
       "output pdf plots file")
      ("norm,n", po::bool_switch(&norm),
       "normalize histograms to unity")
      ("logx", po::bool_switch(&logx), "")
      ("logy", po::bool_switch(&logy), "")
      ("name-title", po::bool_switch(&name_title),
       "replace hist title with name")
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
    return 1;
  }
  // END OPTIONS ****************************************************

  vector< unique_ptr<TFile> > f;
  unordered_map< string, vector<TH1*> > h_;

  for (auto& input_name : fin) {
    string fname = input_name.substr(input_name.find(':')+1);
    if (fname.size()==0) {
      cout << "No separator \':\' in argument: " << input_name << endl;
      exit(1);
    }

    f.push_back( unique_ptr<TFile>(new TFile(fname.c_str(),"read")) );
    if ( f.back()->IsZombie() ) exit(1);

    static TKey *key;
    TIter nextkey(f.back()->GetListOfKeys());
    while (( key = (TKey*)nextkey() )) {
      TH1 *hist = dynamic_cast<TH1*>( key->ReadObj() );
      hist->SetLineColor(f.size()+1);
      hist->SetMarkerColor(f.size()+1);
      //hist->SetLineWidth(2);
      h_[hist->GetName()].push_back(hist);
    }
  }

  pair<Double_t,Double_t> (*TH1RangeFcn)(const TH1*,int);
  if (logy) TH1RangeFcn = &TH1PositiveRange;
  else      TH1RangeFcn = &TH1Range;

  for (auto& _h : h_) {
    const size_t nh = _h.second.size();
    TH1* h = _h.second[0];
    if (norm) h->Scale(1./h->Integral());
    //double ymax = h->GetMaximum();
    auto range = TH1RangeFcn(h,(logx ? h->FindFixBin(0.) : 1));

    for (size_t i=1;i<nh;++i) {
      h = _h.second[i];
      if (norm) h->Scale(1./h->Integral());
      auto _range = TH1RangeFcn(h,(logx ? h->FindFixBin(0.) : 1));
      if (_range.first <range.first ) range.first  = _range.first;
      if (_range.second>range.second) range.second = _range.second;
    }

    //if (logx) TH1SetLogXRange(_h.second[0]);

    if (logy) {
      _h.second[0]->SetMinimum(pow(range.first,0.98));
      _h.second[0]->SetMaximum(pow(range.second,1.025));
    } else {
      _h.second[0]->SetMinimum(range.first*0.98);
      _h.second[0]->SetMaximum(range.second*1.025);
    }
  }

  TCanvas canv;
  gStyle->SetOptStat(0);
  if (logx) canv.SetLogx();
  if (logy) canv.SetLogy();

  canv.SaveAs((fout+'[').c_str());

  for (auto& _h : h_) {
    auto& hv = _h.second;
    const size_t nh = hv.size();

    TLegend leg(0.72,0.92-0.06*nh,0.95,0.92);
    leg.SetFillColor(0);

    for (size_t i=0;i<nh;++i) {
      TH1* h = hv[i];
      leg.AddEntry(h,Form("%s N=%.2fe6 #int=%.2fe6",
        fin[i].substr(0,fin[i].find(':')).c_str(),
        h->GetEntries()/1e6,
        h->Integral()/1e6
      ) );
      if (name_title) h->SetTitle(h->GetName());
      if (i==0) h->Draw();
      else h->Draw("same");
    }

    leg.Draw();
    canv.SaveAs(fout.c_str());
  }

  canv.SaveAs((fout+']').c_str());

  return 0;
}
