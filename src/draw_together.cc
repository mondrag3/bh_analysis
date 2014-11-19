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
      ("input,i", po::value< vector<string> >(&fin),
       "label:file, input root files with histograms")
      ("output,o", po::value<string>(&fout),
       "output pdf plots file")
      ("norm,n", po::bool_switch(&norm),
       "normalize histograms to unity")
      ("logx", po::bool_switch(&logx),
       "")
      ("logy", po::bool_switch(&logy),
       "")
      ("name-title", po::bool_switch(&name_title),
       "replace hist title with name")
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

  for (auto it=h_.begin(), end=h_.end(); it!=end; ++it) {
    const size_t nh = it->second.size();
    TH1* h = it->second[0];
    if (norm) h->Scale(1./h->Integral());
    double ymax = h->GetMaximum();

    for (size_t i=1;i<nh;++i) {
      h = it->second[i];
      if (norm) h->Scale(1./h->Integral());
      double max = h->GetMaximum();
      if (max>ymax) ymax = max;
    }

    it->second[0]->SetAxisRange(0.,ymax*1.05,"Y");
  }

  TCanvas canv;
  gStyle->SetOptStat(0);
  if (logx) canv.SetLogx();
  if (logy) canv.SetLogy();

  canv.SaveAs((fout+'[').c_str());

  for (auto it=h_.begin(), end=h_.end(); it!=end; ++it) {
    auto& hv = it->second;
    const size_t nh = hv.size();

    TLegend leg(0.75,0.92-0.06*nh,0.95,0.92);
    leg.SetFillColor(0);

    for (size_t i=0;i<nh;++i) {
      TH1* h = hv[i];
      leg.AddEntry(h,Form("%s N=%.0f",
        fin[i].substr(0,fin[i].find(':')).c_str(),
        h->GetEntries()
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
