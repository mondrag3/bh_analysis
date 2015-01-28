#include <iostream>
#include <string>
#include <sstream>
#include <vector>

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

using namespace std;
namespace po = boost::program_options;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> fin, labels;
  string fout, jet_alg, pdf, title;
  float of_lim;
  bool  of_draw = false;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<vector<string>>(&fin)->required()->multitoken(),
     "*input files with histograms")
    ("output,o", po::value<string>(&fout)->required(),
     "*output pdf plots")
    ("title,t", po::value<string>(&title),
     "string appended to each title")
    ("label,l", po::value<vector<string>>(&labels),
     "add a text label to the plot")
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



  return 0;
}
