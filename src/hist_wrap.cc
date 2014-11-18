#include "hist_wrap.h"

#include <fstream>

#include <boost/regex.hpp>

#include <TH1.h>

using namespace std;

// binning ************************************************

hist::binning::binning(): nbins(100), min(0.), max(100.) { }
hist::binning::~binning() { }

// hist ***************************************************

hist::hist(): b(), underflow(0,0.), overflow(0,0.), h(NULL) { }
hist::~hist() { } // doesn't delete TH1F* because TFile will do that

hist::hist(const string& name,const string& title)
: b(get_binning(name)), underflow(0,b.min), overflow(0,b.max),
  h( new TH1F(name.c_str(),title.c_str(),b.nbins,b.min,b.max) )
{
  all.push_back(this);
}
void hist::Fill(double x, double w) {
  h->Fill(x,w);
  if (x>=b.max) {
    ++overflow.first;
    if (x>overflow.second) overflow.second = x;
  }
  else if (x<b.min) {
    ++underflow.first;
    if (x<underflow.second) underflow.second = x;
  }
}

void hist::FillOverflow(double weight) {
  Fill(h->GetNbinsX()+1,weight);
}

void hist::read_binnings(const char* filename, const char* regex) {
  ifstream binsfile(filename);

  string hname;
  binning b;
  while ( binsfile >> hname ) {
    binsfile >> b.nbins;
    binsfile >> b.min;
    binsfile >> b.max;

    if (hname.compare("default")) binnings[hname] = b;
    else default_binning = b;
  }

  if (regex) binning_name_regex_pattern = regex;
}

ostream& hist::print_overflow(ostream& out) {
  for (vector<const hist*>::iterator it=all.begin(), end=all.end(); it<end; ++it) {
    const hist* h = *it;
    if (h->underflow.first) {
      out << "Underflow in " << h->h->GetName()
          << ": N="<<h->underflow.first << " min="<<h->underflow.second << endl;
    }
    if (h->overflow.first) {
      out << "Overflow in " << h->h->GetName()
          << ": N="<<h->overflow.first << " max="<<h->overflow.second << endl;
    }
  }
  return out;
}

void hist::delete_all() {
  for (vector<const hist*>::iterator it=all.begin(), end=all.end(); it<end; ++it)
    delete *it;
}

const hist::binning hist::get_binning(const string& hist_name) {
  if (binning_name_regex_pattern.size()) {

    static const boost::regex patern(binning_name_regex_pattern);

    boost::smatch result;
    if (boost::regex_search(hist_name, result, patern))
      return binnings[string(result[0].first, result[0].second)];
    else return default_binning;

  } else {
    if (binnings.count(hist_name)) return binnings[hist_name];
    else return default_binning;
  }
}

map<string,hist::binning> hist::binnings;
vector<const hist*> hist::all;
string hist::binning_name_regex_pattern;
hist::binning hist::default_binning;
