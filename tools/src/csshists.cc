#include "csshists.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>

#include <boost/regex.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>

#include <TH1.h>

using namespace std;

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

namespace kiwi {

// Properties *******************************************************

enum prop_t { kClass, kBins, kLineColor, kLineWidth };

/*

virtual void	TAttFill::SetFillColor(Color_t fcolor)
virtual void	TAttFill::SetFillColorAlpha(Color_t fcolor, Float_t falpha)
virtual void	TAttFill::SetFillStyle(Style_t fstyle)
virtual void	TH1::SetLabelColor(Color_t color = 1, Option_t* axis = "X")
virtual void	TH1::SetLabelFont(Style_t font = 62, Option_t* axis = "X")
virtual void	TH1::SetLabelOffset(Float_t offset = 0.005, Option_t* axis = "X")
virtual void	TH1::SetLabelSize(Float_t size = 0.02, Option_t* axis = "X")
virtual void	TAttLine::SetLineColor(Color_t lcolor)
virtual void	TAttLine::SetLineColorAlpha(Color_t lcolor, Float_t lalpha)
virtual void	TAttLine::SetLineStyle(Style_t lstyle)
virtual void	TAttLine::SetLineWidth(Width_t lwidth)
virtual void	TAttMarker::SetMarkerColor(Color_t mcolor = 1)
virtual void	TAttMarker::SetMarkerColorAlpha(Color_t mcolor, Float_t malpha)
virtual void	TAttMarker::SetMarkerSize(Size_t msize = 1)
virtual void	TAttMarker::SetMarkerStyle(Style_t mstyle = 1)
virtual void	TH1::SetMaximum(Double_t maximum = -1111)
virtual void	TH1::SetMinimum(Double_t minimum = -1111)
virtual void	TH1::SetName(const char* name)MENU
virtual void	TH1::SetNameTitle(const char* name, const char* title)
virtual void	TH1::SetNdivisions(Int_t n = 510, Option_t* axis = "X")
virtual void	TH1::SetNormFactor(Double_t factor = 1)
virtual void	TH1::SetOption(Option_t* option = " ")
virtual void	TH1::SetStats(Bool_t stats = kTRUE)
virtual void	TH1::SetTickLength(Float_t length = 0.02, Option_t* axis = "X")
virtual void	TH1::SetTitle(const char* title)MENU
virtual void	TH1::SetTitleFont(Style_t font = 62, Option_t* axis = "X")
virtual void	TH1::SetTitleOffset(Float_t offset = 1, Option_t* axis = "X")
virtual void	TH1::SetTitleSize(Float_t size = 0.02, Option_t* axis = "X")
virtual void	TH1::SetXTitle(const char* title)
virtual void	TH1::SetYTitle(const char* title)

*/

struct prop {
  virtual void apply(TH1* h) const =0;
  virtual ~prop() { }
};

struct prop_Class: public prop {
  string name;
  prop_Class(const string& name): name(name) { }
  virtual ~prop_Class() { }
  // this property is used for construction of TH1
  virtual void apply(TH1* h) const { }
};

struct prop_Bins: public prop {
  bool isrange;
  prop_Bins(bool isrange): isrange(isrange) { }
  virtual ~prop_Bins() { }
  // this property is used for construction of TH1
  virtual void apply(TH1* h) const { }
};

struct prop_Bins_range: public prop_Bins {
  Double_t nbinsx, xlow, xup;
  prop_Bins_range(const string& str): prop_Bins(true) {
    stringstream ss(str);
    ss >> nbinsx >> xlow >> xup;
  }
};

struct prop_Bins_vals: public prop_Bins {
  vector<Double_t> xbins;
  prop_Bins_vals(const string& str): prop_Bins(false) {
    stringstream ss(str);
    Double_t bin;
    while(ss >> bin) xbins.push_back(bin);
  }
};

struct prop_int: public prop {
  int x;
  prop_int(const string& str): x( atoi(str.c_str()) ) { }
  virtual ~prop_int() { }
  virtual void apply(TH1* h) const =0;
};

struct prop_LineColor: public prop_int {
  prop_LineColor(const string& str): prop_int(str) { }
  virtual void apply(TH1* h) const {
    h->SetLineColor(x);
  }
};

struct prop_LineWidth: public prop_int {
  prop_LineWidth(const string& str): prop_int(str) { }
  virtual void apply(TH1* h) const {
    h->SetLineWidth(x);
  }
};

pair<prop_t,const prop*> mkprop(const string& str) {
  prop_t type;
  const prop* p;

  size_t sep = str.find_first_of(":");
  pair<string,string> ps(
    boost::algorithm::trim_copy(str.substr(0,sep)),
    boost::algorithm::trim_copy(str.substr(sep+1))
  );

  if (!ps.first.compare("Class")) {
    type = kClass;
    p = new const prop_Class(ps.second);
  } else if (!ps.first.compare("Bins")) {
    type = kBins;
    p = new const prop_Bins_range(ps.second);
  } else if (!ps.first.compare("Bins*")) {
    type = kBins;
    p = new const prop_Bins_vals(ps.second);
  } else if (!ps.first.compare("LineColor")) {
    type = kLineColor;
    p = new const prop_LineColor(ps.second);
  } else if (!ps.first.compare("LineWidth")) {
    type = kLineWidth;
    p = new const prop_LineWidth(ps.second);
  } else throw runtime_error(
    "undefined property "+ps.first
  );

  return make_pair(type,p);
}

typedef boost::unordered_map<prop_t, const prop*> prop_map;

// PIMPL ************************************************************

struct csshists::impl {
  vector< pair<const boost::regex*,prop_map*> > rules;

  ~impl() {
    for (size_t i=0,n=rules.size();i<n;++i) {
      delete rules[i].first;
      for (prop_map::iterator it=rules[i].second->begin(),
           end=rules[i].second->end();it!=end;++it)
      {
        delete it->second;
      }
      delete rules[i].second;
    }
  }

  typedef vector< pair<const boost::regex*,prop_map*> >::iterator iter;
};

// Constructor ******************************************************

csshists::csshists(const string& cssfilename)
: _impl( new impl )
{
  // Read CSS file
  ifstream css(cssfilename.c_str());
  char c;
  bool brak=false, q1=false, q2=false;
  vector< pair< string, vector<string> > > rule_str(1);
  while ( css.get(c) ) {
    if (brak) {
      if (!q2) if (c=='\'') q1 = !q1;
      if (!q1) if (c=='\"') q2 = !q2;

      if (!(q1||q2)) {
        if (c=='}') {
          brak=false;
          // remove whitespace from last property string
          boost::algorithm::trim(rule_str.back().second.back());
          // ref to last rule string
          string& last_rule = rule_str.back().first;
          // remove whitespace from last rule string
          boost::algorithm::trim(last_rule);
          // unquote last rule string
          if ( (last_rule[0]=='\''&&(*(last_rule.end()-1))=='\'') ||
               (last_rule[0]=='\"'&&(*(last_rule.end()-1))=='\"') )
            last_rule = last_rule.substr(1,last_rule.size()-2);
          // start new rule
          rule_str.push_back(pair< string, vector<string> >());
          continue;
        }
        if (c=='{') { throw runtime_error(
          "consecutive { in file \""+cssfilename+"\""
        ); }
      }

      if (c==';') {
        // remove whitespace from last property string
        boost::algorithm::trim(rule_str.back().second.back());
        // start new property string
        rule_str.back().second.push_back(string());
      } else {
        rule_str.back().second.back() += c;
      }

    } else {
      if (c=='{') {
        brak=true;
        // start first property string
        rule_str.back().second.push_back(string());
        continue;
      }
      if (c=='}') { throw runtime_error(
        "} before { in file \""+cssfilename+"\""
      ); }

      rule_str.back().first += c;
    }
  }

  // Convert rules to regex and props
  _impl->rules.reserve(rule_str.size());
  for (size_t i=0,n=rule_str.size();i<n;++i) {
    if (rule_str[i].second.size()==0) continue; // skip blank rule

    _impl->rules.push_back( make_pair(
      new const boost::regex( rule_str[i].first ),
      new prop_map()
    ) );
    //_impl->rules.back().second.reserve( rule_str[i].second.size() );
    for (size_t j=0,m=rule_str[i].second.size();j<m;++j) {
      const string& prop_str = rule_str[i].second[j];
      if (prop_str.size()==0) continue; // skip blank property

      _impl->rules.back().second->insert( mkprop( prop_str ) );
    }
  }
  

  // Test print
  //for (size_t i=0,n=_impl->rules.size();i<n;++i) {
  //  cout << _impl->rules[i].first->str() << endl;
  //  for (prop_map::iterator j=_impl->rules[i].second.begin(),
  //       m=_impl->rules[i].second.end();j!=m;++j)
  //    cout << j->first <<' '<< j->second->str << endl;
  //    //cout <<j<<". " << _impl->rules[i].second[j]->str <<';'<< endl;
  //  cout << "-----------" << endl;
  //}
}

// Make Historgram **************************************************

TH1* csshists::mkhist(const string& name) const {
  prop_map props;
  TH1* h;

  boost::smatch result;
  for (impl::iter it=_impl->rules.begin(),
       end=_impl->rules.end();it!=end;++it) {
    if (boost::regex_match(name, result, *it->first)) {
      for (prop_map::iterator jt=it->second->begin(),
           end2=it->second->end();jt!=end2;++jt)
        props[jt->first] = jt->second;
    }
  }

  if (!props.size()) throw runtime_error(
    "no rules matched for histogram \""+name+"\""
  );

  if (!props.count(kClass)) throw runtime_error(
    "cannot find class for histogram \""+name+"\""
  );
  string class_name(static_cast<const prop_Class*>(props[kClass])->name);

  if (!props.count(kBins)) throw runtime_error(
    "cannot find binning for histogram \""+name+"\""
  );
  if (static_cast<const prop_Bins*>(props[kBins])->isrange) {
    const prop_Bins_range* b
      = static_cast<const prop_Bins_range*>(props[kBins]);

    if (!class_name.compare("TH1F")) {
      h = new TH1F(name.c_str(),"",b->nbinsx,b->xlow,b->xup);
    } else if (!class_name.compare("TH1D")) {
      h = new TH1D(name.c_str(),"",b->nbinsx,b->xlow,b->xup);
    } else if (!class_name.compare("TH1I")) {
      h = new TH1I(name.c_str(),"",b->nbinsx,b->xlow,b->xup);
    } else if (!class_name.compare("TH1C")) {
      h = new TH1C(name.c_str(),"",b->nbinsx,b->xlow,b->xup);
    } else if (!class_name.compare("TH1S")) {
      h = new TH1S(name.c_str(),"",b->nbinsx,b->xlow,b->xup);
    } else throw runtime_error(
      "undefined histogram class "+class_name
    );

  } else {
    const prop_Bins_vals* b
      = static_cast<const prop_Bins_vals*>(props[kBins]);

    if (!class_name.compare("TH1F")) {
      h = new TH1F(name.c_str(),"",b->xbins.size()-1,&b->xbins[0]);
    } else if (!class_name.compare("TH1D")) {
      h = new TH1D(name.c_str(),"",b->xbins.size()-1,&b->xbins[0]);
    } else if (!class_name.compare("TH1I")) {
      h = new TH1I(name.c_str(),"",b->xbins.size()-1,&b->xbins[0]);
    } else if (!class_name.compare("TH1C")) {
      h = new TH1C(name.c_str(),"",b->xbins.size()-1,&b->xbins[0]);
    } else if (!class_name.compare("TH1S")) {
      h = new TH1S(name.c_str(),"",b->xbins.size()-1,&b->xbins[0]);
    } else throw runtime_error(
      "undefined histogram class "+class_name
    );
  }
  props.erase(kClass);
  props.erase(kBins);

  for (prop_map::iterator it=props.begin(),end=props.end();it!=end;++it)
    it->second->apply(h);

  return h;
}

// Destructor *******************************************************

csshists::~csshists() { delete _impl; }

} // end flock namespace
