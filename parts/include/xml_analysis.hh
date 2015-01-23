#ifndef xml_analysis_hh
#define xml_analysis_hh

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <forward_list>
#include <unordered_map>
#include <stdexcept>
#include <memory>

#include <TDirectory.h>
#include <TLorentzVector.h>
#include <TH1.h>

#include <boost/optional.hpp>

namespace rapidxml {
  template<class T> class xml_node;
}

class weight;
class csshist;

// Particle struct **************************************************

struct kinematics {
  bool use_mass, use_pt, use_eta, use_phi, use_tau;
  double mass, pt, eta, phi, tau;
  kinematics();
};

struct particle {
  std::string name;
  boost::optional<int> pid;
  boost::optional<double> pt_cut, eta_cut;
  std::forward_list<particle> daughters;

  mutable TLorentzVector p;
  mutable kinematics vars;

  static std::unordered_map<std::string,particle*> by_name;

  particle(const rapidxml::xml_node<> *node);
  void operator()
};

std::ostream& operator<<(std::ostream& out, const particle& p) {
  out << p.name;
  if (p.pid) out << " pid=" << *p.pid;
  if (p.pt_cut) out << " pt>" << *p.pt_cut;
  if (p.eta_cut) out << " |eta|<" << *p.eta_cut;
  out << std::endl;
  for (auto& d : p.daughters) cout << "  " << d;
  return out;
}

// Histogram wrapper ************************************************

class hist {
public:
  hist(const double *var);
  virtual void Fill(Double_t x) noexcept =0;
};

class hist {
  const double *var;
  std::unordered_map<const weight*,TH1*> h;

public:
  hist(const std::string& name, const double *var);
  virtual void Fill(Double_t x) noexcept;

  static std::unique_ptr<const csshists> css;
  static std::unordered_map<const weight*,TDirectory*> dirs;
};

// Analysis definition **********************************************

struct analysis_def {
  const BHEvent *event;
  const std::vector<TLorentzVector> *jets_;

  std::forward_list<particle> particles;

  size_t njets;
  boost::optional<double> jet_pt_cut, jet_eta_cut;
  std::vector<kinematics> jets;

  std::forward_list<hist> hists;

  analysis_def(const string& xml_file,
               const BHEvent *event, const std::vector<TLorentzVector> *jets);
};

#endif
