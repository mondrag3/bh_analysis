#include "xml_analysis.hh"

#include "rapidxml-1.13/rapidxml.hpp"

#include "weight.hh"
#include "csshists.hh"

using namespace std;
using namespace rapidxml;
using boost::optional;

inline xml_attribute<char>*
safe_attribute(const xml_node<> *node, const char* name) {
  xml_attribute<char> *attr = node->first_attribute(name);
  if (attr) return attr;
  else {
    stringstream ss;
    ss << "XML: Node " << node->name() << " has no attribute " << name;
    throw runtime_error(ss.str());
  }
}

inline int read_int(const xml_attribute<char> *attr) {
  if (attr) return atoi(attr->value());
  else return 0;
}
inline optional<int> read_opt_int(const xml_attribute<char> *attr) {
  if (attr) return atoi(attr->value());
  else return boost::none;
}
inline optional<double> read_opt_double(const xml_attribute<char> *attr) {
  if (attr) return atof(attr->value());
  else return boost::none;
}

// Particle struct **************************************************

class particle;

void add_particles(const xml_node<> *node, forward_list<particle>& pl) {
  auto it = pl.before_begin();
  for (xml_node<> *p = node->first_node("particle"); p;
       p = p->next_sibling("particle")) it = pl.emplace_after(it,p);
}

kinematics::kinematics()
: use_mass(false), use_pt(false), use_eta(false), use_phi(false), use_tau(false)
{ }

struct particle {

  string name;
  optional<int> pid;
  optional<double> pt_cut, eta_cut;
  forward_list<particle> daughters;
  mutable kinematics vars;

  static unordered_map<string,particle*> by_name;

  particle(const xml_node<> *node)
  : name( node->first_attribute("name")->value() ),
    pid ( read_opt_int(node->first_attribute("pid")) ),
    pt_cut  ( read_opt_double(node->first_attribute("pt_cut")) ),
    eta_cut ( read_opt_double(node->first_attribute("eta_cut")) )
  {
    add_particles(node,daughters);

    if (by_name.count(name)) {
      cerr << "Particle with name " << name << "already exists" << endl;
    } else {
      by_name[name] = this;
    }
  }

};
unordered_map<string,particle*> particle::by_name;

ostream& operator<<(ostream& out, const particle& p) {
  out << p.name;
  if (p.pid) out << ' ' << *p.pid;
  if (p.pt_cut) out << " pt>" << *p.pt_cut;
  if (p.eta_cut) out << " |eta|<" << *p.eta_cut;
  out << endl;
  for (auto& d : p.daughters) cout << "  " << d;
  return out;
}

// Histogram wrappers ***********************************************
hist::hist(const string& name, const double *var): var(var) {
  TH1* hist = css->mkhist(name);
  // hist->Sumw2(false); // in ROOT6 true seems to be the default
  for (auto& wt : weight::all) {
    const weight *w = wt.get();
    dirs[w]->cd();
    h[w] = static_cast<TH1*>( hist->Clone() );
  }
  delete hist;
}

virtual void hist::Fill(Double_t x) noexcept {
  for (auto& wt : weight::all)
    h[wt.get()]->Fill(x,wt->is_float ? wt->w.f : wt->w.d);
}

unique_ptr<const csshists> hist::css;
unordered_map<const weight*,TDirectory*> hist::dirs;

// Analysis definition **********************************************
struct analysis_def {
  optional<double> jet_pt_cut, jet_eta_cut;
  forward_list<particle> particles;
  vector<kinematics> jets;
  vector<hist> hists;

  analysis_def(const string& xml_file,
               const BHEvent *event, const vector<TLorentzVector> *jets)
  : event(event), jets_(jets)
  {
    xml_document<> doc;
  	// Read the xml file into a vector
  	ifstream file(xml_file);
  	vector<char> buffer((istreambuf_iterator<char>(file)),
                        istreambuf_iterator<char>());
  	buffer.push_back('\0');
    file.close();
  	// Parse the buffer using the xml file parsing library into doc
  	doc.parse<0>(buffer.data());
    // Find root node
  	const xml_node<> *analysis_node = doc.first_node("analysis");

    const xml_node<> *particles_node = analysis_node->first_node("particles");
    if (particles_node) add_particles(particles_node,particles);

    const xml_node<> *jets_node = analysis_node->first_node("jets");
    if (!jets_node) {
      cout << "Warning: No jets node in xml file " << xml_file << endl;
    }
    njets = read_int(safe_attribute(jets_node,"num"));
    jets.resize(njets);
    jet_pt_cut  = read_opt_double(jets_node->first_attribute("pt_cut"));
    jet_eta_cut = read_opt_double(jets_node->first_attribute("eta_cut"));

    const xml_node<> *histograms_node = analysis_node->first_node("histograms");
    for (xml_node<> *h = histograms_node->first_node("hist"); h;
         h = h->next_sibling("hist"))
    {
      string var( safe_attribute(h,"var")->value() );
      string pstr( safe_attribute(h,"particles")->value() );
      forward_list<string> p_; // particles
      auto it = p_.before_begin();
      size_t comma;
      while ((comma=pstr.find(','))!=string::npos) {
        it = p_.emplace_after(it,pstr,0,comma);
        pstr.erase(0,comma+1);
      }
      cout << var << ": ";
      for (auto& p : p_) cout << p << ' ';
      cout << endl;
    }
  }
};
