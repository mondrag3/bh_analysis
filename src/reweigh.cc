#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <stdexcept>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>

#include "rapidxml-1.13/rapidxml.hpp"

#include "rew_calc.hh"
#include "timed_counter.hh"

using namespace std;
namespace po = boost::program_options;

// using xml_doc  = rapidxml::xml_document<>;
using xml_node = rapidxml::xml_node<>;
using xml_attr = rapidxml::xml_attribute<char>;

inline const char* get_attr(const xml_node *node, const char* name) {
  xml_attr *attr = node->first_attribute(name);
  if (attr) return attr->value();
  else throw runtime_error(
    string("XML: Node ") + node->name() + " has no attribute " + name
  );
}

#define node_loop_all(root) \
  for (xml_node *node = root->first_node(); node; \
       node = node->next_sibling())

#define node_loop(root,name) \
  for (xml_node *node = root->first_node(name); node; \
       node = node->next_sibling(name))

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

BHEvent event; // extern

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  string BH_file, weights_file, pdf_set, xml_file;
  bool old_bh, counter_newline;
  Long64_t num_events = 0;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
      ("help,h", "produce help message")
      ("bh", po::value<string>(&BH_file)->required(),
       "*input event root file (Blackhat ntuple)")
      ("config,c", po::value<string>(&xml_file)->required(),
       "*configuration XML file")
      ("output,o", po::value<string>(&weights_file)->required(),
       "*output root file with new event weights")
      ("pdf", po::value<string>(&pdf_set)->default_value("CT10nlo"),
       "LHAPDF set name")
      ("num-events,n", po::value<Long64_t>(&num_events),
       "process only this many events. Zero means all.")
      ("old-bh", po::bool_switch(&old_bh),
       "read an old BH tree (no part & alphas_power branches)")
      ("counter-newline", po::bool_switch(&counter_newline),
       "do not overwrite previous counter message")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, all_opt), vm);
    if (argc == 1 || vm.count("help")) {
      cout << all_opt << endl;
      exit(0);
    }
    po::notify(vm);
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Open input event file
  TFile *fin = new TFile(BH_file.c_str(),"READ");
  if (fin->IsZombie()) exit(1);

  cout << "Input BH event file: " << fin->GetName() << endl;

  TTree *tin = (TTree*)fin->Get("t3");

  // Find number of events to process
  if (num_events>0) num_events = min( num_events, tin->GetEntries() );
  else num_events = tin->GetEntries();

  // Set up BlackHat event
  event.SetTree(tin, BHEvent::reweighting, old_bh);

  tin->SetBranchAddress("weight", &event.weight);

  if (old_bh) {

    if      (BH_file.find("born")!=string::npos) event.SetPart('B');
    else if (BH_file.find("real")!=string::npos) event.SetPart('R');
    else if (BH_file.find("loop")!=string::npos) event.SetPart('V');
    else if (BH_file.find("vsub")!=string::npos) event.SetPart('I');
    else {
      cerr << "\033[31mCannot determine part type from file name\033[0m" << endl;
      exit(1);
    }

    const size_t jpos = BH_file.find_first_of('j');
    if (jpos!=string::npos) {
      const char njc = BH_file[jpos-1];
      if (isdigit(njc)) {

        Char_t nj = njc - '0';
        if (event.part[0]=='V' || event.part[0]=='I') ++nj;
        event.SetAlphasPower(nj);

      } else {
        cerr << "\033[31mCannot determine number of jets from file name\033[0m" << endl;
        exit(1);
      }
    } else {
      cerr << "\033[31mCannot determine number of jets from file name\033[0m" << endl;
      exit(1);
    }
  }

  // Load PDF
  cout << endl;
  usePDFset(pdf_set);
  cout << endl;

  // Open output weights file
  TFile *fout = new TFile(weights_file.c_str(),"recreate");
  if (fout->IsZombie()) exit(1);

  cout << "Output weights file: " << fout->GetName() << endl;

  TTree *tree = new TTree("weights","");

  // Setup new weights **********************************************
  unordered_map<string,const mu_fcn*> mu;
  unordered_map<string,const fac_calc*> fac;
  unordered_map<string,const ren_calc*> ren;
  vector<const reweighter*> weights;

  rapidxml::xml_document<> doc;
	// Read the xml file into a vector
	ifstream file(xml_file);
	vector<char> buffer((istreambuf_iterator<char>(file)),
                      istreambuf_iterator<char>());
	buffer.push_back('\0');
  file.close();
	// Parse the buffer using the xml file parsing library into doc
	doc.parse<0>(buffer.data());
  // Find root node
	const xml_node *format_node   = doc.first_node("bh_format");
	const xml_node *energies_node = doc.first_node("energies");
	const xml_node *scales_node   = doc.first_node("scales");
	const xml_node *weights_node  = doc.first_node("weights");

  if (const xml_attr* alphas = format_node->first_attribute("alphas")) {
    if (!strcmp(alphas->value(),"two_mH"))
      ren_calc::bh_alphas = alphas_fcn::two_mH;
  }

  node_loop_all(energies_node) {
    const char* tag_name = node->name();
    const char* name = get_attr(node,"name");
    if (mu.count(name)) {
      cerr << "Warning: already existing energy definition " << name
           << " is replaced" << endl;
      delete mu[name];
    }
    if (!strcmp(tag_name,"Ht"))
      mu[name] = new mu_fHt(atof(get_attr(node,"frac")));
    else if (!strcmp(tag_name,"Ht_Higgs"))
      mu[name] = new mu_fHt_Higgs(atof(get_attr(node,"frac")));
    else if (!strcmp(tag_name,"fixed"))
      mu[name] = new mu_fixed(atof(get_attr(node,"val")));
    else if (!strcmp(tag_name,"fac_default"))
      mu[name] = new mu_fac_default();
    else if (!strcmp(tag_name,"ren_default"))
      mu[name] = new mu_ren_default();
    else
      cerr << "Warning: unrecognized energy definition: " << tag_name << endl;
  }

  node_loop(scales_node,"fac") {
    const char* name = get_attr(node,"name");
    if (fac.count(name)) {
      cerr << "Warning: already existing scale definition " << name
           << " is replaced" << endl;
      delete fac[name];
    }
    const xml_attr* pdfunc = node->first_attribute("pdfunc");
    fac[name] = new fac_calc(
      mu[get_attr(node,"energy")], pdfunc && !strcmp(pdfunc->value(),"true")
    );
  }

  node_loop(scales_node,"ren") {
    const char* name = get_attr(node,"name");
    if (ren.count(name)) {
      cerr << "Warning: already existing scale definition " << name
           << " is replaced" << endl;
      delete ren[name];
    }
    ren[name] = new ren_calc( mu[get_attr(node,"energy")] );
  }

  node_loop(weights_node,"weight") {
    const xml_attr* pdfunc = node->first_attribute("pdfunc");
    const char* fac_name = get_attr(node,"fac");
    const char* ren_name = get_attr(node,"ren");
    weights.push_back( new reweighter(
      make_pair(fac[fac_name],fac_name),
      make_pair(ren[ren_name],ren_name),
      tree,
      pdfunc && !strcmp(pdfunc->value(),"true")
    ) );
  }

  // Reading events from the input ntuple ***************************
  cout << "\nReading " << num_events << " events" << endl;
  timed_counter counter(counter_newline);

  cout << scientific;
  cout.precision(10);

  for (Long64_t ent=0; ent<num_events; ++ent) {
    counter(ent);
    tin->GetEntry(ent);

    // use event id for event number
    event.eid = ent;

    // REWEIGHTING
    for (auto f : fac) f.second->calc();
    for (auto r : ren) r.second->calc();
    for (auto w : weights) w->stitch();
    tree->Fill();
  }
  counter.prt(num_events);
  cout << endl;

  fout->Write();
  cout << "\n\033[32mWrote\033[0m: " << fout->GetName() << endl;
  fout->Close();
  delete fout;

  fin->Close();
  delete fin;

  for (auto m : mu)  delete m.second;
  for (auto f : fac) delete f.second;
  for (auto r : ren) delete r.second;
  for (auto w : weights) delete w;

  return 0;
}
