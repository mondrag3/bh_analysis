#include "rew_calc.h"

#include <iostream>
#include <sstream>
#include <cmath>

#include <TTree.h>

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;

template<typename T> inline T sq(T x) { return x*x; }

// PDF variables
string pdfname("def");
const LHAPDF::PDFSet* pdfset(nullptr);
vector<LHAPDF::PDF*> pdfs;

// PDF garbage collector
struct PDFgc {
  inline void clear() {
    delete pdfset;
    for ( auto pdf : pdfs ) delete pdf;
  }
  ~PDFgc() { clear(); }
} __pdf;

// Function to make PDFs
void usePDFset(const std::string& setname) {
  __pdf.clear();
  pdfset = new LHAPDF::PDFSet(setname);
  pdfname = pdfset->name();
  pdfs = pdfset->mkPDFs();
}

//-----------------------------------------------
// Function classes to get scales values
//-----------------------------------------------

mu_fcn::mu_fcn(const string& str): str(str) { }

mu_const::mu_const(double mu)
: mu_fcn( [&mu] () { // lambda
    stringstream ss;
    ss << mu << "GeV";
    return ss.str();
  } () ),
  _mu(mu)
{ }
double mu_const::mu() const noexcept {
  return _mu;
}

mu_fHt::mu_fHt(double fHt)
: mu_fcn( [&fHt] () { // lambda
    stringstream ss;
    ss << fHt << "Ht";
    return ss.str();
  } () ),
  fHt(fHt)
{ }
double mu_fHt::mu() const noexcept {
  return fHt*event.Ht();
}

mu_fac_default::mu_fac_default(): mu_fcn("def") { }
double mu_fac_default::mu() const noexcept {
  return event.fac_scale;
}

mu_ren_default::mu_ren_default(): mu_fcn("def") { }
double mu_ren_default::mu() const noexcept {
  return event.ren_scale;
}

//-----------------------------------------------
// Reweighting computation
//-----------------------------------------------

std::vector<std::unique_ptr<const calc_base>> calc_base::all;

// TODO: Use multiple threads here
void calc_all_scales() noexcept {
  for (auto& c : calc_base::all ) c->calc();
}

// Factorization --------------------------------

const fac_calc*
mk_fac_calc(const mu_fcn* mu_f, bool unc, bool defaultPDF) {
  return new fac_calc(mu_f,unc,defaultPDF);
}

fac_calc::fac_calc(const mu_fcn* mu_f, bool unc, bool defaultPDF)
: calc_base(), mu_f(mu_f), unc(unc), defaultPDF(defaultPDF)
{
  all.push_back( unique_ptr<const calc_base>(this) );
}

fac_calc::~fac_calc() {
  delete mu_f;
}

// TODO: Add a check for branches with the same name

// Renormalization ------------------------------

const ren_calc*
mk_ren_calc(const mu_fcn* mu_r, bool defaultPDF) {
  return new ren_calc(mu_r,defaultPDF);
}

ren_calc::ren_calc(const mu_fcn* mu_r, bool defaultPDF)
: calc_base(), mu_r(mu_r), defaultPDF(defaultPDF)
{
  if (defaultPDF) ar = 1;
  all.push_back( unique_ptr<const calc_base>(this) );
}

ren_calc::~ren_calc() {
  delete mu_r;
}

// Reweighter: combines fac and ren -------------

// Constructor creates branches on tree
reweighter::reweighter(const fac_calc* fac, const ren_calc* ren, TTree* tree)
: fac(fac), ren(ren)
{
  if (!pdfset) {
    cerr << "\033[31mNo PDF loaded\033[0m"  << endl;
    exit(1);
  }

  stringstream ss;
  ss <<  "Fac" << fac->mu_f->str
     << "_Ren" << ren->mu_r->str
     << "_PDF" << pdfname;
  if (fac->unc) ss << "_unc";
  else ss << "_cent";
  string name( ss.str() );

  cout << "Creating branch: " << name << endl;

  tree->Branch(name.c_str(), &weight[0], (name+"/F").c_str());
}

reweighter::~reweighter() { }

//**************************************************************
//
// PHYSICS \/   \/   \/   \/   \/   \/   \/   \/   \/   \/   \/
//
//**************************************************************

//-----------------------------------------------
// Factorization
//-----------------------------------------------

void fac_calc::calc() const noexcept {

  const double mu = mu_f->mu();

  // Born & Real
  if (defaultPDF) m[0] = event.weight;
  else {
    for (short i=0;i<2;++i) {
      f[i][0] = pdfs[0]->xfxQ(event.iid[i], event.x[i], mu)/event.x[i];
    }
    m[0] = event.me_wgt2;
  }

  // TODO: Integrated subtraction


}

//-----------------------------------------------
// Renormalization
//-----------------------------------------------

void ren_calc::calc() const noexcept {
  // There is only one global event variable
  // so these references always points to the right place
  static const char& part = event.part[0];
  static const Double_t& alphas  = event.alphas;
  static const Char_t& n = event.alphas_power;
  static const Double_t& ren_scale = event.ren_scale;

  const double mu = mu_r->mu();

  // Calculate α_s change from renormalization
  if (!defaultPDF) ar = pow( pdfs[0]->alphasQ(mu) / alphas, n );

  // Calculate lr, same as l in Eq (30)
  if (part=='V' || part=='I')
    lr = 2.*log( mu / ren_scale );

}

//-----------------------------------------------
// Stitch
//-----------------------------------------------

void reweighter::stitch() const noexcept {
  // There is only one global event variable
  // so these references always points to the right place
  static const char& part = event.part[0];
  static Double_t* const& usr_wgts = event.usr_wgts;

  s = fac->m[0]; // m0

  if (part=='V' || part=='I')
    s += (ren->lr)*usr_wgts[0] + 0.5*sq(ren->lr)*usr_wgts[1]; // m0

  for (short i=0;i<2;++i) s *= fac->f[i][0]; // m0 becomes s

  // TODO: if part=='I' s += sum over quarks


  weight[0] = s * ren->ar;

  if (!isnormal(weight[0])) {
    cerr << "\033[31mEvent " << event.id << "\033[0m: "
         << "weight=" << weight[0] << endl;
    weight[0] = 0.;
  }
}
