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
LHAPDF::PDF* pdf; // central PDF

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
  pdf = pdfs[0];
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

/////////////////////////////////////////////////////////////////////
//
// Reweighting adopted from arXiv:1310.7439v1
// "Ntuples for NLO Events at Hadron Colliders" pp. 20 - 25
//
/////////////////////////////////////////////////////////////////////

//-----------------------------------------------
// Factorization
//-----------------------------------------------

void fac_calc::calc() const noexcept {
  // There is only one global event variable
  // so these references always points to the right place
  static const char& part = event.part[0];
  static Double_t* const& usr_wgts = event.usr_wgts;

  const double mu = mu_f->mu();

  // Born & Real
  if (defaultPDF) m[0] = event.weight;
  else {
    for (short i=0;i<2;++i) {
      f[i][0] = pdf->xfxQ(event.id[i], event.x[i], mu)/event.x[i];
    }
    m[0] = event.me_wgt2;
  }

  // Integrated subtraction
  if (part=='I') {

    // Calculate terms in Eq. (43)
    lf = 2.*log( mu / event.fac_scale );
    for (short i=1;i<9;++i)
      m[i] = usr_wgts[i+1] + usr_wgts[i+9]*lf;

    // Calculate terms in Eq. (44)
    static Double_t x, xp;
    static Int_t id;
    double si_[2] = { 0., 0. };
    for (short i=0;i<2;++i) {
      id = event.id[i];
      x  = event.x  [i];
      xp = event.xp [i];

      f[i][1] = ( id==21 // Eq. (46)
              ? quark_sum(    x, mu)/x
              : pdf->xfxQ(id, x, mu)/x
      );
      f[i][2] = ( id==21 // Eq. (47)
              ? quark_sum(    x/xp, mu)/x
              : pdf->xfxQ(id, x/xp, mu)/x
      );
      f[i][3] = pdf->xfxQ(21, x,    mu)/x; // Eq. (48)
      f[i][4] = pdf->xfxQ(21, x/xp, mu)/x; // Eq. (49)
    }
    for (short i=0;i<2;++i)
      for (short j=1;j<5;++j)
        si_[i] += f[i][j]*m[j+(i ? 4 : 0)]; // si_[0] = f_1^(j)*ω_j
                                            // si_[1] = f_2^(j)*ω_j+4

    si = f[1][0]*si_[0] + f[0][0]*si_[1]; // f_1^(j)*ω_j*f_2 + f_1*f_2^(j)*ω_j+4

  }
}

Double_t fac_calc::quark_sum(Double_t x, Double_t mu_fac) const {
  static const int nquarks = 10;
  static Int_t q[nquarks] = {
    1, 2, 3, 4, 5,
   -1,-2,-3,-4,-5,
  };

  Double_t f = 0.;
  for (int i=0;i<nquarks;++i) f += pdf->xfxQ(q[i], x, mu_fac);

  return f;
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
  if (!defaultPDF) ar = pow( pdf->alphasQ(mu) / alphas, n );

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

  if (part=='I') s += fac->si;

  weight[0] = s * ren->ar;

  if (!isnormal(weight[0])) {
    cerr << "\033[31mEvent " << event.eid << "\033[0m: "
         << "weight=" << weight[0] << endl;
    weight[0] = 0.;
  } /*else {
    cout << "\033[32mEvent " << event.eid << "\033[0m: "
         << "weight=" << weight[0];
  }
  cout << " fac=" << static_cast<const mu_fHt*>(fac->mu_f)->fHt
       << " ren=" << static_cast<const mu_fHt*>(ren->mu_r)->fHt << endl;*/
}
