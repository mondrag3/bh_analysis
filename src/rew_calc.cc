#include "rew_calc.h"

#include <iostream>
#include <sstream>
#include <cmath>

#include <TTree.h>

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"

using namespace std;

// definitions of externs
unique_ptr<const LHAPDF::PDFSet> pdfset;
vector<unique_ptr<const LHAPDF::PDF>> pdfs;

// Function to make PDFs
void usePDFset(const std::string& setname) {
  using namespace LHAPDF;
  pdfset = unique_ptr<const PDFSet>( new PDFSet(setname) );
  vector<PDF*> _pdfs = pdfset->mkPDFs();
  pdfs.reserve( _pdfs.size() );
  for ( auto pdf : _pdfs )
    pdfs.push_back( unique_ptr<const PDF>( pdf ) );
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
double mu_const::mu() const {
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
double mu_fHt::mu() const {
  return fHt*event.Ht();
}

mu_fac_default::mu_fac_default(): mu_fcn("default") { }
double mu_fac_default::mu() const {
  return event.fac_scale;
}

mu_ren_default::mu_ren_default(): mu_fcn("default") { }
double mu_ren_default::mu() const {
  return event.ren_scale;
}

//-----------------------------------------------
// Reweighting computation
//-----------------------------------------------

std::vector<std::unique_ptr<const calc_base>> calc_base::all;

// TODO: Use multiple threads here
void calc_all_scales() {
  for (auto& c : calc_base::all ) c->calc();
}

// Factorization --------------------------------

const fac_calc*
mk_fac_calc(const mu_fcn* mu_f, bool unc) {
  return new fac_calc(mu_f,unc);
}

fac_calc::fac_calc(const mu_fcn* mu_f, bool unc)
: calc_base(), mu_f(mu_f), unc(unc)
{
  all.push_back( unique_ptr<const calc_base>(this) );
}

fac_calc::~fac_calc() {
  delete mu_f;
}

// Renormalization ------------------------------

const ren_calc*
mk_ren_calc(const mu_fcn* mu_r) {
  return new ren_calc(mu_r);
}

ren_calc::ren_calc(const mu_fcn* mu_r)
: calc_base(), mu_r(mu_r)
{
  all.push_back( unique_ptr<const calc_base>(this) );
}

ren_calc::~ren_calc() {
  delete mu_r;
}

// Reweighter: combines fac and ren -------------

// Constructor creates branches on tree
reweighter::reweighter(const fac_calc* fac, const ren_calc* ren, TTree* tree)
: fac(fac), ren(ren), weight( new Float_t )
{
  stringstream ss;
  ss <<  "Fac" << fac->mu_f->str
     << "_Ren" << ren->mu_r->str
     << "_PDF" << pdfset->name();
  if (fac->unc) ss << "_unc";
  else ss << "_cent";
  string name( ss.str() );

  tree->Branch(name.c_str(), weight, (name+"/F").c_str());
}

reweighter::~reweighter() {
  delete weight;
}

//**************************************************************
//
// PHYSICS \/   \/   \/   \/   \/   \/   \/   \/   \/   \/   \/
//
//**************************************************************

//-----------------------------------------------
// Factorization
//-----------------------------------------------

void fac_calc::calc() const {
  cout << "fac_calc: " << mu_f->str << " = " << mu_f->mu() << endl;
}

//-----------------------------------------------
// Renormalization
//-----------------------------------------------

void ren_calc::calc() const {
  cout << "ren_calc: " << mu_r->str << " = " << mu_r->mu() << endl;
}

/*
void ren_calc::calc() const {
  // There is only one global event variable
  // so these references always points to the right place
  static const char& part = event.part[0];
  static const Double_t& alphas  = event.alphas;
  static const Char_t& n = event.alphas_power;
  static const Double_t& ren_scale = event.ren_scale;

  const double mu = mu_r->mu();

  // Calculate α_s change from renormalization
  ar = pow( pdf->alphasQ(mu) / alphas, n );

  // Calculate lr, same as l in Eq (30)
  if (part=='V' || part=='I')
    lr = 2.*log( mu / ren_scale );

}
*/

//-----------------------------------------------
// Stitch
//-----------------------------------------------

void reweighter::stitch() const {
  cout << "Calculating weight for:"
       << " Fac" << fac->mu_f->str
       << " Ren" << ren->mu_r->str << endl;
}
