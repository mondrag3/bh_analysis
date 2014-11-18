#include "rew_calc.h"

#include <iostream>
#include <sstream>
#include <cmath>

using namespace std;

//-----------------------------------------------
// Function classes to get scales values
//-----------------------------------------------

mu_fcn::mu_fcn(const string& str): str(str) { }

mu_const::mu_const(double mu)
: mu_fcn( [&mu] () {
    stringstream ss;
    ss << "const" << mu;
    return ss.str();
  } () ),
  _mu(mu)
{ }
double mu_const::mu() const {
  return _mu;
}

mu_fHt::mu_fHt(double fHt)
: mu_fcn( [&fHt] () {
    stringstream ss;
    ss << "Ht" << fHt;
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
mk_fac_calc(const mu_fcn* mu_f/*, const std::vector<PDF*>& pdfs, bool unc=false*/) {
  return new fac_calc(mu_f);
}

fac_calc::fac_calc(const mu_fcn* mu_f/*, const vector<PDF*>& pdfs*/)
: calc_base(), mu_f(mu_f) //, pdfs(pdfs), unc(unc)
{
  all.push_back( unique_ptr<const calc_base>(this) );
}

void fac_calc::calc() const {
  cout << "fac_calc: " << mu_f->str << endl;
}

// Renormalization ------------------------------

const ren_calc*
mk_ren_calc(const mu_fcn* mu_r/*, const PDF* pdf*/) {
  return new ren_calc(mu_r);
}

ren_calc::ren_calc(const mu_fcn* mu_r/*, const PDF* pdf*/)
: calc_base(), mu_r(mu_r) //, pdf(pdf)
{
  all.push_back( unique_ptr<const calc_base>(this) );
}

void ren_calc::calc() const {
  cout << "ren_calc: " << mu_r->str << endl;
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

// Reweighter: combines fac and ren -------------

// Constructor creates branches on tree
reweighter::reweighter(const fac_calc* fac, const ren_calc* ren/*, TTree* tree*/)
: fac(fac), ren(ren), weight( new Float_t )
{
  
}

reweighter::~reweighter() {
  delete weight;
}

void reweighter::stitch() const {
  cout << "Calculating weight for:"
       << " Fac" << fac->mu_f->str
       << " Ren" << ren->mu_r->str << endl;
}

