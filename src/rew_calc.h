#ifndef rew_calc_h
#define rew_calc_h

#include <string>
#include <vector>
#include <memory>

#include "BHEvent.h"

extern BHEvent event;

namespace LHAPDF {
  class PDF;
  class PDFSet;
}

//-----------------------------------------------
// Function classes to get scales values
//-----------------------------------------------

struct mu_fcn {
  std::string str;
  mu_fcn(const std::string& str);
  virtual double mu() const =0;
};

struct mu_const: public mu_fcn {
  double _mu;
  mu_const(double mu);
  virtual double mu() const;
};

struct mu_fHt: public mu_fcn {
  double fHt;
  mu_fHt(double fHt);
  virtual double mu() const;
};

struct mu_fac_default: public mu_fcn {
  mu_fac_default();
  virtual double mu() const;
};

struct mu_ren_default: public mu_fcn {
  mu_ren_default();
  virtual double mu() const;
};

//-----------------------------------------------
// Reweighting computation
//-----------------------------------------------

class calc_base {
protected:
  static std::vector<std::unique_ptr<const calc_base>> all;
  virtual void calc() const =0;

friend void calc_all_scales(); // run calc for all calcs
};

void calc_all_scales();

class reweighter;

// Factorization --------------------------------

class fac_calc: public calc_base {
protected:
  const mu_fcn* mu_f;
  //const std::vector<LHAPDF::PDF*>& pdfs;
  //bool unc;

  mutable double f[2][5], m[5], lf;

  fac_calc(const mu_fcn* mu_f/*, const std::vector<PDF*>& pdfs*/);
  virtual void calc() const;

public:

friend const fac_calc* mk_fac_calc(const mu_fcn* mu_f/*, const std::vector<PDF*>& pdfs, bool unc=false*/);

friend class reweighter;
};

const fac_calc* mk_fac_calc(const mu_fcn* mu_f/*, const std::vector<PDF*>& pdfs, bool unc=false*/);

// Renormalization ------------------------------

class ren_calc: public calc_base {
protected:
  const mu_fcn* mu_r;
  //const LHAPDF::PDF* pdf;

  mutable double ar, lr;

  ren_calc(const mu_fcn* mu_r/*, const PDF* pdf*/);
  virtual void calc() const;

public:

friend const ren_calc* mk_ren_calc(const mu_fcn* mu_r/*, const PDF* pdf*/);

friend class reweighter;
};

const ren_calc* mk_ren_calc(const mu_fcn* mu_r/*, const PDF* pdf*/);

// Reweighter: combines fac and ren -------------

class reweighter {
  const fac_calc *fac;
  const ren_calc *ren;

  Float_t *weight;

public:
  // Constructor creates branches on tree
  reweighter(const fac_calc* fac, const ren_calc* ren/*, TTree* tree*/);
  ~reweighter();
  void stitch() const;
};

#endif
