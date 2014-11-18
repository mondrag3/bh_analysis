#ifndef rew_calc_h
#define rew_calc_h

#include <string>
#include <vector>
#include <memory>

#include "BHEvent.h"

extern BHEvent event;

// Function to make PDFs
void usePDFset(const std::string& setname);

//-----------------------------------------------
// Function classes to get scales values
//-----------------------------------------------

struct mu_fcn {
  std::string str;
  mu_fcn(const std::string& str);
  virtual double mu() const noexcept =0;
  virtual ~mu_fcn() { }
};

struct mu_const: public mu_fcn {
  double _mu;
  mu_const(double mu);
  virtual double mu() const noexcept;
  virtual ~mu_const() { }
};

struct mu_fHt: public mu_fcn {
  double fHt;
  mu_fHt(double fHt);
  virtual double mu() const noexcept;
  virtual ~mu_fHt() { }
};

struct mu_fac_default: public mu_fcn {
  mu_fac_default();
  virtual double mu() const noexcept;
  virtual ~mu_fac_default() { }
};

struct mu_ren_default: public mu_fcn {
  mu_ren_default();
  virtual double mu() const noexcept;
  virtual ~mu_ren_default() { }
};

//-----------------------------------------------
// Reweighting computation
//-----------------------------------------------

class calc_base {
protected:
  static std::vector<std::unique_ptr<const calc_base>> all;
  virtual void calc() const noexcept =0;

public:
  virtual ~calc_base() { }

friend void calc_all_scales() noexcept; // run calc for all calcs
};

void calc_all_scales() noexcept;

class reweighter;

//-----------------------------------------------
// Factorization --------------------------------
//-----------------------------------------------

class fac_calc: public calc_base {
protected:
  const mu_fcn* mu_f;
  bool unc;

  mutable double f[2][5], m[5], lf;

  fac_calc(const mu_fcn* mu_f, bool unc);
  virtual void calc() const noexcept;

public:
  virtual ~fac_calc();

friend
const fac_calc* mk_fac_calc(const mu_fcn* mu_f, bool unc);

friend class reweighter;
};

const fac_calc* mk_fac_calc(const mu_fcn* mu_f, bool unc=false);

//-----------------------------------------------
// Renormalization ------------------------------
//-----------------------------------------------

class ren_calc: public calc_base {
protected:
  const mu_fcn* mu_r;

  mutable double ar, lr;

  ren_calc(const mu_fcn* mu_r);
  virtual void calc() const noexcept;

public:
  virtual ~ren_calc();

friend
const ren_calc* mk_ren_calc(const mu_fcn* mu_r);

friend class reweighter;
};

const ren_calc* mk_ren_calc(const mu_fcn* mu_r);

//-----------------------------------------------
// Reweighter: combines fac and ren -------------
//-----------------------------------------------

class reweighter {
  const fac_calc *fac;
  const ren_calc *ren;

  mutable double s;

  mutable Float_t weight[3];

public:
  // Constructor creates branches on tree
  reweighter(const fac_calc* fac, const ren_calc* ren, TTree* tree);
  ~reweighter();
  void stitch() const noexcept;
};

#endif
