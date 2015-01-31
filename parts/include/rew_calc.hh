#ifndef rew_calc_h
#define rew_calc_h

#include <string>
#include <vector>
#include <memory>

#include "BHEvent.hh"

extern BHEvent event;

// Function to make PDFs
void usePDFset(const std::string& setname);

//-----------------------------------------------
// Function classes to get scales values
//-----------------------------------------------

struct mu_fcn {
  std::string str;
  mu_fcn(const std::string& str) noexcept;
  virtual double mu() const noexcept =0;
  virtual ~mu_fcn() { }
};

struct mu_fixed: public mu_fcn {
  double _mu;
  mu_fixed(double mu) noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fixed() { }
};

struct mu_fHt: public mu_fcn {
  double fHt;
  mu_fHt(double fHt) noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fHt() { }
};

struct mu_fHt_Higgs: public mu_fcn {
  double fHt;
  mu_fHt_Higgs(double fHt) noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fHt_Higgs() { }
};

struct mu_fac_default: public mu_fcn {
  mu_fac_default() noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fac_default() { }
};

struct mu_ren_default: public mu_fcn {
  mu_ren_default() noexcept;
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
  bool pdf_unc;
  bool defaultPDF;

  mutable double f[3][2][5], m[9], lf, si[3];

  fac_calc(const mu_fcn* mu_f, bool pdf_unc, bool defaultPDF) noexcept;
  virtual void calc() const noexcept;

  // Double_t quark_sum(Double_t x, Double_t mu_fac) const;

public:
  virtual ~fac_calc();

friend const fac_calc* mk_fac_calc(const mu_fcn*, bool, bool) noexcept;

friend class reweighter;
};

const fac_calc* mk_fac_calc(
  const mu_fcn* mu_f, bool pdf_unc=false, bool defaultPDF=false) noexcept;

//-----------------------------------------------
// Renormalization ------------------------------
//-----------------------------------------------

enum class alphas_fcn: char { all_mu, two_mH };

class ren_calc: public calc_base {
protected:
  const mu_fcn* mu_r;
  alphas_fcn asfcn;
  bool defaultPDF;

  mutable double ar, lr, m0;

  ren_calc(const mu_fcn* mu_r, alphas_fcn asfcn, bool defaultPDF) noexcept;
  virtual void calc() const noexcept;

public:
  virtual ~ren_calc();

friend const ren_calc* mk_ren_calc(const mu_fcn*, alphas_fcn, bool) noexcept;

friend class reweighter;
};

const ren_calc* mk_ren_calc(
  const mu_fcn* mu_r, alphas_fcn asfcn=alphas_fcn::all_mu, bool defaultPDF=false
) noexcept;

//-----------------------------------------------
// Reweighter: combines fac and ren -------------
//-----------------------------------------------

class reweighter {
  const fac_calc *fac;
  const ren_calc *ren;
  bool pdf_unc;
  short nk;

  mutable double s[3];

  mutable Float_t weight[3];

public:
  // Constructor creates branches on tree
  reweighter(const fac_calc* fac, const ren_calc* ren,
             TTree* tree, bool pdf_unc=false);
  ~reweighter();
  void stitch() const noexcept;
};

#endif
