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
// Factorization --------------------------------
//-----------------------------------------------

class reweighter;

class fac_calc {
protected:
  const mu_fcn* mu_f;
  bool pdf_unc;
  bool defaultPDF;

  mutable double f[3][2][5], m[9], lf, si[3];

  fac_calc(const mu_fcn* mu_f, bool pdf_unc, bool defaultPDF) noexcept;
  void calc() const noexcept;

public:
  ~fac_calc();

friend class reweighter;
};

//-----------------------------------------------
// Renormalization ------------------------------
//-----------------------------------------------

enum class alphas_fcn: char { all_mu, two_mH };

class ren_calc {
protected:
  const mu_fcn* mu_r;
  alphas_fcn asfcn;
  bool defaultPDF;

  mutable double ar, lr, m0;

  ren_calc(const mu_fcn* mu_r, alphas_fcn asfcn, bool defaultPDF) noexcept;
  void calc() const noexcept;

public:
  ~ren_calc();

friend class reweighter;
};

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
