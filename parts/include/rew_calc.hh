#ifndef rew_calc_h
#define rew_calc_h

#include <string>
#include <vector>
#include <algorithm>

#include "BHEvent.hh"

extern BHEvent event;

// Function to make PDFs
void usePDFset(const std::string& setname);

//-----------------------------------------------
// Function classes to get scales values
//-----------------------------------------------

class mu_fcn {
public:
  virtual double mu() const noexcept =0;
  virtual ~mu_fcn() { }
};

class mu_fixed: public mu_fcn {
  double _mu;
public:
  mu_fixed(double mu) noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fixed() { }
};

class mu_fHt: public mu_fcn {
  double fHt;
public:
  mu_fHt(double fHt) noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fHt() { }
};

class mu_fHt_Higgs: public mu_fcn {
  double fHt;
public:
  mu_fHt_Higgs(double fHt) noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fHt_Higgs() { }
};

class mu_fac_default: public mu_fcn {
public:
  mu_fac_default() noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fac_default() { }
};

class mu_ren_default: public mu_fcn {
public:
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

public:
  fac_calc(const mu_fcn* mu_f, bool pdf_unc=false, bool defaultPDF=false) noexcept;
  void calc() const noexcept;
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
  bool defaultPDF;

  mutable double ar, lr, m0;

public:
  ren_calc(const mu_fcn* mu_r, bool defaultPDF=false) noexcept;
  void calc() const noexcept;
  ~ren_calc();

  static alphas_fcn bh_alphas;

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
  reweighter(const std::pair<const fac_calc*,std::string>& fac,
             const std::pair<const ren_calc*,std::string>& ren,
             TTree* tree, bool pdf_unc=false);
  ~reweighter();
  void stitch() const noexcept;
};

#endif
