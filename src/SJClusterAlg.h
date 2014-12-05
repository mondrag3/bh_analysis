#ifndef SJClusterAlg_h
#define SJClusterAlg_h

#include <string>
#include <vector>
#include <memory>

#include "TLorentzVector.h"

class TTree;

struct SJClusterAlg {
  Int_t N;
  std::vector<Float_t> *eta, *phi, *e, *mass, *pt, *numC;
  std::vector<Int_t> *ind;
  // these have to be vector pointers
  // that's how data is stored
  // they don't have to be deleted (no memory leak occurs)
  // For some ROOT hocus-pocus reason, these pointers have
  // to be initialized to zero

  const std::string name;

  std::vector<TLorentzVector> jetsByPt(double pt_cut, double eta_cut) const;

  SJClusterAlg(TTree* tree, const std::string& name);
  ~SJClusterAlg();

  static std::vector<std::unique_ptr<const SJClusterAlg>> all;
  static void add(TTree* tree, const std::string& name);
};

#endif
