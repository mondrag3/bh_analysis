#ifndef SJClusterAlg_h
#define SJClusterAlg_h

#include <string>
#include <vector>

#include "TLorentzVector.h"

class TTree;

class SJClusterAlg;
typedef std::vector<const SJClusterAlg*>::iterator SJAlgIter;

class SJClusterAlg {
private:
  Int_t N;
  std::vector<Float_t> *eta, *phi, *e, *mass, *pt, *numC;
  std::vector<Int_t> *ind;
  // these have to be vector pointers
  // that's how data is stored
  // they don't have to be deleted (no memory leak occurs)
  // For some ROOT hocus-pocus reason, these pointers have
  // to be initialized to zero

  const std::string alg;

  static std::vector<const SJClusterAlg*> AllAlgs;

  SJClusterAlg(TTree* tree, const char* algorithm);
  ~SJClusterAlg();

class sortByPt;

public:

  static const SJClusterAlg* AddAlg(TTree* tree, const char* algorithm);
  static void clean();

  const std::string& name() const;

  std::vector<TLorentzVector> jetsByPt(double pt_cut, double eta_cut) const;

  static SJAlgIter begin();
  static SJAlgIter   end();
};

#endif
