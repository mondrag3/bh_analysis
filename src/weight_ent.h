#ifndef weight_ent_h
#define weight_ent_h

#include "TTree.h"

struct weight_ent {

  TTree *tree;

  Int_t num_fHt_fac;
  Float_t  *fHt_fac;
  Int_t num_fHt_ren;
  Float_t  *fHt_ren;

  Float_t *weights;
  Float_t w_pdf_up, w_pdf_cent, w_pdf_down;

  weight_ent();
  weight_ent(TTree* scales_tree, int s=0);
  ~weight_ent();
  weight_ent(const weight_ent& other);
  weight_ent& operator=(const weight_ent& other);

  bool operator==(const weight_ent& other) const;
  bool operator!=(const weight_ent& other) const;

  void set_tree(TTree* weights_tree);

  Float_t get_weight(Int_t fac, Int_t ren) const;

};

#endif
