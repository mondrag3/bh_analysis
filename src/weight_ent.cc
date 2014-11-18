#include "weight_ent.h"

#include <iostream>

using namespace std;

weight_ent::weight_ent()
: tree(NULL), fHt_fac(NULL), fHt_ren(NULL), weights(NULL)
{ }

weight_ent::weight_ent(TTree* scales_tree, int s)
: tree(NULL), fHt_fac(NULL), fHt_ren(NULL), weights(NULL)
{
  scales_tree->SetBranchAddress("num_fHt_fac", &num_fHt_fac);
  scales_tree->SetBranchAddress("num_fHt_ren", &num_fHt_ren);

  scales_tree->GetEntry(s);

  fHt_fac = new Float_t[num_fHt_fac];
  fHt_ren = new Float_t[num_fHt_ren];

  scales_tree->SetBranchAddress("fHt_fac", fHt_fac);
  scales_tree->SetBranchAddress("fHt_ren", fHt_ren);

  scales_tree->GetEntry(s);

  // TEST: print header
  cout << "fHt_fac:";
  for (Int_t i=0;i<num_fHt_fac;++i) cout <<' '<< fHt_fac[i];
  cout << endl;
  cout << "fHt_ren:";
  for (Int_t i=0;i<num_fHt_ren;++i) cout <<' '<< fHt_ren[i];
  cout << endl;
}

void weight_ent::set_tree(TTree* weights_tree) {
  weights = new Float_t[num_fHt_fac*num_fHt_ren];

  tree = weights_tree;

  tree->SetBranchAddress("pdf_up",&w_pdf_up);
  tree->SetBranchAddress("pdf_cent",&w_pdf_cent);
  tree->SetBranchAddress("pdf_down",&w_pdf_down);
  tree->SetBranchAddress("weights",weights);
}

Float_t weight_ent::get_weight(Int_t fac, Int_t ren) const {
  return weights[fac*num_fHt_ren+ren];
}

weight_ent::~weight_ent() {
  delete[] fHt_fac;
  delete[] fHt_ren;
  delete[] weights;
}

weight_ent::weight_ent(const weight_ent& other)
: num_fHt_fac(other.num_fHt_fac),
  fHt_fac(new Float_t[num_fHt_fac]),
  num_fHt_ren(other.num_fHt_ren),
  fHt_ren(new Float_t[num_fHt_ren]),
  weights(new Float_t[num_fHt_fac*num_fHt_ren]),
  w_pdf_up  (other.w_pdf_up),
  w_pdf_cent(other.w_pdf_cent),
  w_pdf_down(other.w_pdf_down)
{
  for (Int_t i=0;i<num_fHt_fac;++i)
    fHt_fac[i] = other.fHt_fac[i];
  for (Int_t i=0;i<num_fHt_ren;++i)
    fHt_ren[i] = other.fHt_ren[i];
  for (Int_t i=0,n=num_fHt_fac*num_fHt_ren;i<n;++i)
    weights[i] = other.weights[i];
}

weight_ent& weight_ent::operator=(const weight_ent& other) {
  delete[] fHt_fac;
  delete[] fHt_ren;
  delete[] weights;

  num_fHt_fac = other.num_fHt_fac;
  fHt_fac = new Float_t[num_fHt_fac];
  num_fHt_ren = other.num_fHt_ren;
  fHt_ren = new Float_t[num_fHt_ren];
  weights = new Float_t[num_fHt_fac*num_fHt_ren];
  w_pdf_up   = other.w_pdf_up;
  w_pdf_cent = other.w_pdf_cent;
  w_pdf_down = other.w_pdf_down;

  for (Int_t i=0;i<num_fHt_fac;++i)
    fHt_fac[i] = other.fHt_fac[i];
  for (Int_t i=0;i<num_fHt_ren;++i)
    fHt_ren[i] = other.fHt_ren[i];
  for (Int_t i=0,n=num_fHt_fac*num_fHt_ren;i<n;++i)
    weights[i] = other.weights[i];

  return *this;
}

bool weight_ent::operator==(const weight_ent& other) const {
  if (num_fHt_fac != other.num_fHt_fac) return false;
  if (num_fHt_ren != other.num_fHt_ren) return false;

  for (Int_t i=0;i<num_fHt_fac;++i)
    if (fHt_fac[i] != other.fHt_fac[i]) return false;
  for (Int_t i=0;i<num_fHt_ren;++i)
    if (fHt_ren[i] != other.fHt_ren[i]) return false;

  return true;
}

bool weight_ent::operator!=(const weight_ent& other) const {
  return !( operator==(other) );
}
