#include <iostream>
#include <memory>

#include <TFile.h>
#include <TTree.h>

#include "rew_calc.h"

using namespace std;

// definitions of externs
BHEvent event;

int main(int argc, char** argv)
{
  event.fac_scale = 50;

  usePDFset("CT10");

  auto FacHt4 = mk_fac_calc(new mu_fHt(0.25));
  auto RenHt2 = mk_ren_calc(new mu_fHt(0.5));
  auto FacDef = mk_fac_calc(new mu_fac_default());

  TFile *f = new TFile("out.root","recreate");
  TTree *tree = new TTree("weights","");

  reweighter rew1(FacHt4,RenHt2,tree);
  reweighter rew2(FacDef,RenHt2,tree);

  calc_all_scales();

  rew1.stitch();
  rew2.stitch();

  tree->Fill();

  f->Write();
  f->Close();
  delete f;

  return 0;
}
