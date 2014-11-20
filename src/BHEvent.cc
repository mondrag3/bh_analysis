#include "BHEvent.h"

#include <cmath>

#include <TTree.h>

void BHEvent::SetTree(TTree* tree) {
  Init(tree);
}

void BHEvent::SetTree(TTree* tree, BHEvent::select_t branches) {
  Init(tree);
  TurnAllOff();
  TurnOn(branches);
}

void BHEvent::Init(TTree* tree) {
  this->tree = tree;

  tree->SetBranchAddress("id", &eid);
  tree->SetBranchAddress("nparticle", &nparticle);
  tree->SetBranchAddress("px", px);
  tree->SetBranchAddress("py", py);
  tree->SetBranchAddress("pz", pz);
  tree->SetBranchAddress("E", E);
  tree->SetBranchAddress("alphas", &alphas);
  tree->SetBranchAddress("kf", kf);
  tree->SetBranchAddress("weight", &weight);
  tree->SetBranchAddress("weight2", &weight2);
  tree->SetBranchAddress("me_wgt", &me_wgt);
  tree->SetBranchAddress("me_wgt2", &me_wgt2);
  tree->SetBranchAddress("x1", &x[0]);
  tree->SetBranchAddress("x2", &x[1]);
  tree->SetBranchAddress("x1p", &xp[0]);
  tree->SetBranchAddress("x2p", &xp[1]);
  tree->SetBranchAddress("id1", &id[0]);
  tree->SetBranchAddress("id2", &id[1]);
  tree->SetBranchAddress("fac_scale", &fac_scale);
  tree->SetBranchAddress("ren_scale", &ren_scale);
  tree->SetBranchAddress("nuwgt", &nuwgt);
  tree->SetBranchAddress("usr_wgts", usr_wgts);
  tree->SetBranchAddress("alphasPower", &alphas_power);
  tree->SetBranchAddress("part", part);
}

void BHEvent::TurnOn(const char* branch) {
  tree->SetBranchStatus(branch,1);
}

void BHEvent::TurnOn(BHEvent::select_t branches) {
  switch (branches) {
    case kinematics: {

      tree->SetBranchStatus("nparticle",1);
      tree->SetBranchStatus("px",1);
      tree->SetBranchStatus("py",1);
      tree->SetBranchStatus("pz",1);
      tree->SetBranchStatus("E",1);
      tree->SetBranchStatus("kf",1);
      tree->SetBranchStatus("weight",1);

    } break;
    case reweighting: {

      tree->SetBranchStatus("nparticle",1); // for px & py
      tree->SetBranchStatus("px",1); // for Ht
      tree->SetBranchStatus("py",1); // for Ht
      tree->SetBranchStatus("weight",1);
      tree->SetBranchStatus("weight2",1);
      tree->SetBranchStatus("me_wgt",1);
      tree->SetBranchStatus("me_wgt2",1);
      tree->SetBranchStatus("x1",1);
      tree->SetBranchStatus("x2",1);
      tree->SetBranchStatus("x1p",1);
      tree->SetBranchStatus("x2p",1);
      tree->SetBranchStatus("id1",1);
      tree->SetBranchStatus("id2",1);
      tree->SetBranchStatus("fac_scale",1);
      tree->SetBranchStatus("ren_scale",1);
      tree->SetBranchStatus("usr_wgts",1);
      tree->SetBranchStatus("part",1);
      tree->SetBranchStatus("alphasPower",1);
      tree->SetBranchStatus("alphas",1);

    } break;
  }
}

void BHEvent::TurnAllOff() {
  tree->SetBranchStatus("*",0);
};
void BHEvent::TurnOff(const char* branch) {
  tree->SetBranchStatus(branch,0);
};

void BHEvent::SetPart(Char_t part) { this->part[0] = part; }
void BHEvent::SetAlphasPower(Char_t n) { this->alphas_power = n; }

template<typename T> T sq(T x) { return x*x; }

Double_t BHEvent::Ht() const {
  Double_t _Ht = 0.;
  for (Int_t i=0;i<nparticle;++i) {
    _Ht += sqrt( sq(px[i]) + sq(py[i]) ); // <-- Need to make sure that this is really Ht
  }
  return _Ht;
}
