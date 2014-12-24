#include "BHEvent.h"

#include <cmath>

#include <TTree.h>

void BHEvent::SetTree(TTree* tree, select_t branches, bool old) {
  this->tree = tree;

  switch (branches) {
    case kinematics: {

      tree->SetBranchAddress("nparticle", &nparticle);
      tree->SetBranchAddress("px", px);
      tree->SetBranchAddress("py", py);
      tree->SetBranchAddress("pz", pz);
      tree->SetBranchAddress("E" , E );
      tree->SetBranchAddress("kf", kf);

    } break;
    case reweighting: {

      tree->SetBranchAddress("nparticle", &nparticle);
      tree->SetBranchAddress("px", px);
      tree->SetBranchAddress("py", py);
      tree->SetBranchAddress("alphas", &alphas);
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
      tree->SetBranchAddress("usr_wgts", usr_wgts);
      if (!old) {
        tree->SetBranchAddress("alphasPower", &alphas_power);
        tree->SetBranchAddress("part", part);
      }

    } break;
    case cross_section: {

      // tree->SetBranchAddress("nparticle", &nparticle);
      // tree->SetBranchAddress("kf", kf);
      tree->SetBranchAddress("weight", &weight);

    } break;
    default: {

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
  }

}

void BHEvent::SetPart(Char_t part) { this->part[0] = part; }
void BHEvent::SetAlphasPower(Char_t n) { this->alphas_power = n; }

template<typename T> T sq(T x) { return x*x; }

Double_t BHEvent::Ht() const noexcept {
  Double_t _Ht = 0.;
  for (Int_t i=0;i<nparticle;++i)
    _Ht += sqrt( sq(px[i]) + sq(py[i]) );
  return _Ht;
}
Double_t BHEvent::Ht_Higgs() const noexcept {
  Double_t _Ht = 0.;
  for (Int_t i=0;i<nparticle;++i) {
    if (kf[i]==25) _Ht += sqrt( sq(E[i]) - sq(pz[i]) ); // sqrt(m^2+pt^2)
    else _Ht += sqrt( sq(px[i]) + sq(py[i]) ); // pt
  }
  return _Ht;
}
