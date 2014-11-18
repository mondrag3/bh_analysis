#include "rew_calc.h"

#include <iostream>
#include <memory>

using namespace std;

BHEvent event;

int main(int argc, char** argv)
{
  event.fac_scale = 50;

  auto FacHt4 = mk_fac_calc(new mu_fHt(0.25));
  auto RenHt2 = mk_ren_calc(new mu_fHt(0.5));
  auto FacDef = mk_fac_calc(new mu_fac_default());

  calc_all_scales();

  reweighter rew1(FacHt4,RenHt2);
  reweighter rew2(FacDef,RenHt2);

  rew1.stitch();
  rew2.stitch();

  return 0;
}
