#include "rew_calc.h"

#include <iostream>
#include <memory>

using namespace std;

BHEvent event;

int main(int argc, char** argv)
{
  event.fac_scale = 50;

  unique_ptr<mu_fcn> fd(new mu_fac_default());
  unique_ptr<mu_fcn> rd(new mu_ren_default());

  unique_ptr<mu_fcn> Ht4(new mu_fHt(0.25));
  unique_ptr<mu_fcn> Ht2(new mu_fHt(0.5));

  cout << "Event default fac_scale = " << fd->mu() << endl;

  auto FacHt4 = mk_fac_calc(Ht4.get());
  auto RenHt2 = mk_ren_calc(Ht2.get());
  auto FacDef = mk_fac_calc(fd.get());

  calc_all_scales();

  return 0;
}
