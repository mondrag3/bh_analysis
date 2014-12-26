#include "selector_base.h"

class hist_defs;

class selector: public selector_base {
  hist_defs *h;

public:
  selector(timed_counter *counter);
  virtual ~selector();

  virtual void SlaveBegin();
  virtual Bool_t Process(Long64_t entry);
};
