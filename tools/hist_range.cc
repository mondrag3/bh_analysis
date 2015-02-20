#include "hist_range.hh"

#include <iostream>
#include <limits>
#include <tuple>

#include <TH1.h>

using namespace std;

// constructor
hist_range::hist_range()
: min(numeric_limits<double>::max()), max(-min), applied(false) { }

// add values
void hist_range::operator()(double min, double max) noexcept {
  if (applied) {
    cerr << "This range has been applied" << endl;
    return;
  }
  if (min < this->min) this->min = min;
  if (max > this->max) this->max = max;
}

// apply to a histogram
TH1* hist_range::operator()(TH1* h) noexcept {
  if (!applied) {
    applied = true;
    bool both = false;
    if (min > 0.) {
      if (min/max < 0.25) {
        min = 0.;
        max /= 0.95;
      } else both = true;
    } else if (max < 0.) {
      if (min/max < 0.25) {
        max = 0.;
        min /= 0.95;
      } else both = true;
    } else if (min==0.) {
      max /= 0.95;
    } else if (max==0.) {
      min /= 0.95;
    } else both = true;
    if (both) {
      tie(min,max) =
      forward_as_tuple(1.05556*min - 0.05556*max, 1.05556*max - 0.05556*min);
    }
  }
  
  h->SetMinimum(min);
  h->SetMaximum(max);
  
  return h;
}
