#ifndef hist_range_hh
#define hist_range_hh

class TH1;

class hist_range {
  double min, max;
  bool applied;
public:
  hist_range();
  // add values
  void operator()(double min, double max) noexcept;
  // apply to a histogram
  TH1* operator()(TH1* h) noexcept;
};

#endif
