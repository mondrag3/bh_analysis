#ifndef timed_counter_h
#define timed_counter_h

#include <ctime>

class timed_counter {
  unsigned short seconds, minutes, hours;
  time_t last_time, cur_time;
  bool newline;
  typedef long num_t;

public:
  timed_counter(bool newline=false);
  void prt(const num_t& ent) const noexcept;
  void operator()(const num_t& ent) noexcept;
};

#endif
