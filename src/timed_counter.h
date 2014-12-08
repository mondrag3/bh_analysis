#ifndef timed_counter_h
#define timed_counter_h

#include <ctime>
#include <iostream>
#include <iomanip>

class timed_counter {
  unsigned short seconds = 0, minutes = 0, hours = 0;
  time_t last_time, cur_time;
public:
  timed_counter()
  : seconds(0), minutes(0), hours(0), last_time( time(0) )
  { }
  template<typename T>
  void prt(const T& ent) const {
    using namespace std;

    cout << setw(10) << ent << " | ";
    if (hours) {
      cout << setw(5) << hours   << ':'
      << setw(2) << minutes << ':'
      << setw(2) << seconds;
    } else if (minutes) {
      cout << setw(2) << minutes << ':'
      << setw(2) << seconds;
    } else {
      cout << setw(2) << seconds <<'s';
    }
  }
  template<typename T>
  void operator()(const T& ent) {
    using namespace std;
    
    cur_time = time(0);
    if ( difftime(cur_time,last_time) > 1 ) {
      prt(ent);
      cout.flush();

      if (hours)        for (char i=0;i<24;i++) cout << '\b';
      else if (minutes) for (char i=0;i<18;i++) cout << '\b';
      else              for (char i=0;i<16;i++) cout << '\b';

      last_time = cur_time;
      ++seconds;
      if (seconds==60) {
        seconds = 0;
        ++minutes;
        if (minutes==60) {
          minutes = 0;
          ++hours;
        }
      }

    }
  }
};

#endif
