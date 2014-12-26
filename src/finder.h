#ifndef finder_h
#define finder_h

#include <string>
#include <sstream>
#include <stdexcept>

template <class T>
class finder {
private:
  const T * const arr;
  const size_t n;
public:
  finder(const T* arr, size_t n)
  : arr(arr), n(n) { }
  size_t operator ()(const T& element, size_t occurance=1) {
    size_t cnt = 0;
    for (size_t i=0;i<n;++i) {
      if (arr[i]==element) ++cnt;
      else continue;
      if (cnt==occurance) return i;
    }
    std::stringstream ss;
    ss << "Event doesn't have ";
    occurance==1 ? ss << "any" : ss << occurance;
    ss << " particles with PDG id " << element;
    throw std::runtime_error(ss.str());
  }
};

#endif
