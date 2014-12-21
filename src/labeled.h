//---------------------------
// Developer: Ivan Pogrebnyak
//---------------------------

#ifndef labeled_h
#define labeled_h

#include <string>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <set>
#include <stdexcept>

template<class T>
class labeled {
  typedef T               val_t;
  typedef std::string     lbl_t;
  typedef std::set<lbl_t> set_t;
  typedef std::vector<T>  container_t;
  typedef std::unordered_map<lbl_t,container_t> map_t;

  set_t _labels;
  map_t _values;

  void add(const std::string& str) {
    static labeled<T>::lbl_t label;
    static labeled<T>::val_t value;

    const size_t sep = str.find(':');
    if (sep==std::string::npos) {
      label = std::string();
      std::stringstream(str) >> value;
    } else {
      std::stringstream(str.substr(0,sep)) >> label;
      std::stringstream(str.substr(sep+1)) >> value;
    }
    _labels.insert(label);
    _values[label].push_back(value);
  }

public:
  labeled(const std::vector<std::string>& container) {
    for (const std::string& str : container) add(str);
  }

  const set_t& labels() const noexcept { return _labels; }
  const container_t& values(const lbl_t& label) const {
    try {
      return _values.at(label);
    } catch (std::out_of_range& e) {
      throw std::out_of_range(
        std::string("In ") + __PRETTY_FUNCTION__ + ": " + e.what()
      );
    }
  }

  friend std::istream& operator>> (std::istream& is, labeled<T>& x) {
    static std::string str;
    is >> str;
    x.add(str);
    return is;
  }
};

#endif
