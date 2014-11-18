//---------------------------
// Developer: Ivan Pogrebnyak
//---------------------------

#ifndef propmap_h
#define propmap_h

#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <stdexcept>

//#include <iostream>

#include <boost/unordered_map.hpp>

class prop_base {
public:
  virtual ~prop_base() { }
  virtual bool  l(const prop_base* other) const =0;
  virtual bool eq(const prop_base* other) const =0;
  virtual prop_base* clone() const =0;
  virtual size_t hash() const =0;
  virtual std::string str() const =0;
};

template<typename T>
class prop: public prop_base {
protected:
  T x;
public:
  prop(const prop& other): x(other.x) { }
  prop(const T& other): x(other) { }
  virtual prop& operator=(const prop& other) { x=other.x; return *this; }
  virtual prop& operator=(const T& other) { x=other; return *this; }
  virtual ~prop() { }
  virtual bool  l(const prop_base* other) const
  { return (x <  dynamic_cast< const prop<T>* >(other)->x); }
  virtual bool eq(const prop_base* other) const
  { return (x == dynamic_cast< const prop<T>* >(other)->x); }
  virtual prop_base* clone() const {
    return new prop<T>(*this);
  }
  virtual size_t hash() const {
    return boost::hash_value(x);
  }
  virtual std::string str() const {
    std::stringstream ss;
    ss << x;
    return ss.str();
  }
};
template<> std::string prop<std::string>::str() const {
  return x;
}

class prop_ptr;
size_t hash_value(const prop_ptr& x);

class prop_ptr {
protected:
  const prop_base *ptr;
  friend size_t hash_value(const prop_ptr& x);

  static boost::unordered_map<const prop_base*,size_t> refcount;
  void addref(const prop_base* ptr) {
    if ( refcount.count(ptr)==0 ) refcount[ptr] = 1;
    else ++refcount[ptr];
  }
  void release(const prop_base* ptr) {
    size_t& cnt = refcount[ptr];
    if (cnt==1) delete ptr;
    --cnt;
  }

public:
  prop_ptr(): ptr(NULL) { }
  prop_ptr(const prop_base* other): ptr(other) {
    addref(ptr);
  }
  prop_ptr(const prop_ptr& other): ptr(other.ptr) {
    addref(ptr);
  }
  ~prop_ptr() { release(ptr); }
  prop_ptr& operator= (const prop_base* other) {
    release(ptr);
    ptr = other;
    addref(ptr);
    return *this;
  }
  prop_ptr& operator= (const prop_ptr& other) {
    release(ptr);
    ptr = other.ptr;
    addref(ptr);
    return *this;
  }
  bool operator< (const prop_ptr& other) const {
    return ( ptr->l(other.ptr) );
  }
  bool operator==(const prop_ptr& other) const {
    return ( ptr->eq(other.ptr) );
  }
  const prop_base& operator*() const { return *ptr; }
  const prop_base* operator->() const { return ptr; }
};
boost::unordered_map<const prop_base*,size_t> prop_ptr::refcount;

size_t hash_value(const prop_ptr& x) {
  return x.ptr->hash();
}

typedef std::set<prop_ptr>::iterator propset_iter;

template <typename T>
class propmap {
  size_t nprop;

  typedef std::set<prop_ptr> prop_set;
  std::vector<prop_set> _props;
  boost::unordered_map< std::vector<prop_ptr>, T > _map;

public:
  propmap(size_t nprop): nprop(nprop), _props(nprop) { }

  void insert(const std::vector<prop_ptr>& key, const T& x) {
    if (key.size()==nprop) {
/*
      std::cout << "+ " << x << ' '
                << key[0]->str() << ' ' << key[1]->str() << std::endl;
*/
      for (size_t i=0;i<nprop;++i) {
        prop_set& set = _props[i];
        if ( set.count(key[i])==0 ) {
          set.insert(key[i]->clone());
        }
      }
      _map[key] = x;

    } else throw std::runtime_error("In propmap: Wrong length key");
  }

  bool get(const std::vector<prop_ptr>& key, T& x) {
    if ( _map.count(key)==0 ) {
      return false;
    } else {
      x = _map[key];
      return true;
    }
  }

  propset_iter begin(size_t i) { return _props.at(i).begin(); }
  propset_iter   end(size_t i) { return _props.at(i).end();   }

};

#define pmloop(propmap_var,it,i) \
  for (propset_iter it=propmap_var.begin(i), \
    it##_end=propmap_var.end(i); it!=it##_end; ++it)

#endif
