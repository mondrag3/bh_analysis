//---------------------------
// Developer: Ivan Pogrebnyak
//---------------------------

#ifndef propmap_h
#define propmap_h

#include <string>
#include <set>
#include <array>
#include <unordered_map>
#include <sstream>
#include <stdexcept>

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

  virtual bool operator< (const prop& other) { return x<other.x; }
  virtual bool operator< (const T& other) { return x<other; }
  virtual bool operator==(const prop& other) { return x==other.x; }
  virtual bool operator==(const T& other) { return x==other; }

  virtual bool  l(const prop_base* other) const
  { return (x <  dynamic_cast< const prop<T>* >(other)->x); }
  virtual bool eq(const prop_base* other) const
  { return (x == dynamic_cast< const prop<T>* >(other)->x); }

  virtual prop_base* clone() const {
    return new prop<T>(*this);
  }
  virtual size_t hash() const {
    return std::hash<T>()(x);
  }
  virtual std::string str() const {
    std::stringstream ss;
    ss << x;
    return ss.str();
  }

  friend std::istream& operator>> (std::istream  &in, prop<T> &p) {
    in >> p.x;
    return in;
  }
  friend std::ostream& operator<< (std::istream  &out, const prop<T> &p) {
    out << p.x;
    return out;
  }
};
template<> std::string prop<std::string>::str() const {
  return x;
}

class prop_ptr {
protected:
  const prop_base *ptr;
  friend class std::hash<prop_ptr>;
  template <typename T, size_t N> friend class propmap;

  static std::unordered_map<const prop_base*,size_t> refcount;
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

  friend std::ostream& operator<< (std::ostream  &out, const prop_ptr &p) {
    out << p->str();
    return out;
  }
};
std::unordered_map<const prop_base*,size_t> prop_ptr::refcount;

namespace std {
  template<>
  struct hash<prop_ptr> {
    size_t operator()(const prop_ptr& p) const noexcept {
      return p.ptr->hash();
    }
  };

  template<typename T, size_t N>
  struct hash<array<T, N>> {
    size_t operator()(const array<T, N>& a) const noexcept {
      hash<T> hasher;
      size_t h = 0;
      for (size_t i = 0; i < N; ++i) h = h * 31 + hasher(a[i]);
        return h;
    }
  };
}

template <typename T, size_t N>
class propmap {
public:
  typedef std::array<prop_ptr,N> Key;
  typedef std::set<prop_ptr>     Set;

private:
  std::array<Set,N> _props;
  std::unordered_map<Key,T> _map;

public:
  void insert(const Key& key, const T& x) {
    for (size_t i=0;i<N;++i) {
      if (key[i].ptr) {
        auto& set = _props[i];
        if ( set.count(key[i])==0 ) set.insert(key[i]->clone());
      } else throw std::runtime_error("In propmap: NULL property");
    }
    _map[key] = x;
  }

  bool get(const Key& key, T& x) noexcept {
    if ( _map.count(key)==0 ) {
      return false;
    } else {
      x = _map[key];
      return true;
    }
  }

  T& get(const Key& key) {
    try { return _map.at(key);
    } catch (std::out_of_range& e) { throw std::out_of_range(
      std::string("In ") + __PRETTY_FUNCTION__ + ": " + e.what()
    ); }
  }

  template<size_t i> inline const Set& pset() noexcept {
    return std::get<i>(_props);
  }

};

#endif
