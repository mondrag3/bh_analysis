#ifndef weight_hh
#define weight_hh

#include <string>
#include <vector>
#include <memory>

class TTree;

// Weights collector ************************************************

struct weight {
  std::string name;
  union {
    Double_t d;
    Float_t f;
  } w;
  bool is_float;
  weight(TTree *tree, const std::string& name, bool is_float=true);

  static std::vector<std::unique_ptr<const weight>> all;
  static void add(TTree* tree, const std::string& name, bool is_float=true) noexcept;
};

#endif
