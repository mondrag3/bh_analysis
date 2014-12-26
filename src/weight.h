#ifndef filler_h
#define filler_h

#include <string>
#include <vector>
#include <memory>

#include <Rtypes.h>

class TTree;

// Weights collector ************************************************
struct weight {
  std::string name;
  Float_t w;
  weight(TTree *tree, const std::string& name);

  static std::vector<std::unique_ptr<const weight>> all;

  static void add(TTree* tree, const std::string& name);
};

#endif
