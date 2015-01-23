#include "weight.hh"

#include <TTree.h>

using namespace std;

weight::weight(TTree *tree, const string& name, bool is_float=true)
: name(name), is_float(is_float)
{
  TBranch* const br = tree->GetBranch(name.c_str());
  if (br) {
    if (is_float) br->SetAddress(&w.f);
    else br->SetAddress(&w.d);
  } else exit(1);
}

void weight::add(TTree* tree, const string& name, bool is_float) noexcept {
  all.emplace_back( new weight(tree,name,is_float) );
}

vector<unique_ptr<const weight>> weight::all;
