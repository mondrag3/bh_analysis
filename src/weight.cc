#include "weight.h"

#include <TTree.h>

using namespace std;

// Weights collector ************************************************
weight::weight(TTree *tree, const string& name): name(name) {
  if ( tree->SetBranchAddress(name.c_str(), &w) == TTree::kMissingBranch )
    exit(1);
}

void weight::add(TTree* tree, const string& name) {
  all.emplace_back( new weight(tree,name) );
}

std::vector<std::unique_ptr<const weight>> weight::all;
