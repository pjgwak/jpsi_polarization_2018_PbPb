#include "Data.h"

Data::Data(TTree* tree) : m_tree(tree) {
  // m_tree->SetBranchAddress("XXX", &XXX);
  m_tree->SetBranchAddress("runNb", &runNb);
}