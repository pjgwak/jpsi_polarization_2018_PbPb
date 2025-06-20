#include "Data.h"

Data::Data(std::shared_ptr<TTree> tree) : m_tree(tree) {
  m_tree->SetBranchAddress("ditau_mmc_mlm_m", &ditau_mmc_mlm_m);
  m_tree->SetBranchAddress("tau_0", &tau_0);
  m_tree->SetBranchAddress("tau_1", &tau_1);
  m_tree->SetBranchAddress("tau_0_p4", &tau_0_p4);
  m_tree->SetBranchAddress("tau_1_p4", &tau_1_p4);
  m_tree->SetBranchAddress("weight_mc", &weight_mc);
}