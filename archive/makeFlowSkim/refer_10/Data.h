#pragma once

#include "TTree.h"
#include "TLorentzVector.h"
#include <memory>

class Data {
  public:
    /**
     * @brief Construct a new Data object
     * @param tree - pointer to the TTree (or TChain) class
     */
    Data(std::shared_ptr<TTree> tree);

    /**
     * @brief Tree variables
     */
    Float_t ditau_mmc_mlm_m;
    UInt_t tau_0;
    UInt_t tau_1;
    TLorentzVector* tau_0_p4 = nullptr; // Don' need and have to use smart pinter.
    TLorentzVector *tau_1_p4 = nullptr; // when handle SetBranch
    Double_t weight_mc;

  protected:
    /**
     * @brief pointer to the TTree (or TChain) class
     */
    std::shared_ptr<TTree> m_tree = nullptr;
};