# pragma once

#include "TTree.h"

class Data {
  public:
    /**
     * @brief construct a new Data object
     */
    Data(TTree* tree);

    /**
     * @brief tree variables
     */
    UInt_t runNb;

  protected:
    /**
     * @brief pointer to the TTree (or TChain) class
     */
    TTree* m_tree = nullptr;

    

    /**
     * @brief
     */

    /**
     * @brief
     */

    /**
     * @brief
     */


};