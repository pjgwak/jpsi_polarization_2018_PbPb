#pragma once

#include "TString.h"
#include "TH1D.h"
#include "Data.h"

class Algorithm {
  public:
    /**
     * @brief construct a Alogrithm instance
     */
    Algorithm();

    /**
     * @brief initialize the algorithm
     * @param data - pointer to the instance of Data class
     */
    void initialize(Data* data);

    /**
     * @brief execute logics
     */
    void execute();

    /**
     * @brief pointer to the TH1D class instance
     */
    TH1D* h_runNb_m = nullptr;

  protected:
    /**
     * @brief instance of the Data class
     */
    Data* m_data = 0;

    /**
     * @brief
     */

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