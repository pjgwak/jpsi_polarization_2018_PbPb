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
    void execute(Long64_t irqq);

    /**
     * @brief pointer to the TH1D class instance
     */
    TH1D* h_runNb_m = nullptr;

    /**
     * @brief selection cuts
     */
    bool cut_selectSf = false;

  protected:
    /**
     * @brief apply the selection
     */
    bool passedSelection(Long64_t irqq);
    bool passedSoftMuonCut(Long64_t irqq);
    bool passedMuonAcc2018();
    bool checkMuonAcc2018(double muPt, double muEta);

    /**
     * @brief fill the histograms
     */
    void fillHists();    

    /**
     * @brief instance of the Data class
     */
    Data* m_data = nullptr;



    /**
     * @brief
     */

    /**
     * @brief
     */

    
};