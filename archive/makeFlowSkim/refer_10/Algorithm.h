#pragma once

#include <memory>
#include "TString.h"
#include "TH1D.h"
#include "Data.h"
#include "Cut.h"


class Algorithm {
  public:
    /**
     * @brief Construct a new EventLoop object
     */
    Algorithm();

    /**
     * @brief Initialize the algorithm
     * @param data - pointer to the data-access class instance
     */
    void initialize(std::shared_ptr<Data> data);

    /**
     * @brief Execute.
     */
    void execute();

    /**
     * @brief Pointer to the histograms class instance defined as public attribute
     */
    TH1D *h_ditau_m = nullptr;
    TH1D *h_ditau_visM = nullptr;

    /**
     * @brief Selection cuts
     */
    bool cut_selectSF = false;
    bool cut_selectDF = false;

    /**
     * @brief Selection cuts implemented using the Cut class
     */
    Cut cut_maxVisMass;

    /** 
     * @brief Other algorithm properties
     */
    bool p_isMC = false;
    double p_sampleWeight = 1.;

  protected:
    /**
     * @brief Apply the selection
     * @return true when event passed the selection
     */
    bool passedSelection();

    /**
     * @brief Evaluate the event weight
     */
    void calculateWeight();

    /**
     * @brief Fill the histograms
     */
    void fillPlots();

    /**
     * @brief Instance of the data-access class
     */
    std::shared_ptr<Data> m_data = nullptr;

    /**
     * @brief Total event weight used when filling histograms
     */
    double m_eventWeight;
};