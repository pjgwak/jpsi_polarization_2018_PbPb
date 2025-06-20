#pragma once

#include <vector>
#include <memory>
#include "TString.h"
#include "TChain.h"
#include "Data.h"
#include "Algorithm.h"

class EventLoop {
  public:
    /**
    * @brief Construct a new EventLoop object
    */
    EventLoop();

    /**
     * @brief Initialize the event loop
     */
    void initialize();

    /**
     * @brief Execute the event loop
     */
    void execute();

    /**
     * @brief list of input ROOT file names
     */
    std::vector<TString> inputFiles;

    /**
     * @brief Name of the TTree instance. Must be same in all files
     */
    TString treeName;

    /**
     * @brief List of algorithms to be executed in the event loop
     */
    // std::vector<std::shared_ptr<Algorithm>> algorithms;
    std::vector<Algorithm*> algorithms;

  protected:
    /**
     * @brief Instance of the TChain class used to read the data
     */
    std::shared_ptr<TChain> m_chain = nullptr;

    /**
     * @brief Instance of the data-access class
     */
    std::shared_ptr<Data> m_data = nullptr;
};