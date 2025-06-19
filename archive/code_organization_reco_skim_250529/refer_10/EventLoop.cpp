#include "EventLoop.h"
#include <iostream>
#include <stdexcept>
// #include <vector>

EventLoop::EventLoop() {
  // skip
}

void EventLoop::initialize() {
  // create an instance of the TChain class
  m_chain = std::make_shared<TChain>(treeName);

  // add input files
  for (const auto& inputFile : inputFiles) {
    m_chain->Add(inputFile);
    std::cout << "Add file: " << inputFile << std::endl;
  }

  // create instance of the Data class
  m_data = std::make_shared<Data>(m_chain);

  // initialize the algorithms
  for (const auto& algorithm : algorithms) {
    algorithm->initialize(m_data);
  }
}

void EventLoop::execute() {
  // sanity check.
  if (!m_chain) {
    throw std::runtime_error("Calling execute while the event loop was not initialized.");
  }

  // start event loop
  for (Long64_t i = 0; i < m_chain->GetEntries(); ++i) {
    // event number printout
    if (i%1000==0) {
      std::cout << "Event: " << i << std::endl;
    }

    // read the data for i-th event
    m_chain->GetEntry(i);

    // execute all algorithms attached to this event loop
    for (const auto& algorithm : algorithms) {
      algorithm->execute();
    }
  }
}