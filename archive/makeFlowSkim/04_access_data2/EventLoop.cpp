#include "EventLoop.h"
#include <stdexcept>
#include <iostream>

EventLoop::EventLoop() {
  // noting
}

void EventLoop::initialize() {
  // create an instance of the TChain class
  m_chain = new TChain(treeName);

  // loop through the input files and add them to the chain
  for (const auto& inputFile : inputFiles) {
    m_chain->Add(inputFile);
    std::cout << "Added files: " << inputFile << "\n";
  }

  // create a Data class instance
  m_data = new Data(m_chain);

  // initialize the algorithms
  for (const auto& algorithm : algorithms) {
    algorithm->initialize(m_data);
  }
}

void EventLoop::execute() {
  // sainty check. m_chain must not be zero
  if (!m_chain) {
    throw std::runtime_error("Calling execute() without initializiation.\n");
  }

  // start event loop
  Long64_t numEvents =m_chain->GetEntries();
  numEvents = 10000; 
  std::cout << "Tot evet: " << numEvents << "\n";
  for (Long64_t ir=0; ir < numEvents; ++ir) {
    // event number printout
    if (ir%100000==0&&ir!=0) {
      std::cout << "Event " << ir << "\n";
    }

    //read the data for i-th event
    m_chain->GetEntry(ir);

    // printout test
    // if (ir%100000==0) {
    //   std::cout << "Run number: " << m_data->runNb << "\n";
    // }

    // execute all the algorithms attached to this event loop
    for (const auto &algorithm : algorithms) {
      algorithm->execute();
    }
  }
  std::cout << "EventLoop done\n";
}
