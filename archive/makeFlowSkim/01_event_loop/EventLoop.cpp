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
}

void EventLoop::execute() {
  // sainty check. m_chain must not be zero
  if (!m_chain) {
    throw std::runtime_error("Calling execute() without initializiation.\n");
  }

  // start event loop
  Long64_t numEvents =m_chain->GetEntries();
  numEvents = 10000; 
  for (Long64_t ir=0; ir < numEvents; ++ir) {
    // event number printout
    if (ir%1000==0) {
      std::cout << "Event " << ir << "\n";
    }

    //read the data for i-th event
    m_chain->GetEntry(ir);
  }
  std::cout << "EventLoop done\n";
}