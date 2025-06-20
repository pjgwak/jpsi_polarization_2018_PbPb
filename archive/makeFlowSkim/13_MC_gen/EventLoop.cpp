#include "EventLoop.h"
#include <stdexcept>
#include <iostream>

EventLoop::EventLoop() {
  // nothing
}

void EventLoop::setNumEvent() {
  Long64_t nEvt_temp = m_chain->GetEntries();
  if (nEvent == -1)
    nEvent = nEvt_temp;
  else if (nEvt_temp < nEvent) 
    nEvent = nEvt_temp;
  // Otherwise - don't change nEvent
}

void EventLoop::initialize()
{
  std::cout << "\n\n\n===== Initialize the event loop =====\n";

  m_chain = new TChain(treeName);
  for (const auto &inputFile : inputFiles) {
    m_chain->Add(inputFile);
    std::cout << "Added files: " << inputFile << "\n";
  }

  // create a Data class instance
  if (useRun2)
    m_data = new DataRun2(m_chain, isMC, isGenOnly, isEP, true); // isRun2_flag=true
  else
    m_data = new DataRun3(m_chain, isMC, isGenOnly, isEP, false); // Data run3

  // initialize the algorithms
  for (const auto &algorithm : algorithms)
    algorithm->initialize(m_data, useRun2);

  setNumEvent();
}

void EventLoop::execute()
{
  // sainty check. m_chain must not be zero
  if (!m_chain)
    throw std::runtime_error("Calling execute() without initializiation.\n");

  // start event loop
  std::cout << "Tot evet: " << nEvent << "\n";
  for (Long64_t ir = 0; ir < nEvent; ++ir)
  {
    // event number printout
    if (ir % 100000 == 0 && ir != 0)
      std::cout << "Event " << ir << "\n";
    
    // read the data for i-th event
    m_chain->GetEntry(ir);

    // prepare variables
    Long64_t recoQQSize = m_data->getRecoQQSize();
    
    // ===== reco loop ===== //
    for (Long64_t irqq = 0; irqq < recoQQSize; ++irqq) {
      for (const auto &algorithm : algorithms)
        algorithm->execute(irqq);
    } // reco loop end

    // fill the outTree
    for (const auto &algorithm : algorithms)
      algorithm->fillTree();
  } // end of event loop
  std::cout << "EventLoop done\n";
}

void EventLoop::executeGen()
{
  if (!m_chain)
    throw std::runtime_error("Calling execute() without initializiation.\n");

  // start event loop
  std::cout << "Tot evet: " << nEvent << "\n";
  for (Long64_t ir = 0; ir < nEvent; ++ir)
  {
    // event number printout
    if (ir % 100000 == 0 && ir != 0)
      std::cout << "Event " << ir << "\n";

    m_chain->GetEntry(ir);

    // prepare variables
    Long64_t genQQSize = m_data->getGenQQSize();

    // ===== gen loop ===== //
    for (Long64_t igqq = 0; igqq < genQQSize; ++igqq) {
      for (const auto &algorithm : algorithms)
        algorithm->executeGen(igqq);
    } // gen loop end

    // fill the outTree
    for (const auto &algorithm : algorithms)
      algorithm->fillTree();
  } // end of event loop
  std::cout << "EventLoop done\n";
}