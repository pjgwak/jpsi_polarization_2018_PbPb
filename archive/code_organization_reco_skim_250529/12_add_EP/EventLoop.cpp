#include "EventLoop.h"
#include <stdexcept>
#include <iostream>

EventLoop::EventLoop() {
  // noting
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

  // create an instance of the TChain class
  m_chain = new TChain(treeName);

  // loop through the input files and add them to the chain
  for (const auto &inputFile : inputFiles)
  {
    m_chain->Add(inputFile);
    std::cout << "Added files: " << inputFile << "\n";
  }

  // create a Data class instance
  if (useRun2)
    m_data = new DataRun2(m_chain, isMC, isEP);
  else
    m_data = new DataRun3(m_chain, isMC, isEP);

  // initialize the algorithms
  for (const auto &algorithm : algorithms)
    algorithm->initialize(m_data, useRun2);

  // set a number of events to use
  setNumEvent();
}

void EventLoop::execute()
{
  // sainty check. m_chain must not be zero
  if (!m_chain)
  {
    throw std::runtime_error("Calling execute() without initializiation.\n");
  }

  // start event loop
  std::cout << "Tot evet: " << nEvent << "\n";
  for (Long64_t ir = 0; ir < nEvent; ++ir)
  {
    // event number printout
    if (ir % 100000 == 0 && ir != 0)
    {
      std::cout << "Event " << ir << "\n";
    }

    // read the data for i-th event
    m_chain->GetEntry(ir);

    // printout test
    // if (ir%100000==0) {
    //   std::cout << "Run number: " << m_data->runNb << "\n";
    // }

    // ===== reco loop ===== //
    Long64_t recoQQSize = 0; // may be we can make it as function returning Long64_t
    if (useRun2) {
        DataRun2* data_run2 = dynamic_cast<DataRun2*>(m_data);
        recoQQSize = static_cast<Long64_t>(data_run2->Reco_QQ_size);
    } else {
        DataRun3* data_run3 = dynamic_cast<DataRun3*>(m_data);
        recoQQSize = static_cast<Long64_t>(data_run3->Reco_QQ_size);
    }

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


// ====== Childern classes ===== //
// void EventLoopRun2::initialize()
// {
//   std::cout << "\n\n\n===== Initialize the event loop =====\n";

//   // create an instance of the TChain class
//   m_chain = new TChain(treeName);

//   // loop through the input files and add them to the chain
//   for (const auto &inputFile : inputFiles) {
//     m_chain->Add(inputFile);
//     std::cout << "Added files: " << inputFile << "\n";
//   }

//   // create a Data class instance
//   m_data = new DataRun2(m_chain, isMC);

//   // initialize the algorithms
//   for (const auto &algorithm : algorithms)
//     algorithm->initialize(m_data);

//   // set a number of events to use
//   setNumEvent();
// }

// void EventLoopRun2::execute()
// {
//   // sainty check. m_chain must not be zero
//   if (!m_chain)
//   {
//     throw std::runtime_error("Calling execute() without initializiation.\n");
//   }

//   // start event loop
//   std::cout << "Tot evet: " << nEvent << "\n";
//   for (Long64_t ir = 0; ir < nEvent; ++ir)
//   {
//     // event number printout
//     if (ir % 100000 == 0 && ir != 0) {
//       std::cout << "Event " << ir << "\n";
//     }

//     // read the data for i-th event
//     m_chain->GetEntry(ir);

//     // printout test
//     // if (ir%100000==0) {
//     //   std::cout << "Run number: " << m_data->runNb << "\n";
//     // }

//     // reco loop
//     for (Long64_t irqq = 0; irqq < m_data->Reco_QQ_size; ++irqq) {
//       for (const auto &algorithm : algorithms) {
//         algorithm->execute(irqq);
//       }
//     } // reco loop end

//     // fill the outTree
//     for (const auto &algorithm : algorithms)
//       algorithm->fillTree();
//   } // end of event loop
//   std::cout << "EventLoop done\n";
// }


// // ===== Run3 class ===== //
// void EventLoopRun3::initialize()
// {
//   std::cout << "\n\n\n===== Initialize the event loop =====\n";
  
//   m_chain = new TChain(treeName);
//   for (const auto &inputFile : inputFiles) {
//     m_chain->Add(inputFile);
//     std::cout << "Added files: " << inputFile << "\n";
//   }

//   m_data = new DataRun3(m_chain, isMC);

//   for (const auto &algorithm : algorithms)
//     algorithm->initialize(m_data);

//   setNumEvent();
// }

// void EventLoopRun3::execute()
// {
//   if (!m_chain) {
//     throw std::runtime_error("Calling execute() without initializiation.\n");
//   }

//   std::cout << "Tot evet: " << nEvent << "\n";
//   for (Long64_t ir = 0; ir < nEvent; ++ir)
//   {
//     if (ir % 100000 == 0 && ir != 0) {
//       std::cout << "Event " << ir << "\n";
//     }
//     m_chain->GetEntry(ir);

//     for (Long64_t irqq = 0; irqq < m_data->Reco_QQ_size; ++irqq) {
//       for (const auto &algorithm : algorithms)
//         algorithm->execute(irqq);
//     } // reco loop end

//     for (const auto &algorithm : algorithms)
//       algorithm->fillTree();
//   } // end of event loop
//   std::cout << "EventLoop done\n";
// }