#pragma once

#include <vector>
#include "TString.h"
#include "TChain.h"
#include "Data.h"
#include "Algorithm.h"

class EventLoop {
	public:  
    /**
     * @brief construct a new event loop object
     */
		EventLoop();
		// virtual ~EventLoop() = default;

		/**
		 * @brief initialize the event loop
		 */
		void initialize(); // means pure virtual function. Childern must have it.

		/**
		 * @brief execute the event loop
		 */
		void execute();
		void executeGen();

		/**
		 * @brief list of input ROOT file names
		 */
		std::vector<TString> inputFiles;

		/**
		 * @brief name of the TTree instance. Must be same in all files
		 */
		TString treeName;

		/**
		 * @brief sample type flags
		 */
		bool isMC = false;
		bool isGenOnly = false;
		bool isEP = false;
		bool useRun2 = false;

		/**
		 * @brief lost of algorithms to be executed in the loop
		 */
		std::vector<Algorithm*> algorithms;

		/**
		 * @brief number of events to run a event loop
		 * -1: means using all event
		 */
		Long64_t nEvent = -1;

		/**
		 * @brief set a number of events to use
		 */
		void setNumEvent();

		/**
		 * @brief instances of the Data classes
		 */
		Data *m_data = nullptr;
		DataRun2 *data_run2 = nullptr;
		DataRun3 *data_run3 = nullptr;

	protected:
		/**
		 * @brief instance of the TChain class used to read the data
		 */
		TChain *m_chain = nullptr;
};


// // ====== Childern classes ===== //
// class EventLoopRun2 {
// 	public:
// 		EventLoopRun2();
// 		~EventLoopRun2() override = default;

// 		/**
// 		 * @brief initialize the event loop
// 		 */
// 		void initialize()override;
// 		void execute() override;

// 		/**
// 		 * @brief instance of the Data class
// 		 */
// 		DataRun2 *m_data = nullptr;
// };

// class EventLoopRun3 {
// 	public:
// 		EventLoopRun3();
// 		~EventLoopRun3() override = default;

// 		void initialize() override;
// 		void execute() override;

// 		DataRun3 *m_data = nullptr;
// };