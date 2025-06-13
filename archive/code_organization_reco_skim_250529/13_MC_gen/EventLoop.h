#pragma once

#include <vector>
#include "TString.h"
#include "TChain.h"
#include "Data.h"
#include "DataRun2.h"
#include "DataRun3.h"
#include "Algorithm.h"

class EventLoop {
	public:  
		EventLoop();
		void initialize();
		void execute();
		void executeGen();

		// store input file list and tree name
		std::vector<TString> inputFiles;
		TString treeName;

		// sample type flags
		bool isMC = false;
		bool isGenOnly = false;
		bool isEP = false;
		bool useRun2 = false;

		/**
		 * @brief list of algorithms to be executed in the loop
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

		// Data classes. It can cover both Run2 and Run3
		Data *m_data = nullptr;

	protected:
		// instance of the TChain class used to read the data
		TChain *m_chain = nullptr;
};