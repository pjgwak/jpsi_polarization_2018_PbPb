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

		/**
		 * @brief initialize the event loop
		 */
		void initialize();

		/**
		 * @brief execute the event loop
		 */
		void execute();

		/**
		 * @brief list of input ROOT file names
		 */
		std::vector<TString> inputFiles;

		/**
		 * @brief name of the TTree instance. Must be same in all files
		 */
		TString treeName;

		/**
		 * @brief is the sample MC?
		 */
		bool isMC = false;

		/**
		 * @brief lost of algorithms to be executed in the loop
		 */
		std::vector<Algorithm*> algorithms;
	
	protected:
		/**
		 * @brief instance of the TChain class used to read the data
		 */
		TChain* m_chain = nullptr;

		/**
     * @brief instance of the Data class
     */
		Data* m_data = nullptr;
};