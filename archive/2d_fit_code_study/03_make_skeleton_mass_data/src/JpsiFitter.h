#pragma once

#include <string>
#include <vector>
#include <map>
#include "RooWorkspace.h"
#include "memory" // std::unique_ptr

// SPlot
#include "RooStats/SPlot.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"

class JpsiFitter {
public:
  JpsiFitter(RooWorkspace &ws);
  ~JpsiFitter();

  void processTree(const std::string &filePath, const std::string &dsname, const std::string &obs, const std::string &cut);
  void print(); // show the contents of workspace

  // SPlot functions
  void loadMassResult();
  void doSplot();
  void makeSplotPdfs();
  void drawSplot();

private:
  RooWorkspace &m_ws;
  // TFile *mass_result = nullptr;
};