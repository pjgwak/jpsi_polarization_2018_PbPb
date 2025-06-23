#pragma once

#include <string>
#include <vector>
#include <map>
#include "RooWorkspace.h"
#include "TH1D.h"
#include "RooPlot.h"
#include "memory" // std::unique_ptr

// SPlot
#include "RooStats/SPlot.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"

class JpsiFitter
{
public:
  JpsiFitter(RooWorkspace &ws);
  ~JpsiFitter();

  void processTree(const std::string &filePath, const std::string &dsname, const std::string &obs, const std::string &cut);
  void print(); // show the contents of workspace

  // SPlot functions
  void loadMassResult(const std::string &filePath);
  void doSplot(const std::string &filePath, const std::string &cut);
  void makeSplotPdfs(bool isForcedMax, double ctauErrMax_forced);
  void drawSplot(const std::string &outPlotPath);

private:
  RooWorkspace &m_ws;
  // TFile *mass_result = nullptr;

  // ===== SPlot members =====
  double m_ctauErrMin;
  double m_ctauErrMax;
  double m_minRange;
  double m_maxRange;
  int m_nBins;
  int m_newBins;

  double m_Ydown;
  double m_Yup;

  double m_outTot;
  double m_outRes;

  RooPlot *err_frame1;
  TH1D *m_hTot;
};