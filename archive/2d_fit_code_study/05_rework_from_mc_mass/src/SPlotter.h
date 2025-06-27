#pragma once

#include <string>

class RooRealVar;
class RooAbsPdf;
class RooDataSet;
class RooFitResult;
class RooWorkspace;
class TFile;
class TH1D;
class RooPlot;

namespace RooStats
{
  class SPlot;
}

class SPlotter
{
public:
  SPlotter(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, float cosLow, float cosHigh);
  ~SPlotter();

  void run();

private:
  void setupInputOutput();
  void prepareDataset();
  void performSPlot();
  void buildWeighedObjects();
  void plotResults();
  void saveResults();

  // kinematics & settings
  float ptLow_, ptHigh_;
  float yLow_, yHigh_;
  int cLow_, cHigh_;
  float cosLow_, cosHigh_;

  // I/O Paths and Files
  std::string dataInputPath_ = "../../../files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root";
  std::string massFitInputPath_;
  std::string outputPath_;
  TFile *dataFile_ = nullptr;
  TFile *massFitFile_ = nullptr;

  // RooFit & RooStats Objects
  RooWorkspace *ws_err_ = nullptr;
  RooRealVar *ctau3DErr_ = nullptr;
  RooStats::SPlot *sPlotData_ = nullptr;

  double ctauErrMin_ = 0;
  double ctauErrMax_ = 0.25;
  int nBins_ = 100;
  int newBins_ = 100;
};