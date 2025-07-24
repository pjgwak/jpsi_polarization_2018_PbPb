#pragma once
#include <string>

class TFile;
class RooWorkspace;

class CtauErrFit
{
public:
  CtauErrFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
             float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);
  ~CtauErrFit();

  void run();

  std::string inputFilePath;
  std::string DATE;
  double ctauErrMin = -1; // -1 means auto-range
  double ctauErrMax = -1;
  

private:
  void setLabels();
  void openInputFiles();
  void setupWorkspaceAndData();
  void performSPlot();
  void setVariableRanges();
  void buildPdfAndData();
  void makePlot();
  void saveResults();

  float ptLow, ptHigh, yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR, PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  std::string kineLabel;
  std::string fname;

  int nBins = 100;

  TFile *fInputData = nullptr;
  TFile *fMass = nullptr;
  RooWorkspace *ws = nullptr;
};