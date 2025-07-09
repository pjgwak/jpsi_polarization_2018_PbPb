#pragma once
#include <string>

class TFile;
class RooWorkspace;
class RooFitResult;

class CtauBkgFit
{
public:
  CtauBkgFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
             float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);
  ~CtauBkgFit();
  void run();

private:
  void setLabels();
  void openInputFile();
  void createKinematicCut();
  void setupWorkspaceAndData();
  void setVariableRanges();
  void defineModel();
  void performFit();
  void drawPullPlot();
  void drawRatioPlot();
  void saveResults();

  float ptLow, ptHigh, yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR, PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  std::string kineLabel;
  std::string fname;
  std::string DATE;

  TFile *fMass = nullptr;
  TFile *fCErr = nullptr;
  TFile *fCRes = nullptr;
  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;

  int nGauss = 2;
  double ctauMin = -100.; // -100: use automatic range
  double ctauMax = 100; // 100: use automatic range
};
