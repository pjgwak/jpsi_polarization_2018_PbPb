#pragma once
#include <string>

class TFile;
class RooWorkspace;
class RooFitResult;

class CtauResFit
{
public:
  CtauResFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
             float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);
  ~CtauResFit();
  void run();

private:
  void setLabels();
  void openInputFile();
  void setupWorkspaceAndData();
  void setVariableRanges();
  void defineModel();
  void performFit();
  void makePlot();
  void saveResults();

  float ptLow, ptHigh, yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR, PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  std::string kineLabel;
  std::string fname;
  std::string DATE;
  int nGauss = 2;            
  double ctauResMin = -100; // -100: automatic range
  double ctauResMax = 0.0; // 0: nominal use only Res < 0 region, 100: automatic range

  TFile *fMass = nullptr;
  TFile *fCErr = nullptr;
  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;
};