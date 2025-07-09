#pragma once
#include <string>

class TFile;
class RooWorkspace;
class RooFitResult;

class McMassFit
{
public:
  McMassFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
            float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);
  ~McMassFit();

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
  double massMin = 2.6;
  double massMax = 3.22;

  // 핵심 객체 포인터 (nullptr로 초기화)
  TFile *fInputData = nullptr;
  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;
};