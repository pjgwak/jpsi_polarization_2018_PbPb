#pragma once
#include <string>

class TFile;
class RooWorkspace;
class RooFitResult;

class MassFit 
{
public:
  MassFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
          float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);

  ~MassFit();
  void run();

private:
  void setLabels();
  void openInputFile();
  void setupWorkspaceAndData();
  void setVariableRanges();
  void initializeParameters(); // get MC result
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

  double mc_alpha, mc_n, mc_sigma, mc_x, mc_f;

  TFile *fInputData = nullptr;
  TFile *fMcParams = nullptr;
  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;
};