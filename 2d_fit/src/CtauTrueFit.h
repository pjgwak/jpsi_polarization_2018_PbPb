#pragma once

#include <string>

class TFile;
class RooWorkspace;
class RooFitResult;

class CtauTrueFit
{
public:
  CtauTrueFit(float ptLow, float ptHigh,
             float yLow, float yHigh,
             int cLow, int cHigh,
             float cosLow, float cosHigh,
             int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);
  ~CtauTrueFit();

  void run();

private:
  void setLabels();
  void openInputFile();
  void createKinematicCut();
  void setupWorksapceAndData();
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
  std::string kineCutMC;
  std::string DATE;

  TFile *inputFile = nullptr;
  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;

  double ctau3DMin = 1e-5; // should be > 1e-5
  double ctau3DMax = 10;
  int nExp = 2;
};