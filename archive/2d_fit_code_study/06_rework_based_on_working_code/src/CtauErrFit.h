#pragma once

#include <string>
#include <map>
#include <memory>

#include "TString.h"

class RooWorkspace;
class RooAbsArg;
namespace RooStats {
  class SPlot;
}

class CtauErrFit {
public:
  CtauErrFit(float ptLow, float ptHigh,
             float yLow, float yHigh,
             int cLow, int cHigh,
             float cosLow, float cosHigh,
             int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);
  ~CtauErrFit();

  void run();

private:
  void setLabels();
  void prepareDataset();
  void performSPlot();
  void buildCtauErrModel();
  void drawResult();
  void saveResult();

  float ptLow, ptHigh;
  float yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR;
  int PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  TString DATE;
  TString fname;
  TString kineLabel;

  double ctauErrMin = -10;
  double ctauErrMax = 10;
  int newBins = 100;

  RooWorkspace *ws = nullptr;
  RooStats::SPlot *sData = nullptr;
};