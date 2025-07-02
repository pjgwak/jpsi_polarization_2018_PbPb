#pragma once

#include <string>
#include <map>
#include <memory>

#include "TString.h"

class RooWorkspace;
class RooFitResult;
class RooAbsArg;

class CtauResFit
{
public:
  CtauResFit(
      float ptLow = 3, float ptHigh = 6.5,
      float yLow = 1.6, float yHigh = 2.4,
      int cLow = 60, int cHigh = 180,
      float cosLow = 0, float cosHigh = 0.1,
      int PR = 0,
      int PRw = 1, bool fEffW = false, bool fAccW = false, bool isPtW = false, bool isTnP = false);
  ~CtauResFit();

  void run();

private:
  void setLabels();
  void reduceDataset();
  void buildModel();
  void fitModel();
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

  RooWorkspace* ws = nullptr;
  RooFitResult* fitCtauRes = nullptr;

  int nGauss = 2;
  double ctauResMin, ctauResMax;

  std::map<std::string, std::shared_ptr<RooAbsArg>> components;
};