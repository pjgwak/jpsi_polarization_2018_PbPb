#pragma once

#include <string>
#include <map>
#include <memory>

#include "TString.h"

class RooWorkspace;
class RooFitResult;
class RooAbsArg; // RooFit container
class RooFitResult;

class McMassFit
{
public:
  McMassFit(
      float ptLow = 3, float ptHigh = 6.5,
      float yLow = 1.6, float yHigh = 2.4,
      int cLow = 60, int cHigh = 180,
      float cosLow = 0, float cosHigh = 0.1,
      int PR = 0,
      int PRw = 1, bool fEffW = false, bool fAccW = false, bool isPtW = false, bool isTnP = false);
  ~McMassFit();

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

  // // IO labels
  // std::string nameTag; // user custom name tag
  // std::string PrTag; // Prompt, nonprompt

  // // paths
  TString DATE;
  TString fname;

  // // RooFit instances
  RooWorkspace *ws = nullptr;
  RooFitResult *fitMass = nullptr;

  // double massFitMin = 2.6;
  double fit_limit = 3.25;
  
  std::map<std::string, std::unique_ptr<RooAbsArg>> components;
};