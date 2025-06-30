#pragma once

#include <string>
#include <map>
#include <memory>

#include "TString.h"

class RooWorkspace;
class RooFitResult;
class RooAbsArg; // RooFit container
class RooFitResult;

class MassFit
{
public:
  MassFit(
      float ptLow = 3, float ptHigh = 6.5,
      float yLow = 1.6, float yHigh = 2.4,
      int cLow = 60, int cHigh = 180,
      float cosLow = 0, float cosHigh = 0.1,
      int PR = 0,
      int PRw = 1, bool fEffW = false, bool fAccW = false, bool isPtW = false, bool isTnP = false);
  ~MassFit();

  void run();
  void allInOne();

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
  TString kineLabel;

  // // RooFit instances
  RooWorkspace *ws = nullptr;
  RooFitResult *fitMass = nullptr;

  std::map<std::string, std::unique_ptr<RooAbsArg>> components;
};