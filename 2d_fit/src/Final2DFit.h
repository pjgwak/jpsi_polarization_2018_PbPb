#pragma once

#include <string>

class TFile;
class RooWorkspace;
class TString;
class RooFitResult;

class Final2DFit
{
public:
  Final2DFit(float ptLow, float ptHigh,
             float yLow, float yHigh,
             int cLow, int cHigh,
             float cosLow, float cosHigh,
             int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);
  ~Final2DFit();

  void run();
  void makePlots();

private:
  void setLabels();
  void createKinematicCut();
  void openInputFiles();
  void setupWorksapceAndData();
  void setVariableRanges();
  void defineAndFitMassModel();
  void defineCtauModel();
  void define2DFitModel();
  void perform2DFit();
  void plotMass();
  void plotCtau();
  void saveResults();

  RooWorkspace *ws = nullptr;
  TFile *f1 = nullptr; // original RooDataSet - Though it's bad nameing, following name of original codes.
  TFile *fMass = nullptr;
  TFile *fCErr = nullptr;
  TFile *fCRes = nullptr;
  TFile *fCBkg = nullptr;
  TFile *fCTrue = nullptr;

  std::string kineCut;
  std::string kineLabel;

  float ptLow, ptHigh;
  float yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR;
  int PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  RooFitResult *fitMass = nullptr;
  RooFitResult *fit2DFit = nullptr;

  std::string DATE;
  std::string fname;
};