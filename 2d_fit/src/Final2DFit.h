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

  void init();
  void run();
  void makePlots();

  std::string DATE;
  std::string inputFilePath;
  std::string pdfTypeMassSig;
  std::string pdfTypeMassBkg;

  double nGauss, nExpTrue, nExpBkg_L, nExpBkg_R;

  void initVar(const std::string &varName, double init, double low, double high);
  void setConstVar(const std::string &varName, bool isConst, double value = 3096);

private:
  void setLabels();
  void openInputFiles();
  void createKinematicCut();
  void setupWorksapceAndData();
  void setVariableRanges();
  void buildMassModel();
  void buildCtauModel();
  void build2DFitModel();
  void define2DFitModel();
  void defineAndFitMassModel();
  void defineCtauModel();
  void performMassFit();
  void perform2DFit();
  void plotMass();
  void plotCtau();
  void saveResults();

  // mass pdfs
  void buildDoubleCB();
  void buildCBG();
  void buildExpo();
  void buildCheby(); // cheby 1 ~ 6

  // ctau pdfs
  void buildCtauResModel();
  void buildCtauBkgLeftModel();
  void buildCtauBkgRightModel();
  void combineCtauBkgModels();
  void buildCtauTrueModel();

  // parameter print
  void drawTextVar(const char *varName, const char *label, float xp, float yp, int textColor, int textSize);
  void drawTextVarInt(const char *varName, const char *label, float xp, float yp, int textColor, int textSize);

  RooWorkspace *ws2dFit = nullptr;
  RooWorkspace *wsInput = nullptr;
  RooWorkspace *ws = nullptr; //test
  TFile *fInputData = nullptr;
  TFile *f1 = nullptr; //test
  TFile *fMass = nullptr;
  TFile *fCErr = nullptr;
  TFile *fCRes = nullptr;
  TFile *fCBkg = nullptr;
  TFile *fCTrue = nullptr;
  TFile *fMcParams = nullptr;

  std::string kineCut;
  std::string kineLabel;

  std::string fname;

  float ptLow, ptHigh;
  float yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR;
  int PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  RooFitResult *fitMass = nullptr;
  RooFitResult *fit2DFit = nullptr;
};