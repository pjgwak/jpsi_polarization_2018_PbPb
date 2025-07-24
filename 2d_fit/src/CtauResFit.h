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
  void init();
  void run();

  std::string DATE;
  int nGauss = 2;
  double ctauResMin = -100; // -100: automatic range
  double ctauResMax = 0.0;  // 0: nominal use only Res < 0 region, 100: automatic range

  void initVar(const std::string &varName, double init, double low, double high);

private:
  void setLabels();
  void openInputFile();
  void setupWorkspaceAndData();
  void setVariableRanges();
  void defineModel();
  void performFit();
  void makePlot();
  void makeRatioPlot();
  void saveResults();

  // pdfs
  void buildCtauRes1Gaus();
  void buildCtauRes2Gaus();
  void buildCtauRes3Gaus();
  void buildCtauRes4Gaus();

  // parameter print
  void drawTextVar(const char *varName, const char *label, float xp, float yp, int textColor, int textSize);
  void drawTextVarInt(const char *varName, const char *label, float xp, float yp, int textColor, int textSize);

  float ptLow, ptHigh, yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR, PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  std::string kineLabel;
  std::string fname;

  TFile *fMass = nullptr;
  TFile *fCErr = nullptr;
  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;
};