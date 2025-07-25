#pragma once
#include <string>

class TFile;
class RooWorkspace;
class RooFitResult;

class CtauBkgFit
{
public:
  CtauBkgFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
             float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);
  ~CtauBkgFit();
  void init();
  void run();

  std::string DATE;

  int nGauss = 2, nExp_L = 2, nExp_C = 1, nExp_R = 2;
  double ctauMin = -100.; // -100: use automatic range
  double ctauMax = 100;   // 100: use automatic range

  void initVar(const std::string &varName, double init, double low, double high);
  void setConstVar(const std::string &varName, bool isConst, double value = 3096);

private:
  void setLabels();
  void openInputFile();
  void createKinematicCut();
  void setupWorkspaceAndData();
  void setVariableRanges();
  void defineModel();
  void performFit();
  void drawPullPlot();
  void drawRatioPlot();
  void saveResults();

  // pdfs
  void buildResModel();
  void buildLeftModel();
  void buildCenterModel();
  void buildRightModel();
  void combineDecayModels();

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
  TFile *fCRes = nullptr;
  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;
};
