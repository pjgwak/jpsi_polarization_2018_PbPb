#pragma once
#include <string>

class TFile;
class RooWorkspace;
class RooFitResult;

class McMassFit
{
public:
  McMassFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
            float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);
  ~McMassFit();

  void init();
  void run();

  std::string inputFilePath;
  std::string DATE;
  std::string pdfType;
  bool isWeighted = false;

  double massMin = 2.6;
  double massMax = 3.22;

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
  void drawTextVar(const char *varName, const char *label, float xp, float yp, int textColor, int textSize);
  void drawTextVarInt(const char *varName, const char *label, float xp, float yp, int textColor, int textSize);
  void saveResults();

  // pdfs
  void buildDoubleCB();
  void buildCBG();

  float ptLow, ptHigh, yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR, PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  // lablels and paths
  std::string kineLabel;
  std::string fname;

  TFile *fInputData = nullptr;
  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;
};