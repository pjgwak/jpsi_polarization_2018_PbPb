#pragma once
#include <string>

class TFile;
class RooWorkspace;
class RooFitResult;

class MassFit
{
public:
  MassFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
          float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);

  ~MassFit();

  void init();
  void run();

  std::string inputFilePath;
  std::string inputMcPath;
  std::string DATE;
  std::string pdfTypeSig;
  std::string pdfTypeBkg;
  bool isWeighted = false;

  void initVar(const std::string &varName, double init, double low, double high);

private:
  void setLabels();
  void openInputFile();
  void setupWorkspaceAndData();
  void setVariableRanges();
  void defineModel();
  void performFit();
  void makePlot();
  void saveResults();

  // pdfs
  void buildDoubleCB();
  void buildCBG();
  void buildExpo();
  void buildCheby(int order); // cheby 1 ~ 6


  float ptLow, ptHigh, yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR, PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  std::string kineLabel;
  std::string fname;

  TFile *fInputData = nullptr;
  TFile *fMcParams = nullptr;
  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;
};