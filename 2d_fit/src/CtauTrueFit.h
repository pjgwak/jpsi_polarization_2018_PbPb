#pragma once

#include <string>

class TFile;
class RooWorkspace;
class RooFitResult;

class CtauTrueFit
{
public:
  CtauTrueFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
              float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);
  ~CtauTrueFit();

  void init();
  void run();
  
  std::string inputFilePath;
  std::string DATE;
  int nExp = 2;
  double ctau3DMin = 1e-5; // should be > 1e-5
  double ctau3DMax = 10;

  void initVar(const std::string &varName, double init, double low, double high);

private:
  void setLabels();
  void openInputFile();
  void createKinematicCut();
  void setupWorksapceAndData();
  void setVariableRanges();
  void defineModel();
  void performFit();
  void makePlot();
  void saveResults();

  // pdfs
  void buildTrue1Expo();
  void buildTrue2Expo();
  void buildTrue3Expo();
  void buildTrue4Expo();

  // parameter print
  void drawTextVar(const char *varName, const char *label, float xp, float yp, int textColor, int textSize);
  void drawTextVarInt(const char *varName, const char *label, float xp, float yp, int textColor, int textSize);

  float ptLow, ptHigh, yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR, PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  std::string kineLabel;
  std::string kineCutMC;

  TFile *inputFile = nullptr;
  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;
};