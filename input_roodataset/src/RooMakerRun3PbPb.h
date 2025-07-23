#pragma once
#include <string>

class TFile;
class TChain;
class RooDataSet;
class RooRealVar;
class RooArgSet;

class RooMakerRun3PbPb
{
public:
  RooMakerRun3PbPb(int cLow, int cHigh, float massLow, float massHigh,
                   bool isMC, bool isAccW, bool isEffW,
                   bool isTnP, bool isPtW, int hiHFBinEdge, int mcType,
                   std::string fInputPath, std::string inputChainName, std::string userTag);
  ~RooMakerRun3PbPb();

  void init();
  void run();

  // locate in public region for test
  TChain *inChain = nullptr; // input chain

private:
  void setLabels();
  void openInputFile();
  void setupTree();
  void setupRooDataSet();
  void setOutfile();
  void loadEffHist();
  void loadAccHist();
  void defineRooVariables();
  void getEffWeight();
  void getAccWeight();
  void loopAndFillData();
  void saveResults();

  int cLow, cHigh;
  float massLow, massHigh;
  bool isMC, isAccW, isEffW, isTnP, isPtW;
  int hiHFBinEdge, mcType;

  // labels
  std::string sampleTag, mcTag, centHfTag;

  // paths
  std::string fInputPath, fEffPath, fAccPath;
  std::string inputChainName;
  std::string userTag;

  // file
  TFile *fEff = nullptr;
  TFile *fAcc = nullptr;
  TFile *outfile = nullptr;
  
  RooDataSet *dataSet = nullptr;

  // setupTree variables
  int event, cBin, nDimu;
  static const int nMaxDimu = 1000;
  int recoQQsign[nMaxDimu];
  float mass[nMaxDimu], pt[nMaxDimu], y[nMaxDimu];
  float pt1[nMaxDimu], pt2[nMaxDimu], eta1[nMaxDimu], eta2[nMaxDimu];
  float ctau3D[nMaxDimu], ctau3DErr[nMaxDimu];
  float cosHX[nMaxDimu], phiHX[nMaxDimu], cosCS[nMaxDimu], phiCS[nMaxDimu];
  float cosEP[nMaxDimu], phiEP[nMaxDimu];
  // float ctau3D2S[nMaxDimu], ctau3DErr2S[nMaxDimu], ctau3DRes2S[nMaxDimu];
  float vz;
  double weight;

  // RooFit variables
  // if you add new variables, please delete meories at destructor()
  RooRealVar *massVar, *ptVar, *yVar;
  RooRealVar *pt1Var, *pt2Var, *eta1Var, *eta2Var;
  RooRealVar *cBinVar, *evtWeight, *recoQQ;
  RooRealVar *ctau3DVar, *ctau3DErrVar, *ctau3DResVar;
  RooRealVar *ctau3D2SVar, *ctau3DErr2SVar, *ctau3DRes2SVar;
  RooRealVar *NumDimu;
  RooArgSet *argSet;
  
  // angles
  RooRealVar *cosHXVar, *phiHXVar, *cosCSVar, *phiCSVar, *cosEPVar, *phiEPVar;
};