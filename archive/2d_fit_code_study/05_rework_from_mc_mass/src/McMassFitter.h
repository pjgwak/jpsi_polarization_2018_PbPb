#pragma once

#include <string>

// forward declarations
class RooRealVar;
class RooAbsPdf;
class RooDataSet;
class RooFitResult;
class TFile;

class McMassFitter {
public:
  McMassFitter(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, float cosLow, float cosHigh);
  ~McMassFitter();

  void run();

private:
  void setupInputOutput();
  void prepareDataset();
  void buildModel();
  void performFit();
  void plotResults();
  void saveResults();

  const double massLow_ = 2.6,
               massHigh_ = 3.5;
  const int nMassBin_ = 34;

  // kinematics
  float ptLow_, ptHigh_;
  float yLow_, yHigh_;
  int cLow_, cHigh_;
  float cosLow_, cosHigh_;
  float massFitMin_ = 2.6;
  float massFitMax_ = 3.25;

  std::string inputPath_ = "../../../files_roodata/RooDataSet_miniAOD_isMC1_PR_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root";
  std::string outpath_;
  TFile *f_mc_ = nullptr;

  // RooFit dataset and meta variables
  RooDataSet *ds_tmp_ = nullptr;
  RooDataSet *ds_mc_ = nullptr;
  RooAbsPdf *mass_pdf_mc_ = nullptr;
  RooFitResult *fit_result_mc_ = nullptr;
  RooRealVar *mass_mc_ = nullptr;

  // Model Parameters
  RooRealVar *mean_mc_ = nullptr;
  RooRealVar *sigma1_mc_ = nullptr;
  RooRealVar *sigma2_mc_ = nullptr;
  RooRealVar *alpha1_mc_ = nullptr;
  RooRealVar *power1_mc_ = nullptr;
  RooRealVar *fracG1_mc_ = nullptr;
  RooRealVar *nSig_mc_ = nullptr;
};