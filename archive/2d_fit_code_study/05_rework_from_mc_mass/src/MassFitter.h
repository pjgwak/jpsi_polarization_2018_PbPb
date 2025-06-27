#pragma once

#include <string>

// Forward declarations
class RooRealVar;
class RooAbsPdf;
class RooDataSet;
class RooFitResult;
class RooWorkspace;
class TFile;

class MassFitter
{
public:
  MassFitter(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, float cosLow, float cosHigh);
  ~MassFitter();

  void run();

private:
  void setupInputOutput();
  void prepareDataset();
  void loadMcFitResults();
  void buildModel();
  void performFit();
  void plotResults();
  void saveResults();

  // kinematics & settings
  float ptLow_, ptHigh_;
  float yLow_, yHigh_;
  int cLow_, cHigh_;
  float cosLow_, cosHigh_;
  const double massLow_ = 2.6, massHigh_ = 3.5;
  const int nMassBin_ = 70; // 34 for MC

  // I/O Paths and Files
  std::string dataInputPath_ = "../../../files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root";
  std::string mcFitInputPath_;
  std::string outputPath_;
  TFile *dataFile_ = nullptr;
  TFile *mcFitFile_ = nullptr;

  // RooFit objects
  RooWorkspace *ws_mass_ = nullptr;
  RooRealVar *mass_ = nullptr;
  RooDataSet *ds_mass_ = nullptr;
  // RooAbsPdf *mass_pdf_ = nullptr;
  // RooAbsPdf *mass_sig_ = nullptr;
  // RooAbsPdf *mass_bkg_ = nullptr;
  // RooAbsPdf *G1Sig = nullptr;
  // RooAbsPdf *CB1Sig = nullptr;
  RooFitResult *fit_mass_ = nullptr;

  // model parameters
  // signal parameters - bring from MC fit
  RooRealVar *mean_mass_ = nullptr;
  RooRealVar *sigma1_mass_ = nullptr;
  RooRealVar *sigma2_mass_ = nullptr;
  RooRealVar *alpha1_mass_ = nullptr;
  RooRealVar *power1_mass_ = nullptr;
  RooRealVar *fracG1_mass_ = nullptr;

  // bkg parameters
  RooRealVar *sl1_mass_ = nullptr;
  RooRealVar *sl2_mass_ = nullptr;

  // yields
  RooRealVar *nSig_mass_ = nullptr;
  RooRealVar *nBkg_mass_ = nullptr;
};