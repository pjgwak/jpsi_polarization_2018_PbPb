#pragma once
#include <string>
#include <vector>
#include <map>
#include "../headers/json.hpp"

class TFile;
class RooWorkspace;
class RooFitResult;
class RooPlot;

using json = nlohmann::json;

class MassFit
{
public:
  MassFit(const std::string &global_config_path, const std::string &local_config_path);
  ~MassFit();
  void run();

private:
  void loadConfiguration(const std::string &global_path, const std::string &local_path);
  void setupWorkspaceAndData();
  void setVariableRanges();
  void initializeParameters();
  void defineModel();
  void performFit();
  void makePlot();
  void saveResults();

  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;
  TFile *fInputData = nullptr;
  TFile *fMcParams = nullptr;

  // --- parameters from MC fit ---
  double mc_alpha;
  double mc_n;
  double mc_sigma;
  double mc_x;
  double mc_f;

  json global_config;
  json local_config;
  json massfit_config;

  std::string full_input_path;
  std::string output_root_path;
  std::string output_fig_path;
  std::string output_base_name;

  float ptLow, ptHigh, yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR, PRw;
  bool fEffW, fAccW;
  bool isPtW, isTnP;
  std::string fname;

  json plotting_style;
  int nMassBin;
  double massPlotMin, massPlotMax;
  double massFitMin, massFitMax;
};