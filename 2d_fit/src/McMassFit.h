#pragma once
#include <string>
#include <vector>
#include <map>
#include"../headers/json.hpp"

class TFile;
class RooWorkspace;
class RooFitResult;
class RooPlot;

using json = nlohmann::json;

class McMassFit
{
public:
  McMassFit(const std::string &global_config_path, const std::string &local_config_path);
  ~McMassFit();

  void run();

private:
  void loadConfiguration(const std::string &global_path, const std::string &local_path);
  
  void setupWorkspaceAndData();
  void setVariableRanges();
  void defineModel();
  void performFit();
  void makePlot();
  void drawPlotInfo(RooPlot *frame);
  void saveResults();

  TFile *fInputData = nullptr;
  RooWorkspace *ws = nullptr;
  RooFitResult *fitResult = nullptr;

  json global_config;
  json local_config;
  std::string full_input_path;
  std::string output_root_path;
  std::string output_fig_path;
  std::string output_base_name;
  // int nCPU;
  
  // kinematics and flags
  float ptLow, ptHigh, yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR, PRw;
  bool fEffW, fAccW, isPtW, isTnP;
  std::string fname; // PR, NP

  // model info
  json yields;
  std::string signal_model_type;
  json signal_model_params;
  std::vector<std::string> signal_component_names;
  std::vector<std::string> bkg_component_names;

  // plotting info
  std::vector<std::string> params_to_plot;
  json plotting_style;
  int nMassBin;
  double massPlotMin, massPlotMax;
  double massFitMin, massFitMax;
};