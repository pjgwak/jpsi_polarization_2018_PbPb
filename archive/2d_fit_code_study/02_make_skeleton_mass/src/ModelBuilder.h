#pragma once

#include "RooWorkspace.h"
#include <string>
#include <vector>
#include <map>

class ModelBuilder {
public:
  ModelBuilder(RooWorkspace &ws);

  // get params from python wrapper
  void setParameters(const std::map<std::string, std::vector<double>> &params);

  // build pdfs
  void buildMassSig(std::string pdfType);
  void buildMassBkg(std::string pdfType);
  void buildMassModel(bool isMC);
  // void combineModels();

private:
  RooWorkspace &m_ws;
  std::map<std::string, std::vector<double>> m_params;

  std::string getParamString(const std::string &paramName, const std::string &defaultValue);
};