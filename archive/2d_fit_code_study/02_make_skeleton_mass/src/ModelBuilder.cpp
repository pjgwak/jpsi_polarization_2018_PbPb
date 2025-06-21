#include "ModelBuilder.h"
#include "RooAbsPdf.h"
#include <iostream>
#include <stdexcept>


ModelBuilder::ModelBuilder(RooWorkspace &ws) : m_ws(ws) {};

void ModelBuilder::setParameters(const std::map<std::string, std::vector<double>> &params) {
  m_params = params;
}

std::string ModelBuilder::getParamString(const std::string &paramName, const std::string &defaultValue) {
  if (m_params.count(paramName)) {
    // there is paramName in the config

    const auto &values = m_params.at(paramName); // bring values of a param

    if (values.size() == 3) {
      // [initial, lower, upper]
      return paramName + "[" + std::to_string(values[0]) + ", " +
             std::to_string(values[1]) + ", " + std::to_string(values[2]) + "]";
    } else if (values.size() == 1) {
      // fixed variable
      return paramName + "[" + std::to_string(values[0]) + "]";
    }
  }

  // no paramName in the config - make variable with default values
  return paramName + defaultValue;
}

void ModelBuilder::buildMassSig(std::string pdfType) {
  std::string obsName = "mass";
  std::string instanceName = pdfType + "_sig_model";
  std::string factoryCmd;

  if (pdfType == "gauss") {
    factoryCmd = "RooGaussian::" + instanceName + "(" +
                 obsName + ", " +
                 getParamString("mean_mass", "[3.096, 3.08, 3.11]") + ", " +
                 getParamString("sigma_mass", "[0.03, 0.01, 0.1]") + ")";
  } else if (pdfType == "cb") {
    factoryCmd = "RooCBShape::" + instanceName + "(" +
                 obsName + ", " +
                 getParamString("mean_mass", "[3.096, 3.08, 3.11]") + ", " +
                 getParamString("sigma_mass", "[0.03, 0.01, 0.1]") + ", " +
                 getParamString("alpha_mass", "[1.0, 0.1, 5.0]") + ", " +
                 getParamString("n_mass", "[2.0, 0.5, 10.0]") + ")";
  } else if (pdfType=="gaussCb") {
    m_ws.factory(("RooGaussian::" + instanceName + "_g(" +
                  obsName + ", " +
                  getParamString("mean_mass", "[3.096, 3.08, 3.11]") + ", " +
                  getParamString("sigma_g_mass", "[0.02, 0.001, 1])") + ")").c_str());
    m_ws.factory(("RooCBShape::" + instanceName + "_cb(" +
                  obsName + ", mean_mass, " + // same mean
                  getParamString("sigma_cb_mass", "[0.05, 0.001, 1]") + ", " +
                  getParamString("alpha_cb_mass", "[1.0, 0.1, 5.0]") + ", " +
                  getParamString("n_cb_mass", "[2.0, 0.5, 10.0]") + ")").c_str());
    factoryCmd = "RooAddPdf::" + instanceName + "({" + instanceName + "_g, " + instanceName + "_cb}, {" +
                 getParamString("frac_mass", "[0.5, 0, 1]") + "})";
  } else if (pdfType=="doubleCb") {
    // parameter first
    m_ws.factory(("RooFormulaVar::sigma2_mass_formula(\"@0*@1\", {" +
                  getParamString("sigma1_mass", "[0.05, 0.02, 0.1]") + ", " +
                  getParamString("x_mass", "[1.5, 0.1, 5.0]") + "})").c_str());
    m_ws.factory(("RooCBShape::" + instanceName + "_cb1(" +
                  obsName + ", " +
                  getParamString("mean_mass", "[3.096, 3.08, 3.11]") + ", " +
                  "sigma1_mass, " + 
                  getParamString("alpha1_mass", "[1.0, 0.1, 5.0]") + ", " +
                  getParamString("n1_mass", "[1.0, 0.5, 5.0]") + ")").c_str());
    m_ws.factory(("RooCBShape::" + instanceName + "_cb2(" +
                  obsName + ", mean_mass, " +
                  "sigma2_mass_formula, " + 
                  "alpha1_mass, " + 
                  "n1_mass" + ")").c_str());
    factoryCmd = "RooAddPdf::" + instanceName + "({" + instanceName + "_cb1, " + instanceName + "_cb2}, {" +
                 getParamString("frac_mass", "[0.5, 0, 1]") + "})";
  }
  else {
    std::cerr << "Error: We don't have mass sig model: " << pdfType << "\n";
  }

  m_ws.factory(factoryCmd.c_str());
  m_ws.pdf(instanceName.c_str())->SetName("mass_sig");
}

void ModelBuilder::buildMassBkg(std::string pdfType) {
  std::string obsName = "mass";
  std::string instanceName = pdfType + "_bkg_model";
  std::string factoryCmd;

  if (pdfType == "expo") {
    factoryCmd = "RooExponential::" + instanceName +
                 "(" + obsName + ", " + getParamString("s0_mass", "[-0.1, -5, 0]") + ")";
  }
  else if (pdfType == "cheby1") {
    factoryCmd = "RooChebychev::" + instanceName +
                 "(" + obsName + ", {" + getParamString("s1_mass", "[0.1, -1, 1]") + "})";
  }
  else if (pdfType == "cheby2") {
    factoryCmd = "RooChebychev::" + instanceName +
                 "(" + obsName + ", {" +
                 getParamString("s1_mass", "[0.1, -1, 1]") + ", " +
                 getParamString("s2_mass", "[-0.1, -1, 1]") + "})";
  }
  else if (pdfType == "cheby3") {
    factoryCmd = "RooChebychev::" + instanceName +
                 "(" + obsName + ", {" +
                 getParamString("s1_mass", "[0.1, -1, 1]") + ", " +
                 getParamString("s2_mass", "[-0.1, -1, 1]") + ", " +
                 getParamString("s3_mass", "[0.05, -1, 1]") + "})";
  }
  else {
    std::cerr << "Error: We don't have mass bkg model: " << pdfType << "\n";
  }

  m_ws.factory(factoryCmd.c_str());
  m_ws.pdf(instanceName.c_str())->SetName("mass_bkg");
}

void ModelBuilder::buildMassModel(bool isMC) {
  std::string factoryCmd;
  if (isMC) {
    factoryCmd = "RooAddPdf::massModel({mass_sig}, {" +
                 getParamString("nSig_MC", "[100000, 1, 8000000]") + "})";
  } else {
    factoryCmd = "RooAddPdf::massModel({mass_sig, mass_bkg}, {" +
                 getParamString("nSig", "[1000, 1, 50000]") + ", " +
                 getParamString("nBkg", "[50000, 1, 200000]") + "})";
  }
  
  m_ws.factory(factoryCmd);
}

// void ModelBuilder::buildTimeSig()
// {
//   // lifetime (ctau) functions
//   // 일단 mass 함수부터 만들기.
// }

// void ModelBuilder::buildTimeBkg()
// {
//   //
// }

// void ModelBuilder::combineModels() {
//   // 뭐가 됐든 sig, bkg 가져와서 더하기.
// }