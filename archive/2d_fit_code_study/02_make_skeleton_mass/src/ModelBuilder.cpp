#include "ModelBuilder.h"
#include <iostream>
#include <stdexcept>


ModelBuilder::ModelBuilder(RooWorkspace &ws) : m_ws(ws) {};

void ModelBuilder::buildMassSignal() {
  std::string pdfType = "doubleCb";
  // 원하는 게 뭔지 매개변수로 받아야 한다.
  // Ws 안에서 겹치면 안 되니까. 이름도 유니크하게 지어야 한다.
  // 피팅 단계 + 모델명 정도면 괜찮을 듯? -> 변수 이름 받고 여기서 피팅 단계를 prefix로 붙여준다.

  // 초기값 설정하는 부분은 따로 함수로 만들거나?
  // 기본값으로 만들고 상한, 하한, 초기값은 변경할 수 있으니까?

  // 첫 단계니까 함수 다 구현하고 선택하는 부분까지만 만들고 다음 단계로.
  std::string obsName = "mass"; // ctau 쪽에서는 선택 가능하게 -> GenTrue랑 Reco랑 같은 이름으로 하면 하드 코딩 해도 된다. 이 부분은 사용자가 바꿀 이유가 없어서.
  std::string instanceName = m_ws.GetName() + std::string("_") + obsName + "_sig"; // 테스트용!
  std::string factoryCmd;

  if (pdfType == "gauss") {
    factoryCmd = "RooGaussian::" + instanceName + "(" +
                 obsName + ", mean_" + instanceName + "[3.096, 3.086, 3.106]," +
                 "sigma_" + instanceName + "[0.03, 0.01, 0.1])";
  } else if (pdfType == "cb") {
    factoryCmd = "RooCBShape::" + instanceName + "(" +
                 obsName + ", mean_" + instanceName + "[3.096, 3.086, 3.106]," +
                 "sigma_" + instanceName + "[0.03, 0.01, 0.1], " +
                 "alpha_" + instanceName + "[1.0, 0.1, 5.0], " +
                 "n_" + instanceName + "[2.0, 0.5, 10.0])";
  } else if (pdfType=="gaussCb") {
    m_ws.factory(("RooGaussian::" + instanceName + "_g(" +
                  obsName + ", mean_" + instanceName + "[3.096, 3.086, 3.106], " +
                  "sigma_g_" + instanceName + "[0.02, 0.001, 1])").c_str());
    m_ws.factory(("RooCBShape::" + instanceName + "_cb(" +
                  obsName + ", mean_" + instanceName + ", " + // same mean
                  "sigma_cb_" + instanceName + "[0.05, 0.001, 1], " +
                  "alpha_cb_" + instanceName + "[1.0, 0.1, 5.0], " +
                  "n_cb_" + instanceName + "[2.0, 0.5, 10.0])").c_str());
    factoryCmd = "RooAddPdf::" + instanceName + "({" + instanceName + "_g, " + instanceName + "_cb}, " + "{frac_" + instanceName + "[0.5, 0, 1]})";
  } else if (pdfType=="doubleCb") {
    m_ws.factory(("RooCBShape::" + instanceName + "_cb1(" +
                  obsName + ", mean_" + instanceName + "[3.096, 3.086, 3.106], " +
                  "sigma1_" + instanceName + "[0.05, 0.02, 0.1], " +
                  "alpha1_" + instanceName + "[1.0, 0.1, 5.0], " +
                  "n1_" + instanceName + "[1.0, 0.5, 5.0])").c_str());
    m_ws.factory(("RooCBShape::" + instanceName + "_cb2(" +
                  obsName + ", mean_" + instanceName + ", " + // same mean
                  "sigma2_" + instanceName + "[0.05, 0.02, 0.1], " +
                  "alpha2_" + instanceName + "[1.0, 0.1, 5.0], " +
                  "n2_" + instanceName + "[2.0, 0.5, 5.0])").c_str());
    factoryCmd = "RooAddPdf::" + instanceName + "({" + instanceName + "_cb1, " + instanceName + "_cb2}, " + "{frac_" + instanceName + "[0.5, 0, 1]})";
  }
  else {
    std::cerr << "Error: We don't have mass sig model: " << pdfType << "\n";
  }

  m_ws.factory(factoryCmd.c_str());
  // gaussCb
  // doubleCb
}

void ModelBuilder::buildMassBkg() {
  std::string pdfType = "cheby1";
  std::string obsName = "mass";
  std::string instanceName = m_ws.GetName() + std::string("_") + obsName + "_bkg"; // 테스트용!
  std::string factoryCmd;

  // exp
  if (pdfType == "expo") {
    factoryCmd = "RooExponential::" + instanceName +
                 "(" + obsName + ", {c1_" + instanceName + "[0, -1, 1]})";
  } else if (pdfType == "cheby1") {
    factoryCmd = "RooChebychev::" + instanceName +
                 "(" + obsName + ", {c1_" + instanceName + "[0, -1, 1]})";
  } else if (pdfType == "cheby2") {
    factoryCmd = "RooChebychev::" + instanceName +
                 "(" + obsName + ", {c1_" + instanceName + "[0, -1, 1], " +
                 "c2_" + instanceName + "[0, -1, 1]})";
  } else if (pdfType == "cheby3") {
    factoryCmd = "RooChebychev::" + instanceName +
                 "(" + obsName + ", {c1_" + instanceName + "[0, -1, 1], " +
                 "c2_" + instanceName + "[0, -1, 1], " +
                 "c3_" + instanceName + "[0, -1, 1]})";
  }
  else {
    std::cerr << "Error: We don't have mass bkg model: " << pdfType << "\n";
  }

  // cheby4 - 필요하면?
  // cheby5 

  m_ws.factory(factoryCmd.c_str());
}

void ModelBuilder::buildMassModel() {

  std::string factoryCmd = "RooAddPdf::massModel({ws_mass_sig}, {nSig[100000, 1, 8000000]})";
  // std::string factoryCmd = "RooAddPdf::massModel({ws_mass_sig, ws_mass_bkg}, {nSig[1000, 1, 50000], nBkg[50000,1, 200000]})";
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