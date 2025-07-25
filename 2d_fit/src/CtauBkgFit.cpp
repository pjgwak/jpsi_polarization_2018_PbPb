#include "CtauBkgFit.h"
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH1D.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "../headers/polarizationUtilities.h"

using namespace RooFit;

CtauBkgFit::CtauBkgFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
                       float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP)
    : ptLow(ptLow), ptHigh(ptHigh), yLow(yLow), yHigh(yHigh), cLow(cLow), cHigh(cHigh),
      cosLow(cosLow), cosHigh(cosHigh), PR(PR), PRw(PRw),
      fEffW(fEffW), fAccW(fAccW), isPtW(isPtW), isTnP(isTnP)
{
  gStyle->SetEndErrorSize(0);

  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);

  // to avoid ternimal breaking due to too long Tracing Error message
  // It's okay as long as minimizer find good fitting point
  // Always check fit results.
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Tracing);
  RooMsgService::instance().getStream(1).removeTopic(Tracing);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
}

CtauBkgFit::~CtauBkgFit()
{
  RooMsgService::instance().getStream(1).addTopic(RooFit::ObjectHandling);

  RooMsgService::instance().getStream(0).addTopic(Caching);
  RooMsgService::instance().getStream(1).addTopic(Caching);
  RooMsgService::instance().getStream(0).addTopic(Plotting);
  RooMsgService::instance().getStream(1).addTopic(Plotting);
  RooMsgService::instance().getStream(0).addTopic(Tracing);
  RooMsgService::instance().getStream(1).addTopic(Tracing);
  RooMsgService::instance().getStream(0).addTopic(Integration);
  RooMsgService::instance().getStream(1).addTopic(Integration);
  RooMsgService::instance().getStream(1).addTopic(NumIntegration);

  delete fMass;
  delete fCErr;
  delete fCRes;
  delete ws;
  delete fitResult;
}

void CtauBkgFit::init()
{
  setLabels();
  openInputFile();
  setupWorkspaceAndData();
  setVariableRanges();
  defineModel();
}

void CtauBkgFit::run()
{
  performFit();
  drawPullPlot();
  drawRatioPlot();
  saveResults();
}

void CtauBkgFit::setLabels()
{
  std::cout << "===== setLabels() =====\n\n";
  DATE = "No_Weight_2";
  gSystem->mkdir(Form("roots/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("../figs/2DFit_%s/CtauBkg", DATE.c_str()), kTRUE);

  kineLabel = std::string(getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh).Data());

  if (PRw == 1) fname = "PR";
  else if (PRw == 2) fname = "NP";
}

void CtauBkgFit::openInputFile()
{
  std::cout << "===== openInputFile() =====\n\n";

  fMass = new TFile(Form("roots/2DFit_%s/MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  fCErr = new TFile(Form("roots/2DFit_%s/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  fCRes = new TFile(Form("roots/2DFit_%s/CtauResResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  if (!fMass || fMass->IsZombie())
  {
    std::cerr << "CRITICAL: Mass fit result file could not be opened. Aborting.\n";
    exit(1);
  }
  if (!fCErr || fCErr->IsZombie())
  {
    std::cerr << "CRITICAL: Ctau error result file could not be opened. Aborting.\n";
    exit(1);
  }
  if (!fCRes || fCRes->IsZombie())
  {
    std::cerr << "CRITICAL: Ctau resolution result file could not be opened. Aborting.\n";
    exit(1);
  }
}

void CtauBkgFit::setupWorkspaceAndData()
{
  std::cout << "===== setupWorkspaceAndData() =====\n\n";

  RooDataSet *datasetMass = (RooDataSet *)fMass->Get("datasetMass");
  RooAddPdf *pdfMASS_Tot = (RooAddPdf *)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Bkg = (RooDataSet *)fCErr->Get("dataw_Bkg");
  RooAddPdf *GaussModel_Tot = (RooAddPdf *)fCRes->Get("GaussModel_Tot");
  RooAddPdf *pdfCTAUERR_Bkg = (RooAddPdf *)fCErr->Get("pdfCTAUERR_Bkg");

  ws = new RooWorkspace("workspace");

  ws->import(*datasetMass);
  ws->import(*pdfMASS_Tot);
  ws->import(*dataw_Bkg);
  ws->import(*GaussModel_Tot);
  ws->import(*pdfCTAUERR_Bkg);
}

void CtauBkgFit::setVariableRanges()
{
  std::cout << "===== setVariableRanges() =====\n\n";
  double ctauErrMin, ctauErrMax;
  ws->data("dataw_Bkg")->getRange(*ws->var("ctau3DErr"), ctauErrMin, ctauErrMax);

  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("bkgPlotRange", ctauLow, ctauHigh);
  
  ws->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
  ws->var("ctau3DErr")->setRange("bkgErrRange", ctauErrMin, ctauErrMax);
}

void CtauBkgFit::defineModel()
{
  std::cout << "===== defineModel() =====\n\n";

  // common variables
  ws->factory("zeroMean[0.0]");
  ws->factory("N_BkgCtau[1e5, 10, 1e7]");
  // ws->factory("b_Bkg[0.4, 0.1, 1]"); // NP b fraction of bkg

  // resolution model
  buildResModel();

  // decay submodels
  buildLeftModel();
  // buildCenterModel();
  buildRightModel();

  // combine models
  combineDecayModels();
}

void CtauBkgFit::buildResModel()
{
  std::cout << "===== buildResModel() =====\n\n";
  ws->factory("zeroBias[0.0]");
  // ws->factory("One[1.0]");

  ws->var("ctau1_CtauRes")->setConstant(kTRUE);
  ws->var("s1_CtauRes")->setConstant(kTRUE);

  // gauss1
  if (nGauss == 1)
    ws->factory("GaussModel::pdfCTAURES(ctau3D, ctau1_CtauRes, s1_CtauRes, zeroBias, ctau3DErr)");
  else
    ws->factory("GaussModel::ctauRes1(ctau3D, ctau1_CtauRes, s1_CtauRes, zeroBias, ctau3DErr)");

  if (nGauss >= 2) ws->factory("GaussModel::ctauRes2(ctau3D, ctau2_CtauRes, s2_CtauRes, zeroBias, ctau3DErr)");
  if (nGauss >= 3) ws->factory("GaussModel::ctauRes3(ctau3D, ctau3_CtauRes, s3_CtauRes, zeroBias, ctau3DErr)");
  if (nGauss == 4) ws->factory("GaussModel::ctauRes4(ctau3D, ctau4_CtauRes, s4_CtauRes, zeroBias, ctau3DErr)");

  // add models
  if (nGauss == 1)
  {
    // empty intentionally - nGauss = 1 case was handled above.
  } else if (nGauss == 2)
  {
    ws->var("ctau2_CtauRes")->setConstant(kTRUE);
    ws->var("rS21_CtauRes")->setConstant(kTRUE);
    ws->var("f_CtauRes")->setConstant(kTRUE);
    ws->factory("AddModel::pdfCTAURES({ctauRes1, ctauRes2}, {f_CtauRes})");
  }
  else if (nGauss == 3)
  {
    ws->var("ctau3_CtauRes")->setConstant(kTRUE);
    ws->var("rS32_CtauRes")->setConstant(kTRUE);
    ws->var("f2_CtauRes")->setConstant(kTRUE);
    ws->factory("AddModel::ctauRes32({ctauRes3, ctauRes2}, {f2_CtauRes})");
    ws->factory("AddModel::pdfCTAURES({ctauRes1, ctauRes32}, {f_CtauRes})");
  }
  else if (nGauss == 4)
  {
    ws->var("ctau4_CtauRes")->setConstant(kTRUE);
    ws->var("rS43_CtauRes")->setConstant(kTRUE);
    ws->var("f3_CtauRes")->setConstant(kTRUE);
    ws->factory("AddModel::ctauRes43({ctauRes4, ctauRes3}, {f3_CtauRes})");
    ws->factory("AddModel::ctauRes432({ctauRes43, ctauRes2}, {f2_CtauRes})");
    ws->factory("AddModel::pdfCTAURES({ctauRes1, ctauRes432}, {f_CtauRes})");
  }
  else
  {
    std::cerr << "\n[ERROR] nGauss should be 1 - 4, you typed: " << nGauss << std::endl;
    exit(1);
  }
}

void CtauBkgFit::buildRightModel()
{
  std::cout << "===== buildRightModel() =====\n\n";
  if (nExp_R == 1)
  {
    ws->factory("lambdaDSS_Bkg1[1.0, 0.001, 10.0]");

    ws->factory("Decay::pdfCTAUDSS(ctau3D, lambdaDSS_Bkg1, pdfCTAURES, RooDecay::SingleSided)");
  } else if (nExp_R == 2) {
    ws->factory("lambdaDSS_Bkg1[1.0, 0.001, 10.0]");
    ws->factory("rDSS12[2.0, 1.0, 10.0]");
    ws->factory("expr::lambdaDSS_Bkg2('lambdaDSS_Bkg1 * rDSS12', lambdaDSS_Bkg1, rDSS12)");
    ws->factory("fDSS12[0.5, 0.01, 1.0]");

    ws->factory("Decay::pdfCTAUDSS1(ctau3D, lambdaDSS_Bkg1, pdfCTAURES, RooDecay::SingleSided)");
    ws->factory("Decay::pdfCTAUDSS2(ctau3D, lambdaDSS_Bkg2, pdfCTAURES, RooDecay::SingleSided)");

    ws->factory("SUM::pdfCTAUDSS(fDSS12 * pdfCTAUDSS1, pdfCTAUDSS2)");
  } else {
    std::cerr << "[ERROR] nExp_R shoud be 1 - 2, you typed: " << nExp_R << "\n";
    exit(1);
  }
}

void CtauBkgFit::buildCenterModel()
{
  // Not use - dummy

  // std::cout << "===== buildCenterModel() =====\n\n";
  // if (nExp_C == 1) {
  //   ws->factory("lambdaDDS_Bkg1[0.5, 0.001, 10.0]");

  //   ws->factory("Decay::pdfCTAUDDS(ctau3D, lambdaDDS_Bkg1, pdfCTAURES, RooDecay::DoubleSided)");
  // } else if (nExp_C == 2) {
  //   ws->factory("lambdaDDS_Bkg1[0.5, 0.001, 10.0]");
  //   ws->factory("lambdaDDS_Bkg2[1.5, 0.001, 10.0]");
  //   ws->factory("fDDS12[0.4, 0.01, 1.0]");

  //   ws->factory("Decay::pdfCTAUDDS1(ctau3D, lambdaDDS_Bkg1, pdfCTAURES, RooDecay::DoubleSided)");
  //   ws->factory("Decay::pdfCTAUDDS2(ctau3D, lambdaDDS_Bkg2, pdfCTAURES, RooDecay::DoubleSided)");

  //   ws->factory("SUM::pdfCTAUDDS(fDDS12 * pdfCTAUDDS1, pdfCTAUDDS2)");
  // } else {
  //   std::cerr << "[ERROR] nExp_C shoud be 1 - 2, you typed: " << nExp_C << "\n";
  //   exit(1);
  // }
}

void CtauBkgFit::buildLeftModel()
{
  std::cout << "===== buildLeftModel() =====\n\n";
  if (nExp_L == 1) {
    ws->factory("lambdaDF_Bkg1[1.5, 0.001, 10.0]");

    ws->factory("Decay::pdfCTAUDF(ctau3D, lambdaDF_Bkg1, pdfCTAURES, RooDecay::Flipped)");
  } else if (nExp_L == 2) {
    ws->factory("lambdaDF_Bkg1[1.5, 0.001, 10.0]");
    ws->factory("rDF12[2.0, 1.0, 10.0]");
    ws->factory("expr::lambdaDF_Bkg2('lambdaDF_Bkg1 * rDF12', lambdaDF_Bkg1, rDF12)");
    ws->factory("fDF12[0.4, 0.01, 1.0]");

    ws->factory("Decay::pdfCTAUDF1(ctau3D, lambdaDF_Bkg1, pdfCTAURES, RooDecay::Flipped)");
    ws->factory("Decay::pdfCTAUDF2(ctau3D, lambdaDF_Bkg2, pdfCTAURES, RooDecay::Flipped)");

    ws->factory("SUM::pdfCTAUDF(fDF12 * pdfCTAUDF1, pdfCTAUDF2)");
  } else {
    std::cerr << "[ERROR] nExp_L shoud be 1 - 2, you typed: " << nExp_L << "\n";
    exit(1);
  }
}

void CtauBkgFit::combineDecayModels()
{
  std::cout << "===== combineDecayModels() =====\n\n";
  ws->factory("fDecayP[0.5, 0.0, 1.0]");
  ws->factory("fDecayM[0.25, 0.0, 1.0]");
  ws->factory("expr::fResCore('1.0 - fDecayP - fDecayM', fDecayP, fDecayM)");
  ws->factory("SUM::model(fResCore*pdfCTAURES, fDecayM*pdfCTAUDF, fDecayP*pdfCTAUDSS)");

  ws->factory("RooExtendPdf::pdfCTAU_Bkg_Tot(model, N_BkgCtau)");

  // --- legacy test ---
  // // ===== conditional model =====
  // // To check goodness of fit -> NOT USED FOR FIT HERE.
  // RooProdPdf pdfPR("pdfCTAU_BkgPR", "", *ws->pdf("pdfCTAUERR_Bkg"), Conditional(*ws->pdf("pdfCTAUCOND_BkgPR"), RooArgList(*ws->var("ctau3D"))));
  // ws->import(pdfPR);

  // RooProdPdf pdfNoPR("pdfCTAU_BkgNoPR", "", *ws->pdf("pdfCTAUERR_Bkg"), Conditional(*ws->pdf("pdfCTAUCOND_BkgNoPR"), RooArgList(*ws->var("ctau3D"))));
  // ws->import(pdfNoPR);

  // ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU_Bkg", "b_Bkg", "pdfCTAU_BkgNoPR", "pdfCTAU_BkgPR"));
  // RooAbsPdf *pdfCTAU_Bkg_Tot = new RooAddPdf("pdfCTAU_Bkg_Tot", "pdfCTAU_Bkg_Tot", RooArgList(*ws->pdf("pdfCTAU_Bkg")), RooArgList(*ws->var("N_BkgCtau")));

  // ws->import(*pdfCTAU_Bkg_Tot);

  // delete pdfCTAU_Bkg_Tot;
  // --- end of legacy test ---
}

void CtauBkgFit::initVar(const std::string &varName, double init, double low, double high)
{
  RooRealVar *var = ws->var(varName.c_str());
  if (!var)
  {
    std::cerr << "[ERROR] there is no variable:: " << varName << "\n";
    exit(1);
  }

  if (init < low || init > high)
  {
    std::cerr << "[ERROR] init value out of bounds for variable: " << varName << "\n";
    std::cerr << "        init = " << init << ", range = [" << low << ", " << high << "]\n";
    exit(1);
  }

  var->setVal(init);
  var->setMin(low);
  var->setMax(high);
  // var->setRange(low, high);
}

void CtauBkgFit::setConstVar(const std::string &varName, bool isConst, double value)
{
  if (value != 3096)
    ws->var(varName.c_str())->setVal(value);
  ws->var(varName.c_str())->setConstant(isConst);
}

void CtauBkgFit::performFit()
{
  std::cout << "===== performFit() =====\n\n";
  // set proper fit range
  TH1D *hTot = (TH1D *)ws->data("dataw_Bkg")->createHistogram(("hTot"), *ws->var("ctau3D"), Binning(nCtauBins, ctauLow, ctauHigh));
  if (ctauMin == -100)
    ctauMin = hTot->GetBinLowEdge(hTot->FindFirstBinAbove(1, 1));
  if (ctauMax == 100)
    ctauMax = hTot->GetBinLowEdge(hTot->FindLastBinAbove(2, 1)) + hTot->GetBinWidth(hTot->FindLastBinAbove(2, 1));
  delete hTot;
  ws->var("ctau3D")->setRange("bkgFitRange", ctauMin, ctauMax);

  // make a dataset for fitting
  RooDataSet *dataToFit = (RooDataSet *)ws->data("dataw_Bkg")->reduce(CutRange("bkgFitRange"));
  ws->import(*dataToFit, Rename("dataToFit"));

  // fit
  bool isWeighted = ws->data("dataToFit")->isWeighted();
  fitResult = ws->pdf("pdfCTAU_Bkg_Tot")->fitTo(*ws->data("dataToFit"), Save(), Range("bkgFitRange"), Extended(kTRUE), NumCPU(nCPU), PrintLevel(0), SumW2Error(isWeighted), RecoverFromUndefinedRegions(1.), PrintEvalErrors(-1), Strategy(2));
  fitResult->Print("V");

  delete dataToFit;
}

void CtauBkgFit::drawTextVar(const char *varName, const char *label, float xp, float yp, int textColor, int textSize)
{
  RooRealVar *var = ws->var(varName);
  if (!var)
    return;

  double val = var->getVal();
  double err = var->getError();
  double low = var->getMin();
  double high = var->getMax();
  bool isConst = var->isConstant();

  const double abs_epsilon = 1e-4; // boundary-val < 0.0001
  const double rel_epsilon = 1e-3; // 0.1 %
  const double minErr = 1e-4;

  // if boundary-val <0.0001 and diff/limit < 0.1%, parameter is stuck at the boundary
  bool atLowerLimit = (std::fabs(val - low) < abs_epsilon) || (low != 0 && std::fabs((val - low) / low) < rel_epsilon);
  bool atUpperLimit = (std::fabs(val - high) < abs_epsilon) || (high != 0 && std::fabs((val - high) / high) < rel_epsilon);

  TString text;
  if (isConst)
    text = Form("%s = %.4f (fixed)", label, val);
  else if (atLowerLimit || atUpperLimit)
    text = Form("%s = %.4f (limit)", label, val);
  else if (err < minErr)
    text = Form("%s = %.4f #pm < %.4f", label, val, minErr);
  else if (err > (high - low))
    text = Form("%s = %.4f #pm %.4f (unstable)", label, val, err);
  else
    text = Form("%s = %.4f #pm %.4f", label, val, err);

  drawText(text, xp, yp, textColor, textSize);
}

void CtauBkgFit::drawTextVarInt(const char *varName, const char *label, float xp, float yp, int textColor, int textSize)
{
  RooRealVar *var = ws->var(varName);
  if (!var)
    return;

  double val = var->getVal();
  double err = var->getError();
  double low = var->getMin();
  double high = var->getMax();
  bool isConst = var->isConstant();

  const double abs_epsilon = 1e-4; // boundary-val < 0.0001
  const double rel_epsilon = 1e-3; // 0.1 %
  const double minErr = 1e-4;

  // if boundary-val <0.0001 and diff/limit < 0.1%, parameter is stuck at the boundary
  bool atLowerLimit = (std::fabs(val - low) < abs_epsilon) || (low != 0 && std::fabs((val - low) / low) < rel_epsilon);
  bool atUpperLimit = (std::fabs(val - high) < abs_epsilon) || (high != 0 && std::fabs((val - high) / high) < rel_epsilon);

  TString text;
  if (isConst)
    text = Form("%s = %.f (fixed)", label, val);
  else if (atLowerLimit || atUpperLimit)
    text = Form("%s = %.f (limit)", label, val);
  else if (err < minErr)
    text = Form("%s = %.f #pm < %.f", label, val, minErr);
  else if (err > (high - low))
    text = Form("%s = %.f #pm < %.f (unstable)", label, val, minErr);
  else
    text = Form("%s = %.f #pm %.f", label, val, err);

  drawText(text, xp, yp, textColor, textSize);
}

void CtauBkgFit::drawPullPlot()
{
  std::cout << "===== Start drawPullPlot() =====\n";
  TCanvas *myCanvas = new TCanvas("myCanvas", "", 800, 800);
  myCanvas->cd();
  TPad *ctauPad = new TPad("ctauPad", "Main plot pad", 0.0, 0.3, 1.0, 1.0);
  ctauPad->SetTicks(1, 1);
  ctauPad->SetLeftMargin(0.15);
  ctauPad->SetRightMargin(0.07);
  ctauPad->SetTopMargin(0.1);
  ctauPad->SetBottomMargin(0.02);
  ctauPad->Draw();
  ctauPad->cd();
  gPad->SetLogy();

  RooPlot *ctauFrame = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(ctauLow, ctauHigh)); // bins
  ctauFrame->SetTitle("");

  // Todo: remove legacy?
  double normDSTot = 1.0;
  if (ws->data("dataToFit"))
  {
    normDSTot = ws->data("dataToFit")->sumEntries() / ws->data("dataToFit")->sumEntries();
  }
  double normBkg = 1.0;
  if (ws->data("dataw_Bkg"))
  {
    normBkg = ws->data("dataToFit")->sumEntries() * normDSTot / ws->data("dataw_Bkg")->sumEntries();
  }

  ws->pdf("pdfCTAU_Bkg_Tot")->setNormRange("ctaubkgFitRange");
  ctauFrame->updateNormVars(RooArgSet(*ws->var("ctau3D")));
  // ctauFrame->updateNormVars(RooArgSet(*ws->var("mass"), *ws->var("ctau3D"), *ws->var("ctau3DErr")));

  // plotting
  ws->data("dataToFit")->plotOn(ctauFrame, Name("data_points"), DataError(RooAbsData::SumW2)); // for normalization
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("total_pdf"), ProjWData(*ws->data("dataw_Bkg")), Range("bkgPlotRange"), NormRange("bkgFitRange"));
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("gaussCore"), Components(RooArgSet(*ws->pdf("pdfCTAURES"))), LineStyle(kDashed), LineColor(kGreen + 2), LineWidth(3), Range("bkgPlotRange"), NormRange("bkgFitRange"));

  if (nExp_L == 1)
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyL"), Components(RooArgSet(*ws->pdf("pdfCTAUDF"))), LineStyle(kDashed), LineColor(kMagenta + 1), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
  else if (nExp_L == 2) {
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyL1"), Components(RooArgSet(*ws->pdf("pdfCTAUDF1"))), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyL2"), Components(RooArgSet(*ws->pdf("pdfCTAUDF2"))), LineStyle(kDotted), LineColor(kMagenta + 3), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
  }

  if (nExp_R == 1)
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyR"), Components(RooArgSet(*ws->pdf("pdfCTAUDSS"))), LineStyle(kDashed), LineColor(kOrange + 1), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
  else if (nExp_R == 2)
  {
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyR1"), Components(RooArgSet(*ws->pdf("pdfCTAUDSS1"))), LineStyle(kDashed), LineColor(kOrange + 2), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyR2"), Components(RooArgSet(*ws->pdf("pdfCTAUDSS2"))), LineStyle(kDotted), LineColor(kOrange + 3), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
  }

  ws->data("dataToFit")->plotOn(ctauFrame, Name("data_over_pdfs"), DataError(RooAbsData::SumW2)); // draw dataset again over pdf (cosmetic)

  // set Y range
  TH1 *h_tmp = ws->data("dataw_Bkg")->createHistogram("h_tmp", *ws->var("ctau3D"), Binning(ctauFrame->GetNbinsX(), ctauFrame->GetXaxis()->GetXmin(), ctauFrame->GetXaxis()->GetXmax()));
  Double_t YMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
    if (h_tmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, h_tmp->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  // Yup = YMax*TMath::Power((YMax/0.1), 0.5);
  Yup = YMax * TMath::Power((YMax / 0.01), 0.5);
  Ydown = 0.001;

  // ctauFrame->GetYaxis()->SetRangeUser(10e-2, 10e7);
  ctauFrame->GetYaxis()->SetRangeUser(Ydown, Yup);
  ctauFrame->GetXaxis()->SetRangeUser(ctauLow, ctauHigh);
  ctauFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  ctauFrame->SetFillStyle(4000);
  ctauFrame->GetYaxis()->SetTitleOffset(1.43);
  ctauFrame->GetXaxis()->SetLabelSize(0);
  ctauFrame->GetXaxis()->SetTitleSize(0);
  ctauFrame->Draw();

  TLine *minline = new TLine(ctauMin, 0.0, ctauMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  ctauFrame->addObject(minline);
  TLine *maxline = new TLine(ctauMax, 0.0, ctauMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  ctauFrame->addObject(maxline);

  // ctauFrame->Draw();
  TLegend *leg_E = new TLegend(text_x + 0.3, text_y + 0.04, text_x + 0.43, text_y - 0.15);
  leg_E->SetTextSize(text_size);
  leg_E->SetTextFont(43);
  leg_E->SetBorderSize(0);
  leg_E->AddEntry(ctauFrame->findObject("data_points"), "Data_Bkg", "pe");
  leg_E->AddEntry(ctauFrame->findObject("total_pdf"), "Total PDF", "fl");
  leg_E->AddEntry(ctauFrame->findObject("gaussCore"), "Gauss Resolution", "fl");

  if (nExp_L == 1)
    leg_E->AddEntry(ctauFrame->findObject("deacyL"), "Decay Left", "fl");
  else if (nExp_L == 2) {
    leg_E->AddEntry(ctauFrame->findObject("deacyL1"), "Decay Left1", "fl");
    leg_E->AddEntry(ctauFrame->findObject("deacyL2"), "Decay Left2", "fl");
  }

  if (nExp_R == 1)
    leg_E->AddEntry(ctauFrame->findObject("deacyR"), "Decay Right", "fl");
  else if (nExp_R == 2) {
    leg_E->AddEntry(ctauFrame->findObject("deacyR1"), "Decay Right1", "fl");
    leg_E->AddEntry(ctauFrame->findObject("deacyR2"), "Decay Right2", "fl");
  }
  leg_E->Draw("same");

  // --- print parameters ---
  // left: kinematics
  int yCountLeft = 0;
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x+0.05, text_y - y_diff * yCountLeft++, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x+0.05, text_y - y_diff * yCountLeft++, text_color, text_size);
  else
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x+0.05, text_y - y_diff * yCountLeft++, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x+0.05, text_y - y_diff * yCountLeft++, text_color, text_size);

  // right: fit parameters
  int yCountRight = 0;

  // resolution
  drawTextVarInt("N_BkgCtau", "N_{Bkg}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);

  drawTextVar("s1_CtauRes", "s1_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  if (nGauss >= 2)
  {
    drawTextVar("rS21_CtauRes", "(s2/s1)_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f_CtauRes", "f1_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  }
  if (nGauss >= 3)
  {
    drawTextVar("rS32_CtauRes", "(s3/s2)_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f2_CtauRes", "f2_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  }
  if (nGauss == 4)
  {
    drawTextVar("rS43_CtauRes", "(s4/s3)_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f3_CtauRes", "f3_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  }

  // bkg decay models
  drawTextVar("fDecayP", "f_{DecayR}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  drawTextVar("fDecayM", "f_{DecayL}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);

  // left side
  drawTextVar("lambdaDF_Bkg1", "#lambda_{DL1}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  if (nExp_L >= 2)
  {
    drawTextVar("fDF12", "f_{DL12}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("rDF12", "(#lambda_{DL2}/#lambda_{DL1})", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  }

  // right
  drawTextVar("lambdaDSS_Bkg1", "#lambda_{DR1}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  if (nExp_R >= 2)
  {
    drawTextVar("fDSS12", "f_{DR12}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("rDSS12", "(#lambda_{DR2}/#lambda_{DR1})", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  }

  TPad *pullPad = new TPad("pullPad", "Pull plot pad", 0.0, 0.0, 1.0, 0.3);
  myCanvas->cd();
  pullPad->SetLeftMargin(0.15);
  pullPad->SetRightMargin(0.07);
  pullPad->SetTopMargin(0); // Upper and lower plot are joined
  // pullPad->SetTopMargin(0.05);
  pullPad->SetBottomMargin(0.35);
  pullPad->SetTicks(1, 1);
  pullPad->SetFillStyle(4000);
  pullPad->SetFrameFillStyle(4000);
  pullPad->Draw();
  pullPad->cd();

  RooPlot *frameTMP = (RooPlot *)ctauFrame->Clone("frameTMP");
  RooHist *pullHist = frameTMP->pullHist("data_points", "total_pdf", true);
  pullHist->SetMarkerSize(0.8);
  RooPlot *pullFrame = ws->var("ctau3D")->frame(Title("Pull Distribution"));
  pullFrame->addPlotable(pullHist, "P");
  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.08);
  pullFrame->GetYaxis()->SetLabelSize(0.08);
  pullFrame->GetYaxis()->SetRangeUser(-10, 10);
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  pullFrame->GetXaxis()->SetTitleOffset(1.05);
  pullFrame->GetXaxis()->SetLabelOffset(0.04);
  pullFrame->GetXaxis()->SetLabelSize(0.08);
  pullFrame->GetXaxis()->SetTitleSize(0.08);
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw();

  TLine *l1 = new TLine(ctauLow, 0, ctauHigh, 0);
  l1->SetLineStyle(1);
  l1->Draw("same");

  printChi2(ws, pullPad, frameTMP, fitResult, "ctau3D", "data_points", "total_pdf", nCtauBins, false);

  Double_t theNLL = fitResult->minNll();
  std::cout << " *** NLL : " << std::setprecision(15) << theNLL << "\n";

  myCanvas->Update();
  myCanvas->Draw();
  myCanvas->SaveAs(Form("figs/2DFit_%s/Bkg_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  myCanvas->SaveAs(Form("../figs/2DFit_%s/CtauBkg/Bkg_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  delete h_tmp;
  delete frameTMP;
  delete myCanvas;
}

void CtauBkgFit::drawRatioPlot()
{
  std::cout << "===== Start drawRatioPlot() =====\n";
  TCanvas *myCanvas = new TCanvas("myCanvas", "", 800, 800);
  myCanvas->cd();
  TPad *ctauPad = new TPad("ctauPad", "Main plot pad", 0.0, 0.3, 1.0, 1.0);
  ctauPad->SetTicks(1, 1);
  ctauPad->SetLeftMargin(0.15);
  ctauPad->SetRightMargin(0.07);
  ctauPad->SetTopMargin(0.1);
  ctauPad->SetBottomMargin(0.02);
  ctauPad->Draw();
  ctauPad->cd();
  gPad->SetLogy();

  RooPlot *ctauFrame = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(ctauLow, ctauHigh)); // bins
  ctauFrame->SetTitle("");

  // Todo: remove legacy?
  double normDSTot = 1.0;
  if (ws->data("dataToFit"))
  {
    normDSTot = ws->data("dataToFit")->sumEntries() / ws->data("dataToFit")->sumEntries();
  }
  double normBkg = 1.0;
  if (ws->data("dataw_Bkg"))
  {
    normBkg = ws->data("dataToFit")->sumEntries() * normDSTot / ws->data("dataw_Bkg")->sumEntries();
  }

  ws->pdf("pdfCTAU_Bkg_Tot")->setNormRange("ctaubkgFitRange");
  ctauFrame->updateNormVars(RooArgSet(*ws->var("ctau3D")));
  // ctauFrame->updateNormVars(RooArgSet(*ws->var("mass"), *ws->var("ctau3D"), *ws->var("ctau3DErr")));

  // plotting
  ws->data("dataToFit")->plotOn(ctauFrame, Name("data_points"), DataError(RooAbsData::SumW2)); // for normalization
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("total_pdf"), ProjWData(*ws->data("dataw_Bkg")), Range("bkgPlotRange"), NormRange("bkgFitRange"));
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("gaussCore"), Components(RooArgSet(*ws->pdf("pdfCTAURES"))), LineStyle(kDashed), LineColor(kGreen + 2), LineWidth(3), Range("bkgPlotRange"), NormRange("bkgFitRange"));

  if (nExp_L == 1)
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyL"), Components(RooArgSet(*ws->pdf("pdfCTAUDF"))), LineStyle(kDashed), LineColor(kMagenta + 1), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
  else if (nExp_L == 2)
  {
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyL1"), Components(RooArgSet(*ws->pdf("pdfCTAUDF1"))), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyL2"), Components(RooArgSet(*ws->pdf("pdfCTAUDF2"))), LineStyle(kDotted), LineColor(kMagenta + 3), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
  }

  if (nExp_R == 1)
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyR"), Components(RooArgSet(*ws->pdf("pdfCTAUDSS"))), LineStyle(kDashed), LineColor(kOrange + 1), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
  else if (nExp_R == 2)
  {
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyR1"), Components(RooArgSet(*ws->pdf("pdfCTAUDSS1"))), LineStyle(kDashed), LineColor(kOrange + 2), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("deacyR2"), Components(RooArgSet(*ws->pdf("pdfCTAUDSS2"))), LineStyle(kDotted), LineColor(kOrange + 3), LineWidth(2), Range("bkgPlotRange"), NormRange("bkgFitRange"));
  }

  ws->data("dataToFit")->plotOn(ctauFrame, Name("data_over_pdfs"), DataError(RooAbsData::SumW2)); // draw dataset again over pdf (cosmetic)

  // set Y range
  TH1 *h_tmp = ws->data("dataw_Bkg")->createHistogram("h_tmp", *ws->var("ctau3D"), Binning(ctauFrame->GetNbinsX(), ctauFrame->GetXaxis()->GetXmin(), ctauFrame->GetXaxis()->GetXmax()));
  Double_t YMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
    if (h_tmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, h_tmp->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  // Yup = YMax*TMath::Power((YMax/0.1), 0.5);
  Yup = YMax * TMath::Power((YMax / 0.01), 0.5);
  Ydown = 0.001;

  // ctauFrame->GetYaxis()->SetRangeUser(10e-2, 10e7);
  ctauFrame->GetYaxis()->SetRangeUser(Ydown, Yup);
  ctauFrame->GetXaxis()->SetRangeUser(ctauLow, ctauHigh);
  ctauFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  ctauFrame->SetFillStyle(4000);
  ctauFrame->GetYaxis()->SetTitleOffset(1.43);
  ctauFrame->GetXaxis()->SetLabelSize(0);
  ctauFrame->GetXaxis()->SetTitleSize(0);
  ctauFrame->Draw();

  TLine *minline = new TLine(ctauMin, 0.0, ctauMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  ctauFrame->addObject(minline);
  TLine *maxline = new TLine(ctauMax, 0.0, ctauMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  ctauFrame->addObject(maxline);

  // ctauFrame->Draw();
  TLegend *leg_E = new TLegend(text_x + 0.3, text_y + 0.04, text_x + 0.43, text_y - 0.15);
  leg_E->SetTextSize(text_size);
  leg_E->SetTextFont(43);
  leg_E->SetBorderSize(0);
  leg_E->AddEntry(ctauFrame->findObject("data_points"), "Data_Bkg", "pe");
  leg_E->AddEntry(ctauFrame->findObject("total_pdf"), "Total PDF", "fl");
  leg_E->AddEntry(ctauFrame->findObject("gaussCore"), "Gauss Resolution", "fl");

  if (nExp_L == 1)
    leg_E->AddEntry(ctauFrame->findObject("deacyL"), "Decay Left", "fl");
  else if (nExp_L == 2)
  {
    leg_E->AddEntry(ctauFrame->findObject("deacyL1"), "Decay Left1", "fl");
    leg_E->AddEntry(ctauFrame->findObject("deacyL2"), "Decay Left2", "fl");
  }

  if (nExp_R == 1)
    leg_E->AddEntry(ctauFrame->findObject("deacyR"), "Decay Right", "fl");
  else if (nExp_R == 2)
  {
    leg_E->AddEntry(ctauFrame->findObject("deacyR1"), "Decay Right1", "fl");
    leg_E->AddEntry(ctauFrame->findObject("deacyR2"), "Decay Right2", "fl");
  }
  leg_E->Draw("same");

  // --- print parameters ---
  // left: kinematics
  int yCountLeft = 0;
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x + 0.05, text_y - y_diff * yCountLeft++, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x + 0.05, text_y - y_diff * yCountLeft++, text_color, text_size);
  else
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x + 0.05, text_y - y_diff * yCountLeft++, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x + 0.05, text_y - y_diff * yCountLeft++, text_color, text_size);

  // right: fit parameters
  int yCountRight = 0;

  // resolution
  drawTextVarInt("N_BkgCtau", "N_{Bkg}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);

  drawTextVar("s1_CtauRes", "s1_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  if (nGauss >= 2)
  {
    drawTextVar("rS21_CtauRes", "(s2/s1)_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f_CtauRes", "f1_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  }
  if (nGauss >= 3)
  {
    drawTextVar("rS32_CtauRes", "(s3/s2)_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f2_CtauRes", "f2_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  }
  if (nGauss == 4)
  {
    drawTextVar("rS43_CtauRes", "(s4/s3)_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f3_CtauRes", "f3_{Res}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  }

  // bkg decay models
  drawTextVar("fDecayP", "f_{DecayR}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  drawTextVar("fDecayM", "f_{DecayL}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);

  // left side
  drawTextVar("lambdaDF_Bkg1", "#lambda_{DL1}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  if (nExp_L >= 2)
  {
    drawTextVar("fDF12", "f_{DL12}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("rDF12", "(#lambda_{DL2}/#lambda_{DL1})", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  }

  // right
  drawTextVar("lambdaDSS_Bkg1", "#lambda_{DR1}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  if (nExp_R >= 2)
  {
    drawTextVar("fDSS12", "f_{DR12}", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("rDSS12", "(#lambda_{DR2}/#lambda_{DR1})", text_x + 0.55, text_y - y_diff * yCountRight++, text_color, text_size);
  }

  TPad *ratioPad = new TPad("ratioPad", "", 0.0, 0.0, 1.0, 0.3);
  myCanvas->cd();
  ratioPad->SetLeftMargin(0.15);
  ratioPad->SetRightMargin(0.07);
  ratioPad->SetTopMargin(0); // Upper and lower plot are joined
  // ratioPad->SetTopMargin(0.05);
  ratioPad->SetBottomMargin(0.35);
  ratioPad->SetTicks(1, 1);
  ratioPad->SetFillStyle(4000);
  ratioPad->SetFrameFillStyle(4000);
  ratioPad->Draw();
  ratioPad->cd();

  RooHist *h_dataPoints = (RooHist *)ctauFrame->findObject("data_points");
  RooCurve *curveRatio = (RooCurve *)ctauFrame->findObject("total_pdf");

  RooPlot *ratioFrame = ws->var("ctau3D")->frame(Title("Ratio"), Range(ctauLow, ctauHigh));
  ratioFrame->SetTitle("");

  TGraphAsymmErrors *g_ratio = new TGraphAsymmErrors();
  int point_idx = 0;
  for (int i = 0; i < h_dataPoints->GetN(); ++i)
  {
    double x, y;
    h_dataPoints->GetPoint(i, x, y);
    double model_val = curveRatio->Eval(x);
    if (model_val > 1e-9 && y > 0)
    {
      double ratio = y / model_val;
      g_ratio->SetPoint(point_idx, x, ratio);
      double err_y_low = h_dataPoints->GetErrorYlow(i);
      double err_y_high = h_dataPoints->GetErrorYhigh(i);
      g_ratio->SetPointError(point_idx, 0, 0, err_y_low / model_val, err_y_high / model_val);
      point_idx++;
    }
  }

  ratioFrame->SetTitle("");
  ratioFrame->SetTitleSize(0);
  ratioFrame->GetYaxis()->SetTitleOffset(0.3);
  ratioFrame->GetYaxis()->SetTitle("Data / Fit");
  ratioFrame->GetYaxis()->SetTitleSize(0.08);
  ratioFrame->GetYaxis()->SetLabelSize(0.08);
  ratioFrame->GetYaxis()->SetRangeUser(0, 3);
  ratioFrame->GetYaxis()->CenterTitle();

  ratioFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  ratioFrame->GetXaxis()->SetTitleOffset(1.15);
  ratioFrame->GetXaxis()->SetLabelOffset(0.04);
  ratioFrame->GetXaxis()->SetLabelSize(0.08);
  ratioFrame->GetXaxis()->SetTitleSize(0.08);
  ratioFrame->GetXaxis()->CenterTitle();

  ratioFrame->GetYaxis()->SetTickSize(0.04);
  ratioFrame->GetYaxis()->SetNdivisions(404);
  ratioFrame->GetXaxis()->SetTickSize(0.03);
  ratioFrame->Draw();

  TLine *line_at_1 = new TLine(ctauLow, 1.0, ctauHigh, 1.0);
  line_at_1->SetLineColor(kRed);
  line_at_1->SetLineStyle(kDashed);
  line_at_1->Draw("same");

  g_ratio->SetMarkerStyle(kFullCircle);
  g_ratio->SetMarkerSize(0.7);
  g_ratio->Draw("P SAME");

  myCanvas->Update();
  myCanvas->Draw();
  myCanvas->SaveAs(Form("figs/2DFit_%s/Bkg_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_ratio.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  myCanvas->SaveAs(Form("../figs/2DFit_%s/CtauBkg/Bkg_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_ratio.pdf", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  delete h_tmp;
  delete myCanvas;
}

void CtauBkgFit::saveResults()
{
  std::cout << "===== saveResults() =====\n\n";
  // TFile *outputFile = new TFile(Form("roots/2DFit_%s/CtauBkgResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP), "recreate");

  // fitResult->Write();
  // ws->pdf("pdfCTAU_Bkg_Tot")->Write("pdfCTAU_Bkg_Tot");

  // outputFile->Close();
  // delete outputFile;
}