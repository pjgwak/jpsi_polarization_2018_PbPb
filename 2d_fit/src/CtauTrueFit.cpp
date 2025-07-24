#include "CtauTrueFit.h"
#include <iostream>
// #include <algorithm>
#include "TStyle.h"
#include "RooMsgService.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooTruthModel.h"
#include "RooDecay.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "TCanvas.h"
#include "TPad.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "TH1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLegend.h"
#include "TLine.h"
#include "RooHist.h"

#include "../headers/polarizationUtilities.h"

using std::cout; using std::endl;
using namespace RooFit;

CtauTrueFit::CtauTrueFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
                       float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP)
    : ptLow(ptLow), ptHigh(ptHigh),
      yLow(yLow), yHigh(yHigh),
      cLow(cLow), cHigh(cHigh),
      cosLow(cosLow), cosHigh(cosHigh),
      PR(PR), PRw(PRw),
      fEffW(fEffW), fAccW(fAccW), isPtW(isPtW), isTnP(isTnP)
{
  gStyle->SetEndErrorSize(0);

  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
}

CtauTrueFit::~CtauTrueFit()
{
  RooMsgService::instance().getStream(0).addTopic(Caching);
  RooMsgService::instance().getStream(1).addTopic(Caching);
  RooMsgService::instance().getStream(0).addTopic(Plotting);
  RooMsgService::instance().getStream(1).addTopic(Plotting);
  RooMsgService::instance().getStream(0).addTopic(Integration);
  RooMsgService::instance().getStream(1).addTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);

  delete ws;
  delete fitResult;
  delete inputFile;
}

void CtauTrueFit::init()
{
  std::cout << "===== init() =====\n\n";
  setLabels();
  openInputFile();
  createKinematicCut();
  setupWorksapceAndData();
  setVariableRanges();
  defineModel();
}
void CtauTrueFit::run()
{
  std::cout << "===== run() =====\n\n";
  performFit();
  makePlot();
  saveResults();
}

void CtauTrueFit::setLabels()
{
  std::cout << "===== setLabels() =====\n\n";
  DATE = "No_Weight_2";
  gSystem->mkdir(Form("roots/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("../figs/2DFit_%s/CtauTrue", DATE.c_str()), kTRUE);

  TString kineTmp = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
  kineLabel = std::string(kineTmp.Data());
}


void CtauTrueFit::openInputFile()
{
  std::cout << "===== openInputFile() =====\n\n";
  inputFile = new TFile(inputFilePath.c_str());

  if (!inputFile || inputFile->IsZombie())
  {
    std::cerr << "CRITICAL: Input file could not be opened. Aborting.\n";
    exit(1);
  }
}

void CtauTrueFit::createKinematicCut()
{
  std::cout << "===== createKinematicCut() =====\n\n";
  TString cut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5", ptLow, ptHigh, yLow, yHigh);

  // It doesn't use the accCut
  // TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
  // cut += accCut;

  kineCutMC = std::string(cut.Data());
}

void CtauTrueFit::setupWorksapceAndData()
{
  std::cout << "===== setupWorksapceAndData() =====\n\n";
  ws = new RooWorkspace("workspace");

  RooDataSet *datasetMC = (RooDataSet *)inputFile->Get("dataset");
  // Todo: consider gen weight?
  RooDataSet *reducedDS_MC = (RooDataSet *)datasetMC->reduce(kineCutMC.c_str());
  reducedDS_MC->SetName("reducedDS_MC");

  ws->import(*reducedDS_MC);
  delete reducedDS_MC;
}

void CtauTrueFit::setVariableRanges()
{
  std::cout << "===== setVariableRanges() =====\n\n";
  RooRealVar *ctauVar = ws->var("ctau3D");

  ctauVar->setRange(ctau3DMin, ctau3DMax);
  ctauVar->setRange("truePlotRange", ctau3DMin, ctauHigh);
  ctauVar->setRange("trueFitRange", ctau3DMin, ctau3DMax);
}

void CtauTrueFit::defineModel()
{
  std::cout << "===== defineModel() =====\n\n";

  if (nExp == 1) buildTrue1Expo();
  else if (nExp == 2) buildTrue2Expo();
  else if (nExp == 3) buildTrue3Expo();
  else if (nExp == 4) buildTrue4Expo();
  else {
    std::cout << "[ERROR] failed to build a " << nExp << " exponential decay model\n";
    std::cout << "possible nExp: 1 - 4\n";
    exit(1);
  }
}

void CtauTrueFit::buildTrue1Expo()
{
  std::cout << "===== buildTrue1Expo() =====\n\n";
  RooRealVar *ctau = ws->var("ctau3D");
  RooTruthModel deltaFcn("deltaFcn_CTAUTRUE", "", *ctau);

  RooRealVar lambdaDSS("lambdaDSS", "", 0.1, 0.01, 5);
  RooDecay pdfCTAUTRUEDSS1("pdfCTAUTRUEDSS1", "", *ctau, lambdaDSS, deltaFcn, RooDecay::SingleSided);

  RooRealVar N_Jpsi_MC("N_Jpsi_MC", "", 500000, 100000, 1000000);
  RooExtendPdf total("TrueModel_Tot", "", pdfCTAUTRUEDSS1, N_Jpsi_MC);
  ws->import(total);
}

void CtauTrueFit::buildTrue2Expo()
{
  std::cout << "===== buildTrue2Expo() =====\n\n";

  RooRealVar *ctau = ws->var("ctau3D");
  RooTruthModel deltaFcn("deltaFcn_CTAUTRUE", "", *ctau);

  // RooRealVar lambdaDSS("lambdaDSS", "", 0.1, 0.01, 5);
  // RooRealVar lambdaDSS2("lambdaDSS2", "", 0.2, 0.01, 5); // old style - two correlated parameters

  RooRealVar lambdaDSS("lambdaDSS", "base decay const", 0.1, 0.01, 5.0);
  RooRealVar r_lambda2("r_lambda2", "log(lambda2/lambda1)", 0.0, -2.0, 2.0); // ln(lambda2 / lambda1)
  RooFormulaVar lambdaDSS2("lambdaDSS2", "@0 * exp(@1)", RooArgList(lambdaDSS, r_lambda2));

  RooDecay pdfCTAUTRUEDSS1("pdfCTAUTRUEDSS1", "", *ctau, lambdaDSS, deltaFcn, RooDecay::SingleSided);
  RooDecay pdfCTAUTRUEDSS2("pdfCTAUTRUEDSS2", "", *ctau, lambdaDSS2, deltaFcn, RooDecay::SingleSided);

  RooRealVar fDSS("fDSS", "fraction of decay1", 0.6, 0.01, 1.0);
  RooAddPdf pdfCTAUTRUE("pdfCTAUTRUE", "", RooArgList(pdfCTAUTRUEDSS1, pdfCTAUTRUEDSS2), RooArgList(fDSS));

  RooRealVar N_Jpsi_MC("N_Jpsi_MC", "", 500000, 100000, 1000000);
  RooExtendPdf total("TrueModel_Tot", "", pdfCTAUTRUE, N_Jpsi_MC);

  ws->import(total);
}

void CtauTrueFit::buildTrue3Expo()
{
  std::cout << "===== buildTrue3Expo() =====\n\n";
  RooRealVar *ctau = ws->var("ctau3D");
  RooTruthModel deltaFcn("deltaFcn_CTAUTRUE", "", *ctau);

  // RooRealVar lambdaDSS("lambdaDSS", "", 0.1, 0.01, 5);
  // RooRealVar lambdaDSS2("lambdaDSS2", "", 0.2, 0.01, 5);
  // RooRealVar lambdaDSS3("lambdaDSS3", "", 0.3, 0.01, 5);

  RooRealVar lambdaDSS("lambdaDSS", "", 0.1, 0.01, 5.0);
  RooRealVar r_lambda2("r_lambda2", "", 0.0, 0, 2.0);
  RooRealVar r_lambda3("r_lambda3", "", 0.0, 0, 2.0);

  RooFormulaVar lambdaDSS2("lambdaDSS2", "@0 * exp(@1)", RooArgList(lambdaDSS, r_lambda2));
  RooFormulaVar lambdaDSS3("lambdaDSS3", "@0 * exp(@1)", RooArgList(r_lambda2, r_lambda3));

  RooDecay pdf1("pdfCTAUTRUEDSS1", "", *ctau, lambdaDSS, deltaFcn, RooDecay::SingleSided);
  RooDecay pdf2("pdfCTAUTRUEDSS2", "", *ctau, lambdaDSS2, deltaFcn, RooDecay::SingleSided);
  RooDecay pdf3("pdfCTAUTRUEDSS3", "", *ctau, lambdaDSS3, deltaFcn, RooDecay::SingleSided);

  RooRealVar fDSS("fDSS", "fraction of pdf1", 0.4, 0.01, 1.0);
  RooRealVar f2_DSS("f2_DSS", "fraction of pdf2", 0.3, 0.01, 1.0);

  RooAddPdf sum23("pdfCTAUTRUE_23", "", RooArgList(pdf2, pdf3), RooArgList(f2_DSS));
  RooAddPdf pdfCTAUTRUE("pdfCTAUTRUE", "", RooArgList(pdf1, sum23), RooArgList(fDSS));

  RooRealVar N_Jpsi_MC("N_Jpsi_MC", "", 500000, 100000, 1000000);
  RooExtendPdf total("TrueModel_Tot", "", pdfCTAUTRUE, N_Jpsi_MC);
  ws->import(total);
}

void CtauTrueFit::buildTrue4Expo()
{
  std::cout << "===== buildTrue4Expo() =====\n\n";
  RooRealVar *ctau = ws->var("ctau3D");
  RooTruthModel deltaFcn("deltaFcn_CTAUTRUE", "", *ctau);

  // RooRealVar lambdaDSS("lambdaDSS", "", 0.1, 0.01, 5);
  // RooRealVar lambdaDSS2("lambdaDSS2", "", 0.2, 0.01, 5);
  // RooRealVar lambdaDSS3("lambdaDSS3", "", 0.3, 0.01, 5);
  // RooRealVar lambdaDSS4("lambdaDSS4", "", 0.4, 0.01, 5);

  RooRealVar lambdaDSS("lambdaDSS", "", 0.1, 0.01, 5.0);
  RooRealVar r_lambda2("r_lambda2", "", 0.0, -2.0, 2.0);
  RooRealVar r_lambda3("r_lambda3", "", 0.5, -2.0, 2.0);
  RooRealVar r_lambda4("r_lambda4", "", 1.0, -2.0, 2.0);

  RooFormulaVar lambdaDSS2("lambdaDSS2", "@0 * exp(@1)", RooArgList(lambdaDSS, r_lambda2));
  RooFormulaVar lambdaDSS3("lambdaDSS3", "@0 * exp(@1)", RooArgList(lambdaDSS2, r_lambda3));
  RooFormulaVar lambdaDSS4("lambdaDSS4", "@0 * exp(@1)", RooArgList(lambdaDSS3, r_lambda4));

  RooDecay pdf1("pdfCTAUTRUEDSS1", "", *ctau, lambdaDSS, deltaFcn, RooDecay::SingleSided);
  RooDecay pdf2("pdfCTAUTRUEDSS2", "", *ctau, lambdaDSS2, deltaFcn, RooDecay::SingleSided);
  RooDecay pdf3("pdfCTAUTRUEDSS3", "", *ctau, lambdaDSS3, deltaFcn, RooDecay::SingleSided);
  RooDecay pdf4("pdfCTAUTRUEDSS4", "", *ctau, lambdaDSS4, deltaFcn, RooDecay::SingleSided);

  RooRealVar fDSS("fDSS", "", 0.3, 0.01, 1.0);
  RooRealVar f2_DSS("f2_DSS", "", 0.3, 0.01, 1.0);
  RooRealVar f3_DSS("f3_DSS", "", 0.3, 0.01, 1.0);

  RooAddPdf sum34("pdfCTAUTRUE_34", "", RooArgList(pdf3, pdf4), RooArgList(f3_DSS));
  RooAddPdf sum234("pdfCTAUTRUE_234", "", RooArgList(pdf2, sum34), RooArgList(f2_DSS));
  RooAddPdf pdfCTAUTRUE("pdfCTAUTRUE", "", RooArgList(pdf1, sum234), RooArgList(fDSS));

  RooRealVar N_Jpsi_MC("N_Jpsi_MC", "", 500000, 100000, 1000000);
  RooExtendPdf total("TrueModel_Tot", "", pdfCTAUTRUE, N_Jpsi_MC);
  ws->import(total);
}

void CtauTrueFit::initVar(const std::string &varName, double init, double low, double high)
{
  if (init < low || init > high)
  {
    std::cerr << "[ERROR] init value out of bounds for variable: " << varName << "\n";
    std::cerr << "        init = " << init << ", range = [" << low << ", " << high << "]\n";
    exit(1);
  }

  ws->var(varName.c_str())->setVal(init);
  ws->var(varName.c_str())->setRange(low, high);
}

void CtauTrueFit::performFit()
{
  std::cout << "===== performFit() =====\n\n";
  RooDataSet *dataToFit = (RooDataSet *)ws->data("reducedDS_MC")->reduce(Form("ctau3D>=%.f&&ctau3D<%.f", ctau3DMin, ctau3DMax));
  ws->import(*dataToFit, Rename("dataToFit"));

  // bool isWeighted = ws->data("dataToFit")->isWeighted();
  fitResult = ws->pdf("TrueModel_Tot")->fitTo(*ws->data("dataToFit"), Save(), SumW2Error(false), Extended(kTRUE), NumCPU(nCPU), Strategy(2), Range("trueFitRange"), PrintLevel(-1));
  // Range("trueFitRange") makes fit be better

  fitResult->Print("V");
}

void CtauTrueFit::drawTextVar(const char *varName, const char *label, float xp, float yp, int textColor, int textSize)
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

void CtauTrueFit::drawTextVarInt(const char *varName, const char *label, float xp, float yp, int textColor, int textSize)
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

void CtauTrueFit::makePlot()
{
  std::cout << "===== makePlot() =====\n\n";
  ws->pdf("TrueModel_Tot")->setNormRange("trueFitRange");

  // add local label
  TString bCont;
  if (PR == 0) bCont = "Prompt";
  else if (PR == 1) bCont = "NonPrompt";
  else if (PR == 2) bCont = "Inclusive";


  TCanvas *myCanvas = new TCanvas("myCanvas", "", 800, 800);
  myCanvas->cd();

  // make true pad
  TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
  pad1->SetBottomMargin(0);
  pad1->SetTicks(1, 1);
  pad1->Draw();

  pad1->cd();

  // true frame
  RooPlot *ctauFrame = ws->var("ctau3D")->frame(Bins(nCtauTrueBins), Range("truePlotRange"));

  // draw dataset and pdfs
  ws->data("dataToFit")->plotOn(ctauFrame, Name("dataToFit"), DataError(RooAbsData::SumW2));

  ws->pdf("TrueModel_Tot")->plotOn(ctauFrame, Range("truePlotRange"), NormRange("trueFitRange"), LineColor(kBlack), Name("TrueModel_Tot"));
  ws->pdf("TrueModel_Tot")->plotOn(ctauFrame, Name("comp1"), Components(*ws->pdf("pdfCTAUTRUEDSS1")), LineStyle(kDashed), LineColor(kGreen + 2), Range("truePlotRange"), NormRange("trueFitRange"));
  if (nExp >= 2)
    ws->pdf("TrueModel_Tot")->plotOn(ctauFrame, Name("comp2"), Components(*ws->pdf("pdfCTAUTRUEDSS2")), LineStyle(kDashed), LineColor(kAzure - 4), Range("truePlotRange"), NormRange("trueFitRange"));
  if (nExp >= 3)
    ws->pdf("TrueModel_Tot")->plotOn(ctauFrame, Name("comp3"), Components(*ws->pdf("pdfCTAUTRUEDSS3")), LineStyle(kDashed), LineColor(kOrange + 1), Range("truePlotRange"), NormRange("trueFitRange"));
  if (nExp >= 4)
    ws->pdf("TrueModel_Tot")->plotOn(ctauFrame, Name("comp4"), Components(*ws->pdf("pdfCTAUTRUEDSS4")), LineStyle(kDashed), LineColor(kViolet + 1), Range("truePlotRange"), NormRange("trueFitRange"));

  // set y-axis range
  gPad->SetLogy();
  TH1 *h_tmp = ws->data("dataToFit")->createHistogram("h_tmp", *ws->var("ctau3D"), Binning(ctauFrame->GetNbinsX(), ctauFrame->GetXaxis()->GetXmin(), ctauFrame->GetXaxis()->GetXmax()));
  Double_t YMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
    if (h_tmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, h_tmp->GetBinContent(i));
  double Ydown;
  double Yup;
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
  if (Ydown <= 0)
    Ydown = 0.1;
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));
  ctauFrame->GetYaxis()->SetRangeUser(Ydown, Yup);
  delete h_tmp;

  // cosmetics
  ctauFrame->SetFillStyle(4000);
  ctauFrame->GetYaxis()->SetTitleOffset(1.43);
  ctauFrame->GetXaxis()->SetLabelSize(0); // to hide x-axis label
  ctauFrame->GetXaxis()->SetTitleSize(0);
  ctauFrame->Draw();

  // legend
  TLegend *leg = new TLegend(text_x + 0.29, text_y + 0.03, text_x + 0.39, text_y - 0.17);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(ctauFrame->findObject("dataToFit"), "Gen ctau3D", "pe");
  leg->AddEntry(ctauFrame->findObject("TrueModel_Tot"), "Total fit", "l");
  leg->AddEntry(ctauFrame->findObject("comp1"), "Decay Exp1", "l");
  if (nExp >= 2) leg->AddEntry(ctauFrame->findObject("comp2"), "Decay Exp2", "l");
  if (nExp >= 3) leg->AddEntry(ctauFrame->findObject("comp3"), "Decay Exp3", "l");
  if (nExp >= 4) leg->AddEntry(ctauFrame->findObject("comp4"), "Decay Exp4", "l");
  leg->Draw();

  // print parameters
  int yCountLeft = 0;
  int yCountRight = 0;

  // kinematics (left side)
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);

  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);
  
  // parameters
  drawTextVarInt("N_Jpsi_MC", "N_{J/#psi}^{MC}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  drawTextVar("lambdaDSS", "#lambda_{1}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  if (nExp >= 2) {
    drawTextVar("r_lambda2", "ln(#lambda_{2}/#lambda_{1})", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("fDSS", "f_{1}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  }
  if (nExp >= 3) {
    drawTextVar("r_lambda3", "ln(#lambda_{3}/#lambda_{2})", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f2_DSS", "f_{2}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  }
  if (nExp >= 4) {
    drawTextVar("r_lambda3", "ln(#lambda_{4}/#lambda_{3})", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f3_DSS", "f_{3}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  }



  // pull part
  TPad *pullPad = new TPad("pullPad", "", 0.0, 0.0, 1.0, 0.3);
  myCanvas->cd();
  pullPad->Draw();
  pullPad->cd();
  pullPad->SetTopMargin(0); // Upper and lower plot are joined
  pullPad->SetBottomMargin(0.67);
  pullPad->SetBottomMargin(0.4);
  pullPad->SetFillStyle(4000);
  pullPad->SetFrameFillStyle(4000);
  pullPad->SetTicks(1, 1);

  // to ignore empty residual bin warning - it's not a matter.
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);

  RooPlot *pullTmp = (RooPlot *)ctauFrame->Clone("pullTmp");
  RooHist *pullHist = pullTmp->pullHist("dataToFit", "TrueModel_Tot", true);
  pullHist->SetMarkerSize(0.8);
  RooPlot *pullFrame = ws->var("ctau3D")->frame(Title(""), Bins(nCtauTrueBins), Range(ctau3DMin, ctauHigh));
  pullFrame->addPlotable(pullHist, "PX");

  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.08);
  pullFrame->GetYaxis()->SetLabelSize(0.08);
  pullFrame->GetYaxis()->SetRangeUser(-9.5, 9.5);
  pullFrame->GetYaxis()->CenterTitle();
  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);

  pullFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  pullFrame->GetXaxis()->SetTitleOffset(1.55);
  pullFrame->GetXaxis()->SetLabelOffset(0.04);
  pullFrame->GetXaxis()->SetLabelSize(0.08);
  pullFrame->GetXaxis()->SetTitleSize(0.08);
  pullFrame->GetXaxis()->CenterTitle();
  pullFrame->Draw();

  // ===== horizontoal lines on pull hist ===== //
  TLine *line_at_2 = new TLine(pullFrame->GetXaxis()->GetXmin(), 2, pullFrame->GetXaxis()->GetXmax(), 2);
  line_at_2->SetLineColor(kBlack);
  line_at_2->SetLineStyle(3);
  line_at_2->SetLineWidth(1);
  line_at_2->Draw();

  TLine *line_at_4 = new TLine(pullFrame->GetXaxis()->GetXmin(), 4, pullFrame->GetXaxis()->GetXmax(), 4);
  line_at_4->SetLineColor(kBlack);
  line_at_4->SetLineStyle(3);
  line_at_4->SetLineWidth(2);
  line_at_4->Draw();

  TLine *line_at_8 = new TLine(pullFrame->GetXaxis()->GetXmin(), 8, pullFrame->GetXaxis()->GetXmax(), 8);
  line_at_8->SetLineColor(kBlack);
  line_at_8->SetLineStyle(3);
  line_at_8->SetLineWidth(1);
  line_at_8->Draw();

  // negative lines
  TLine *line_at_n2 = new TLine(pullFrame->GetXaxis()->GetXmin(), -2, pullFrame->GetXaxis()->GetXmax(), -2);
  line_at_n2->SetLineColor(kBlack);
  line_at_n2->SetLineStyle(3);
  line_at_n2->SetLineWidth(1);
  line_at_n2->Draw();

  TLine *line_at_n4 = new TLine(pullFrame->GetXaxis()->GetXmin(), -4, pullFrame->GetXaxis()->GetXmax(), -4);
  line_at_n4->SetLineColor(kBlack);
  line_at_n4->SetLineStyle(3);
  line_at_n4->SetLineWidth(2);
  line_at_n4->Draw();

  TLine *line_at_n8 = new TLine(pullFrame->GetXaxis()->GetXmin(), -8, pullFrame->GetXaxis()->GetXmax(), -8);
  line_at_n8->SetLineColor(kBlack);
  line_at_n8->SetLineStyle(3);
  line_at_n8->SetLineWidth(1);
  line_at_n8->Draw();

  // ===== chi2 ===== //
  printChi2(ws, pullPad, pullTmp, fitResult, "ctau3D", "dataToFit", "TrueModel_Tot", nCtauTrueBins, false);
  pullPad->Update();

  // ===== draw main canvas ===== //
  myCanvas->Draw();
  myCanvas->SaveAs(Form("figs/2DFit_%s/CtauTrue_%s_%s.png", DATE.c_str(), bCont.Data(), kineLabel.c_str()));
  myCanvas->SaveAs(Form("../figs/2DFit_%s/CtauTrue/CtauTrue_%s_%s.pdf", DATE.c_str(), bCont.Data(), kineLabel.c_str()));

  delete pullTmp;
  delete myCanvas;
}

void CtauTrueFit::saveResults()
{
  std::cout << "===== saveResults() =====\n\n";
  // TFile *outputFile = new TFile(Form("roots/2DFit_%s/CtauTrueResult_Inclusive_%s.root", DATE.c_str(), kineLabel.c_str()), "RECREATE");

  // fitResult->Write();
  // ws->pdf("TrueModel_Tot")->Write();

  // outputFile->Close();
  // delete outputFile;
}
