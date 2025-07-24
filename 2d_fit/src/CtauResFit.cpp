#include "CtauResFit.h"
#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"

#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooGaussModel.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "RooGaussModel.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "../headers/polarizationUtilities.h"

using namespace RooFit;

CtauResFit::CtauResFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
                       float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP)
    : ptLow(ptLow), ptHigh(ptHigh), yLow(yLow), yHigh(yHigh), cLow(cLow), cHigh(cHigh),
      cosLow(cosLow), cosHigh(cosHigh), PR(PR), PRw(PRw), fEffW(fEffW), fAccW(fAccW), isPtW(isPtW), isTnP(isTnP)
{
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
  gStyle->SetEndErrorSize(0);
}

CtauResFit::~CtauResFit()
{
  RooMsgService::instance().getStream(1).addTopic(RooFit::ObjectHandling);
  if (fMass) { fMass->Close(); delete fMass; }
  if (fCErr) { fCErr->Close(); delete fCErr; }
  delete ws;
  delete fitResult;
}

void CtauResFit::init()
{
  std::cout << "===== init() =====\n\n";
  setLabels();
  openInputFile();
  setupWorkspaceAndData();
  setVariableRanges();
  defineModel();
}

void CtauResFit::run()
{
  std::cout << "===== run() =====\n\n";
  performFit();
  makePlot();
  makeRatioPlot();
  saveResults();
}

void CtauResFit::setLabels()
{
  std::cout << "===== setLabels() =====\n\n";
  DATE = "No_Weight_2";
  gSystem->mkdir(Form("roots/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("../figs/2DFit_%s/CtauRes", DATE.c_str()), kTRUE);

  kineLabel = std::string(getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh).Data());

  if (PRw == 1) fname = "PR";
  else if (PRw == 2) fname = "NP";
}

void CtauResFit::openInputFile()
{
  std::cout << "===== openInputFile() =====\n\n";
  fMass = new TFile(Form("roots/2DFit_%s/MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  fCErr = new TFile(Form("roots/2DFit_%s/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

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
}

void CtauResFit::setupWorkspaceAndData()
{
  std::cout << "===== setupWorkspaceAndData() =====\n\n";
  RooAddPdf *pdfMASS_Tot = (RooAddPdf *)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Sig = (RooDataSet *)fCErr->Get("dataw_Sig");

  ws = new RooWorkspace("workspace");
  ws->import(*pdfMASS_Tot);
  ws->import(*dataw_Sig);

  double ctauErrDataMin, ctauErrDataMax;
  dataw_Sig->getRange(*ws->var("ctau3DErr"), ctauErrDataMin, ctauErrDataMax);

  TString cutTmp = Form("ctau3DRes<0 && ctau3DErr>%f && ctau3DErr<%f", ctauErrDataMin, ctauErrDataMax);
  RooDataSet *ctauResCutDS = (RooDataSet *)dataw_Sig->reduce(cutTmp.Data());
  ctauResCutDS->SetName("ctauResCutDS");
  ws->import(*ctauResCutDS);
  delete ctauResCutDS;
}

void CtauResFit::setVariableRanges()
{
  std::cout << "===== setVariableRanges() =====\n\n";
  TH1D *hTot = (TH1D *)ws->data("dataw_Sig")->createHistogram(("hTot"), *ws->var("ctau3DRes"), Binning(nCtauResBins, ctauResLow, ctauResHigh));

  if (ctauResMin == -100) {
    ctauResMin = hTot->GetBinLowEdge(hTot->FindFirstBinAbove(1, 1));
    if (ctauResMin < -10)
      ctauResMin = -10;
  }

  if (ctauResMax == 100) {
    ctauResMax = hTot->GetBinLowEdge(hTot->FindLastBinAbove(0, 1)) + hTot->GetBinWidth(hTot->FindLastBinAbove(0, 1));
    if (ctauResMax > -10)
      ctauResMax = 10;
  }
  delete hTot;

  ws->var("ctau3DRes")->setRange("resFitRange", ctauResMin, ctauResMax);
}

void CtauResFit::defineModel()
{
  std::cout << "===== defineModel() =====\n\n";

  // common constant
  ws->import(RooRealVar("One", "", 1.0));

  if (nGauss == 1) buildCtauRes1Gaus();
  else if (nGauss == 2) buildCtauRes2Gaus();
  else if (nGauss == 3) buildCtauRes3Gaus();
  else if (nGauss == 4) buildCtauRes4Gaus();
  else {
    std::cout << "[ERROR] failed to build a " << nGauss << " gaussian resolution model\n";
    std::cout << "possible nGuass: 1 - 4\n";
    exit(1);
  }
}

void CtauResFit::buildCtauRes1Gaus(){
  std::cout << "===== buildCtauRes1Gaus() =====\n\n";
  RooRealVar ctauMean("ctauRes_mean", "mean of resolution", 0.0);
  RooRealVar sigma1("s1_CtauRes", "sigma of 1st Gaussian", 0.7, 0.1, 1.1);
  RooRealVar ctau1("ctau1_CtauRes", "bias of 1st Gaussian", 0.0);

  // build model
  RooRealVar *ctauVar = ws->var("ctau3DRes");
  RooGaussModel gm1("GaussModel1_ctauRes", "1Gaus Res Model",
                    *ctauVar, ctau1, sigma1, *ws->var("One"), *ws->var("One"));

  // yield
  RooRealVar nJpsi("N_Jpsi", "Number of J/psi", 100000, 1000, 1e6);

  // Final model
  RooExtendPdf total("GaussModel_Tot", "Extended Total Model",
                     gm1, nJpsi);
  ws->import(total);
}

void CtauResFit::buildCtauRes2Gaus(){
  std::cout << "===== buildCtauRes2Gaus() =====\n\n";
  RooRealVar ctauMean("ctauRes_mean", "mean of resolution", 0.0);

  RooRealVar ctau1("ctau1_CtauRes", "bias of 1st Gaussian", 0.0);
  RooRealVar ctau2("ctau2_CtauRes", "bias of 2nd Gaussian", 0.0);

  RooRealVar sigma1("s1_CtauRes", "sigma of 1st Gaussian", 0.7, 0.1, 1.1);
  RooRealVar rS21("rS21_CtauRes", "ratio sigma2/sigma1", 1.5, 1.0, 3.0);
  RooFormulaVar sigma2("s2_CtauRes", "@0*@1", RooArgList(rS21, sigma1));

  RooRealVar f("f_CtauRes", "fraction of 1st Gaussian", 0.5, 0.01, 1.0);


  // build model
  RooRealVar *ctauVar = ws->var("ctau3DRes");

  RooGaussModel gm1("GaussModel1_ctauRes", "1st Gaussian",
                    *ctauVar, ctau1, sigma1, *ws->var("One"), *ws->var("One"));
  RooGaussModel gm2("GaussModel2_ctauRes", "2nd Gaussian",
                    *ctauVar, ctau2, sigma2, *ws->var("One"), *ws->var("One"));

  RooAddPdf sumModel("GaussModelCOND_ctauRes", "nGaus 2",
                     RooArgList(gm1, gm2), RooArgList(f));

  // yield
  RooRealVar nJpsi("N_Jpsi", "Number of J/psi", 100000, 1000, 1e8);

  RooExtendPdf total("GaussModel_Tot", "Extended Total Model",
                     sumModel, nJpsi);
  ws->import(total);
}

void CtauResFit::buildCtauRes3Gaus(){
  std::cout << "===== buildCtauRes3Gaus() =====\n\n";
  // f2*gm2 + gm3 = gm23
  // total = f*gm1 + gm23

  RooRealVar ctauMean("ctauRes_mean", "mean of resolution", 0.0);

  RooRealVar ctau1("ctau1_CtauRes", "bias of 1st Gaussian", 0.0);
  RooRealVar ctau2("ctau2_CtauRes", "bias of 2nd Gaussian", 0.0);
  RooRealVar ctau3("ctau3_CtauRes", "bias of 3rd Gaussian", 0.0);

  RooRealVar sigma1("s1_CtauRes", "sigma of 1st Gaussian", 0.7, 0.1, 1.1);
  RooRealVar rS21("rS21_CtauRes", "ratio sigma2/sigma1", 1.5, 1.0, 3.0); // should be > 1.0
  RooRealVar rS32("rS32_CtauRes", "ratio sigma3/sigma2", 1.5, 1.0, 3.0);

  RooFormulaVar sigma2("s2_CtauRes", "@0*@1", RooArgList(rS21, sigma1));
  RooFormulaVar sigma3("s3_CtauRes", "@0*@1", RooArgList(rS32, sigma2));

  RooRealVar f("f_CtauRes", "fraction of 1st Gaussian", 0.3, 0.01, 1.0);
  RooRealVar f2("f2_CtauRes", "fraction of 2nd Gaussian", 0.3, 0.01, 1.0);

  // build res model
  RooRealVar *ctauVar = ws->var("ctau3DRes");

  RooGaussModel gm1("GaussModel1_ctauRes", "1st Gaussian",
                    *ctauVar, ctau1, sigma1, *ws->var("One"), *ws->var("One"));
  RooGaussModel gm2("GaussModel2_ctauRes", "2nd Gaussian",
                    *ctauVar, ctau2, sigma2, *ws->var("One"), *ws->var("One"));
  RooGaussModel gm3("GaussModel3_ctauRes", "3rd Gaussian",
                    *ctauVar, ctau3, sigma3, *ws->var("One"), *ws->var("One"));

  RooAddPdf gm23("GaussModel23_ctauRes", "f2*G2 + G3",
                 RooArgList(gm3, gm2), RooArgList(f2));

  RooAddPdf sumModel("GaussModelCOND_ctauRes", "f*G1 + (f2*G2 + G3)",
                 RooArgList(gm1, gm23), RooArgList(f));

  // yield
  RooRealVar nJpsi("N_Jpsi", "Number of J/psi", 100000, 1000, 1e8);

  RooExtendPdf total("GaussModel_Tot", "Extended Total Model",
                     sumModel, nJpsi);
  ws->import(total);
}

void CtauResFit::buildCtauRes4Gaus(){
  std::cout << "===== buildCtauRes4Gaus() =====\n\n";
  RooRealVar ctauMean("ctauRes_mean", "mean of resolution", 0.0);

  RooRealVar ctau1("ctau1_CtauRes", "bias of 1st Gaussian", 0.0);
  RooRealVar ctau2("ctau2_CtauRes", "bias of 2nd Gaussian", 0.0);
  RooRealVar ctau3("ctau3_CtauRes", "bias of 3rd Gaussian", 0.0);
  RooRealVar ctau4("ctau4_CtauRes", "bias of 4th Gaussian", 0.0);

  RooRealVar sigma1("s1_CtauRes", "sigma of 1st Gaussian", 0.7, 0.1, 1.1);
  RooRealVar rS21("rS21_CtauRes", "ratio sigma2/sigma1", 1.5, 1.0, 3.0);
  RooRealVar rS32("rS32_CtauRes", "ratio sigma3/sigma2", 1.5, 1.0, 3.0);
  RooRealVar rS43("rS43_CtauRes", "ratio sigma4/sigma3", 1.5, 1.0, 3.0);

  RooFormulaVar sigma2("s2_CtauRes", "@0*@1", RooArgList(rS21, sigma1));
  RooFormulaVar sigma3("s3_CtauRes", "@0*@1", RooArgList(rS32, sigma2));
  RooFormulaVar sigma4("s4_CtauRes", "@0*@1", RooArgList(rS43, sigma3));

  RooRealVar f("f_CtauRes", "fraction of 1st Gaussian", 0.3, 0.01, 1.0);
  RooRealVar f2("f2_CtauRes", "fraction of 2nd Gaussian", 0.3, 0.01, 1.0);
  RooRealVar f3("f3_CtauRes", "fraction of 3rd Gaussian", 0.5, 1e-6, 1.0);

  // build model
  RooRealVar *ctauVar = ws->var("ctau3DRes");

  RooGaussModel gm1("GaussModel1_ctauRes", "1st Gaussian",
                    *ctauVar, ctau1, sigma1, *ws->var("One"), *ws->var("One"));
  RooGaussModel gm2("GaussModel2_ctauRes", "2nd Gaussian",
                    *ctauVar, ctau2, sigma2, *ws->var("One"), *ws->var("One"));
  RooGaussModel gm3("GaussModel3_ctauRes", "3rd Gaussian",
                    *ctauVar, ctau3, sigma3, *ws->var("One"), *ws->var("One"));
  RooGaussModel gm4("GaussModel4_ctauRes", "4th Gaussian",
                    *ctauVar, ctau4, sigma4, *ws->var("One"), *ws->var("One"));

  // (G3 + G4)
  RooAddPdf gm34("GaussModel43_ctauRes", "G4 + G3",
                 RooArgList(gm4, gm3), RooArgList(f3));

  // (G2 + (G3 + G4))
  RooAddPdf gm234("GaussModel23_ctauRes", "G2 + G34",
                  RooArgList(gm34, gm2), RooArgList(f2));

  // (G1 + (G2 + G3 + G4))
  RooAddPdf sumModel("GaussModelCOND_ctauRes", "G1 + G234",
                    RooArgList(gm1, gm234), RooArgList(f));

  // yield
  RooRealVar nJpsi("N_Jpsi", "Number of J/psi", 100000, 1000, 1e8);
  ws->import(nJpsi);

  RooExtendPdf total("GaussModel_Tot", "Extended Total Model",
                     sumModel, nJpsi);
  ws->import(total);
}

void CtauResFit::initVar(const std::string &varName, double init, double low, double high)
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

void CtauResFit::performFit()
{
  std::cout << "===== performFit() =====\n\n";
  RooAbsData *initial_data = ws->data("ctauResCutDS");
  RooDataSet *dataToFit = (RooDataSet *)initial_data->reduce(CutRange("resFitRange"));
  ws->import(*dataToFit, Rename("dataToFit"));
  delete dataToFit;

  RooAbsPdf *model = ws->pdf("GaussModel_Tot");
  RooAbsData *finalData = ws->data("dataToFit"); // don't use the dataToFit local variable.

  bool isWeighted = finalData->isWeighted();
  // isWeighted = false; // test
  fitResult = model->fitTo(*finalData,
                           Save(),
                           Extended(true),
                           NumCPU(nCPU),
                           PrintLevel(-1),
                           SumW2Error(isWeighted),
                           Strategy(2));

  fitResult->Print("V");
}

void CtauResFit::drawTextVar(const char *varName, const char *label, float xp, float yp, int textColor, int textSize)
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

void CtauResFit::drawTextVarInt(const char *varName, const char *label, float xp, float yp, int textColor, int textSize)
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

void CtauResFit::makePlot()
{
  std::cout << "===== makePlot() =====\n\n";
  TCanvas *myCanvas = new TCanvas("canvas_C", "", 800, 800);

  TPad *pad1 = new TPad("pad1", "", 0, 0.16, 0.98, 1.0);
  pad1->SetTicks(1, 1);
  pad1->Draw();
  pad1->cd();
  gPad->SetLogy();

  RooPlot *resFrame = ws->var("ctau3DRes")->frame(Bins(nCtauResBins), Range(ctauResLow, ctauResHigh)); // bins
  resFrame->SetTitle("");

  ws->data("dataToFit")->plotOn(resFrame, Name("dataToFit"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(0)); // must use it to calculate pull and Chi2
  ws->data("ctauResCutDS")->plotOn(resFrame, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), LineColor(kRed + 2), MarkerColor(kRed + 2)); // to plot full ranges

  ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_ctauRes"), NormRange("resFitRange"), LineColor(kBlack));

  // components
  ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm1"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel1_ctauRes")), LineColor(kGreen + 2));
  if (nGauss >= 2)
    ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm2"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel2_ctauRes")), LineColor(kRed + 2));
  if (nGauss >=3)
    ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm3"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel3_ctauRes")), LineColor(kBlue + 2));
  if (nGauss >= 4)
    ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm4"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel4_ctauRes")), LineColor(kMagenta + 2));

  // y-axis range
  TH1 *hTmp = ws->data("ctauResCutDS")->createHistogram("hTmp", *ws->var("ctau3DRes"), Binning(resFrame->GetNbinsX(), resFrame->GetXaxis()->GetXmin(), resFrame->GetXaxis()->GetXmax()));
  double YMax = hTmp->GetBinContent(hTmp->GetMaximumBin());
  double YMin = 1e99;
  for (int i = 1; i <= hTmp->GetNbinsX(); i++)
    if (hTmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, hTmp->GetBinContent(i));
  delete hTmp;

  double Yup(0.), Ydown(0.);
  Yup = YMax * TMath::Power((YMax / YMin), (0.5 / (1.0 - 0.5 - 0.2)));
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.2 / (1.0 - 0.5 - 0.2))));
  resFrame->GetYaxis()->SetRangeUser(Ydown, Yup);

  double outTot = ws->data("ctauResCutDS")->numEntries();
  double outRes = ws->data("dataToFit")->numEntries();

  std::cout << "\n\nTot evt: " << outTot << "\n";
  std::cout << "Res evt: " << outRes << "\n";
  std::cout << "lost evt: " << (outTot - outRes) << " events (" << (outTot - outRes) * 100 / outTot << " %)\n";

  if (outRes > 0.0)
  {
    // TLine   *minline = new TLine(rangeErr[0], 0.0, rangeErr[0], Ydown*TMath::Power((Yup/Ydown),0.4));
    TLine *minline = new TLine(ctauResMin, 0.0, ctauResMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
    minline->SetLineStyle(2);
    minline->SetLineColor(1);
    minline->SetLineWidth(3);
    resFrame->addObject(minline);
    TLine *maxline = new TLine(ctauResMax, 0.0, ctauResMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
    maxline->SetLineStyle(2);
    maxline->SetLineColor(1);
    maxline->SetLineWidth(3);
    resFrame->addObject(maxline);
  }
  resFrame->GetXaxis()->CenterTitle();
  resFrame->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
  resFrame->SetFillStyle(4000);
  resFrame->GetYaxis()->SetTitleOffset(1.43);
  resFrame->GetXaxis()->SetLabelSize(0);
  resFrame->GetXaxis()->SetTitleSize(0);
  resFrame->Draw();

  TLegend *leg_C = new TLegend(text_x + 0.29, text_y + 0.03, text_x + 0.39, text_y - 0.17);
  leg_C->SetTextSize(text_size);
  leg_C->SetTextFont(43);
  leg_C->SetBorderSize(0);
  leg_C->AddEntry(resFrame->findObject("dataHist_ctauRes"), "Data", "pe");
  leg_C->AddEntry(resFrame->findObject("modelHist_ctauRes"), "Total PDF", "l");
  leg_C->AddEntry(resFrame->findObject("modelHist_gm1"), "Gauss 1", "l");
  if (nGauss >= 2)
    leg_C->AddEntry(resFrame->findObject("modelHist_gm2"), "Gauss 2", "l");
  if (nGauss >= 3)
    leg_C->AddEntry(resFrame->findObject("modelHist_gm3"), "Gauss 3", "l");
  if (nGauss >= 4)
    leg_C->AddEntry(resFrame->findObject("modelHist_gm4"), "Gauss 4", "l");
  leg_C->Draw("same");


  // print parameters
  int yCountLeft = 0;
  int yCountRight = 0;

  // kinematics (left side)
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);

  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff*yCountLeft++, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);
  drawText(Form("Loss: %.f evts (%.4f %s)", outTot - outRes, (outTot - outRes) * 100 / outTot, "%"), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);

  // parameters (right side)
  drawTextVarInt("N_Jpsi", "N_{J/#psi}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  drawTextVar("s1_CtauRes", "s1_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);

  if (nGauss >= 2)
  {
    drawTextVar("rS21_CtauRes", "(s2/s1)_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f_CtauRes", "f1_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  }
  if (nGauss >= 3)
  {
    drawTextVar("rS32_CtauRes", "(s3/s2)_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f2_CtauRes", "f2_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  }
  if (nGauss >= 4)
  {
    drawTextVar("rS43_CtauRes", "(s4/s3)_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f3_CtauRes", "f3_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  }

  // pull pad
  TPad *pad2 = new TPad("pad2", "", 0, 0.006, 0.98, 0.227);
  myCanvas->cd();
  pad2->Draw();
  pad2->cd();
  pad2->SetTopMargin(0); // Upper and lower plot are joined
  pad2->SetBottomMargin(0.67);
  pad2->SetBottomMargin(0.4);
  pad2->SetFillStyle(4000);
  pad2->SetFrameFillStyle(4000);
  pad2->SetTicks(1, 1);

  // ignore empty bin warning of residual pull. - comment them out when you want to test the codes.
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);

  RooPlot *frameTMP = (RooPlot *)resFrame->Clone("frameTMP");
  RooHist *pullHist = frameTMP->pullHist("dataToFit", "modelHist_ctauRes", true);
  pullHist->SetMarkerSize(0.8);
  // RooPlot* pullFrame = ws->var("ctau3DRes")->frame(Title("Pull Distribution"), Bins(nCtauResBins), Range(ctauResLow, ctauResHigh)) ;
  RooPlot *pullFrame = ws->var("ctau3DRes")->frame(Title("Pull Distribution"), Bins(resFrame->GetNbinsX()), Range(resFrame->GetXaxis()->GetXmin(), resFrame->GetXaxis()->GetXmax()));
  pullFrame->addPlotable(pullHist, "PX");
  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.15);
  pullFrame->GetYaxis()->SetLabelSize(0.15);
  // pullFrame->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame->GetYaxis()->SetRangeUser(-15, 15);
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
  pullFrame->GetXaxis()->SetTitleOffset(1.05);
  pullFrame->GetXaxis()->SetLabelOffset(0.04);
  pullFrame->GetXaxis()->SetLabelSize(0.15);
  pullFrame->GetXaxis()->SetTitleSize(0.15);
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw();

  TLine *lC = new TLine(ctauResLow, 0, ctauResHigh, 0);
  lC->SetLineStyle(1);
  lC->Draw("same");

  printChi2(ws, pad2, frameTMP, fitResult, "ctau3DRes", "dataToFit", "modelHist_ctauRes", nCtauResBins, false);
  pad2->Update();

  myCanvas->Update();
  myCanvas->SaveAs(Form("figs/2DFit_%s/CtauRes_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  myCanvas->SaveAs(Form("../figs/2DFit_%s/CtauRes/CtauRes_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  // turn on alarm
  RooMsgService::instance().getStream(0).addTopic(Plotting);
  RooMsgService::instance().getStream(1).addTopic(Plotting);

  delete frameTMP;
  delete resFrame;
  delete pullFrame;
  delete myCanvas;
}

void CtauResFit::makeRatioPlot()
{
  std::cout << "===== makeRatioPlot() =====\n";
  TCanvas *myCanvas = new TCanvas("canvas_C", "", 800, 800);

  TPad *pad1 = new TPad("pad1", "", 0, 0.16, 0.98, 1.0);
  pad1->SetTicks(1, 1);
  pad1->Draw();
  pad1->cd();
  gPad->SetLogy();

  RooPlot *resFrame = ws->var("ctau3DRes")->frame(Bins(nCtauResBins), Range(ctauResLow, ctauResHigh)); // bins
  resFrame->SetTitle("");

  ws->data("dataToFit")->plotOn(resFrame, Name("dataToFit"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(0));                                                        // must use it to calculate pull and Chi2
  ws->data("ctauResCutDS")->plotOn(resFrame, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), LineColor(kRed + 2), MarkerColor(kRed + 2)); // to plot full ranges

  ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_ctauRes"), NormRange("resFitRange"), LineColor(kBlack));

  // components
  ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm1"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel1_ctauRes")), LineColor(kGreen + 2));
  if (nGauss >= 2)
    ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm2"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel2_ctauRes")), LineColor(kRed + 2));
  if (nGauss >= 3)
    ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm3"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel3_ctauRes")), LineColor(kBlue + 2));
  if (nGauss >= 4)
    ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm4"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel4_ctauRes")), LineColor(kMagenta + 2));

  // y-axis range
  TH1 *hTmp = ws->data("ctauResCutDS")->createHistogram("hTmp", *ws->var("ctau3DRes"), Binning(resFrame->GetNbinsX(), resFrame->GetXaxis()->GetXmin(), resFrame->GetXaxis()->GetXmax()));
  double YMax = hTmp->GetBinContent(hTmp->GetMaximumBin());
  double YMin = 1e99;
  for (int i = 1; i <= hTmp->GetNbinsX(); i++)
    if (hTmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, hTmp->GetBinContent(i));
  delete hTmp;

  double Yup(0.), Ydown(0.);
  Yup = YMax * TMath::Power((YMax / YMin), (0.5 / (1.0 - 0.5 - 0.2)));
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.2 / (1.0 - 0.5 - 0.2))));
  resFrame->GetYaxis()->SetRangeUser(Ydown, Yup);

  double outTot = ws->data("ctauResCutDS")->numEntries();
  double outRes = ws->data("dataToFit")->numEntries();

  std::cout << "\n\nTot evt: " << outTot << "\n";
  std::cout << "Res evt: " << outRes << "\n";
  std::cout << "lost evt: " << (outTot - outRes) << " events (" << (outTot - outRes) * 100 / outTot << " %)\n";

  if (outRes > 0.0)
  {
    // TLine   *minline = new TLine(rangeErr[0], 0.0, rangeErr[0], Ydown*TMath::Power((Yup/Ydown),0.4));
    TLine *minline = new TLine(ctauResMin, 0.0, ctauResMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
    minline->SetLineStyle(2);
    minline->SetLineColor(1);
    minline->SetLineWidth(3);
    resFrame->addObject(minline);
    TLine *maxline = new TLine(ctauResMax, 0.0, ctauResMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
    maxline->SetLineStyle(2);
    maxline->SetLineColor(1);
    maxline->SetLineWidth(3);
    resFrame->addObject(maxline);
  }
  resFrame->GetXaxis()->CenterTitle();
  resFrame->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
  resFrame->SetFillStyle(4000);
  resFrame->GetYaxis()->SetTitleOffset(1.43);
  resFrame->GetXaxis()->SetLabelSize(0);
  resFrame->GetXaxis()->SetTitleSize(0);
  resFrame->Draw();

  TLegend *leg_C = new TLegend(text_x + 0.29, text_y + 0.03, text_x + 0.39, text_y - 0.17);
  leg_C->SetTextSize(text_size);
  leg_C->SetTextFont(43);
  leg_C->SetBorderSize(0);
  leg_C->AddEntry(resFrame->findObject("dataHist_ctauRes"), "Data", "pe");
  leg_C->AddEntry(resFrame->findObject("modelHist_ctauRes"), "Total PDF", "l");
  leg_C->AddEntry(resFrame->findObject("modelHist_gm1"), "Gauss 1", "l");
  if (nGauss >= 2)
    leg_C->AddEntry(resFrame->findObject("modelHist_gm2"), "Gauss 2", "l");
  if (nGauss >= 3)
    leg_C->AddEntry(resFrame->findObject("modelHist_gm3"), "Gauss 3", "l");
  if (nGauss >= 4)
    leg_C->AddEntry(resFrame->findObject("modelHist_gm4"), "Gauss 4", "l");
  leg_C->Draw("same");

  // print parameters
  int yCountLeft = 0;
  int yCountRight = 0;

  // kinematics (left side)
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);

  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);
  drawText(Form("Loss: %.f evts (%.4f %s)", outTot - outRes, (outTot - outRes) * 100 / outTot, "%"), text_x, text_y - y_diff * yCountLeft++, text_color, text_size);

  // parameters (right side)
  drawTextVarInt("N_Jpsi", "N_{J/#psi}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  drawTextVar("s1_CtauRes", "s1_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);

  if (nGauss >= 2)
  {
    drawTextVar("rS21_CtauRes", "(s2/s1)_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f_CtauRes", "f1_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  }
  if (nGauss >= 3)
  {
    drawTextVar("rS32_CtauRes", "(s3/s2)_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f2_CtauRes", "f2_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  }
  if (nGauss >= 4)
  {
    drawTextVar("rS43_CtauRes", "(s4/s3)_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
    drawTextVar("f3_CtauRes", "f3_{Res}", text_x + 0.5, text_y - y_diff * yCountRight++, text_color, text_size);
  }

  TPad *ratioPad = new TPad("ratioPad", "", 0, 0.006, 0.98, 0.227);
  myCanvas->cd();
  // ratioPad->SetLeftMargin(0.15);
  // ratioPad->SetRightMargin(0.07);
  ratioPad->SetTopMargin(0); // Upper and lower plot are joined
  // ratioPad->SetTopMargin(0.05);
  ratioPad->SetBottomMargin(0.35);
  ratioPad->SetTicks(1, 1);
  ratioPad->SetFillStyle(4000);
  ratioPad->SetFrameFillStyle(4000);
  ratioPad->Draw();
  ratioPad->cd();

  RooHist *h_dataPoints = (RooHist *)resFrame->findObject("dataToFit");
  RooCurve *curveRatio = (RooCurve *)resFrame->findObject("modelHist_ctauRes");

  RooPlot *ratioFrame = ws->var("ctau3DRes")->frame(Bins(nCtauResBins), Range(ctauResLow, ctauResHigh));
  ratioFrame->SetTitle(" ");

  TGraphAsymmErrors *g_ratio = new TGraphAsymmErrors();
  int point_idx = 0;
  double x_min = ws->var("ctau3DRes")->getMin("resFitRange");
  double x_max = ws->var("ctau3DRes")->getMax("resFitRange");
  for (int i = 0; i < h_dataPoints->GetN(); ++i)
  {
    double x, y;
    h_dataPoints->GetPoint(i, x, y);
    double model_val = curveRatio->Eval(x);
    if (x < x_min || x > x_max)
      continue;

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
  ratioFrame->GetYaxis()->SetRangeUser(0, 5); // 0, 3
  ratioFrame->GetYaxis()->CenterTitle();

  ratioFrame->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
  ratioFrame->GetXaxis()->SetTitleOffset(1.4);
  ratioFrame->GetXaxis()->SetLabelOffset(0.04);
  ratioFrame->GetXaxis()->SetLabelSize(0.08);
  ratioFrame->GetXaxis()->SetTitleSize(0.08);
  ratioFrame->GetXaxis()->CenterTitle();

  ratioFrame->GetYaxis()->SetTickSize(0.04);
  ratioFrame->GetYaxis()->SetNdivisions(404);
  ratioFrame->GetXaxis()->SetTickSize(0.03);
  ratioFrame->Draw();

  g_ratio->SetMarkerStyle(kFullCircle);
  g_ratio->SetMarkerSize(0.7);
  g_ratio->Draw("P SAME");

  double x_min_line = ratioFrame->GetXaxis()->GetXmin();
  double x_max_line = ratioFrame->GetXaxis()->GetXmax();
  TLine *line_at_1 = new TLine(x_min_line, 1.0, x_max_line, 1.0);
  line_at_1->SetLineColor(kRed);
  line_at_1->SetLineStyle(kDashed);
  line_at_1->Draw("same");

  myCanvas->Draw();
  myCanvas->SaveAs(Form("figs/2DFit_%s/CtauRes_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_ratio.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  myCanvas->SaveAs(Form("../figs/2DFit_%s/CtauRes/CtauRes_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_ratio.pdf", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  delete myCanvas;
}

void CtauResFit::saveResults()
{
  std::cout << "===== saveResults() =====\n\n";
  TFile *outputFile = new TFile(Form("roots/2DFit_%s/CtauResResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP), "recreate");

  ws->pdf("GaussModel_Tot")->Write("GaussModel_Tot");
  fitResult->Write();

  outputFile->Close();
  delete outputFile;
}