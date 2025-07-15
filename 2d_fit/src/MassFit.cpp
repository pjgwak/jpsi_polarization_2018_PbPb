#include "MassFit.h"
#include <iostream>
#include <algorithm> // std::min
#include <utility> // std::make_pair
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TCanvas.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooMsgService.h"
#include "../headers/polarizationUtilities.h"

using namespace RooFit;

MassFit::MassFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
                 float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP)
    : ptLow(ptLow), ptHigh(ptHigh), yLow(yLow), yHigh(yHigh), cLow(cLow), cHigh(cHigh),
      cosLow(cosLow), cosHigh(cosHigh), PR(PR), PRw(PRw),
      fEffW(fEffW), fAccW(fAccW), isPtW(isPtW), isTnP(isTnP)
{
  gStyle->SetEndErrorSize(0);

  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
}

MassFit::~MassFit()
{
  RooMsgService::instance().getStream(1).addTopic(RooFit::ObjectHandling);

  delete fInputData;
  delete fMcParams;
  delete ws;
  delete fitResult;
}

void MassFit::init()
{
  std::cout << "===== Start init() =====\n\n";
  setLabels();
  openInputFile();
  setupWorkspaceAndData();
  setVariableRanges();
  defineModel();
}

void MassFit::run()
{
  std::cout << "===== Start run() =====\n\n";
  performFit();
  makePlot();
  saveResults();
}
void MassFit::setLabels()
{
  std::cout << "===== Start setLabels() =====\n\n";
  gSystem->mkdir(Form("roots/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("../figs/2DFit_%s/Mass", DATE.c_str()), kTRUE); // summary plots

  kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh).Data();

  fname = (PRw == 1) ? "PR" : "NP";
}

void MassFit::openInputFile()
{
  std::cout << "===== openInputFile() =====\n\n";
  fInputData = new TFile(inputFilePath.c_str());

  if (!fInputData || fInputData->IsZombie())
  {
    std::cerr << "[ERROR]: Failed to open input file: " << inputFilePath << "\n";
    exit(1);
  }
}

void MassFit::setupWorkspaceAndData()
{
  std::cout << "===== Start setupWorkspaceAndData() =====\n\n";
  ws = new RooWorkspace("workspace");
  RooDataSet *dataset = (RooDataSet *)fInputData->Get("dataset");
  if (!dataset)
  {
    std::cerr << "[ERROR] Failed to load RooDataSet 'dataset' from file.\n";
    exit(1);
  }
  ws->import(*dataset);

  TString kineCutStr = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);
  TString accCutStr = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018
  TString osCutStr = "recoQQsign==0";
  TString angleCutStr = Form("cos_ep>%.2f && cos_ep<%.2f", cosLow, cosHigh);
  TString finalCut = TString::Format("%s && %s && %s && %s", osCutStr.Data(), accCutStr.Data(), kineCutStr.Data(), angleCutStr.Data());

  // is weighted?
  RooDataSet *datasetW = nullptr;
  if (isWeighted)
    datasetW = new RooDataSet("datasetW", "A sample", *dataset->get(), Import(*dataset), WeightVar(*ws->var("weight")));
  else
    datasetW = new RooDataSet("datasetW", "A sample", *dataset->get(), Import(*dataset));

  // check variables
  std::vector<std::string> varNames = {"ctau3DRes", "ctau3D", "ctau3DErr", "mass", "pt", "y"};
  RooArgSet selectedVars;

  for (const auto &varName : varNames)
  {
    RooRealVar *var = ws->var(varName.c_str());
    if (!var)
    {
      std::cerr << "[ERROR] Variable '" << varName << "' not found in workspace.\n";
      exit(1);
    }
    selectedVars.add(*var);
  }

  RooDataSet *dsAB = (RooDataSet *)datasetW->reduce(selectedVars, finalCut.Data());
  ws->import(*dsAB, Rename("dsAB"));

  delete datasetW;
}

void MassFit::setVariableRanges()
{
  std::cout << "===== setVariableRanges() =====\n\n";
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->setRange("massFitRange", massLow, massHigh);
  // ws->var("mass")->Print();
}

void MassFit::buildDoubleCB()
{
  std::cout << "===== buildDoubleCB() =====\n\n";
  // --- get MC fit result ---
  // open MC fit file
  TString mcFilePath = Form("roots/2DFit_%s/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root",
                            DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP);
  fMcParams = new TFile(mcFilePath.Data());
  if (!fMcParams || fMcParams->IsZombie())
  {
    std::cerr << "[ERROR] Can't open MC fit result file: " << mcFilePath.Data() << "\n";
    exit(1);
  }

  // load MC dataset to get the fit result
  RooDataSet *dataset_fit = (RooDataSet *)fMcParams->Get("datasetMass");
  RooWorkspace *ws_fit = new RooWorkspace("ws_mc_fit");
  ws_fit->import(*dataset_fit);

  // extract fit results
  double sigma_index = 5; // +- 5 sigma

  // use lambda function
  auto getValErr = [&](const char *name, double &val, double &err)
  {
    val = ws_fit->var(name)->getVal();
    err = ws_fit->var(name)->getError();
    double lower = val - sigma_index * err;
    double upper = val + sigma_index * err;
    if (std::string(name)=="f") {
      if (upper > 1.0) upper = 1.0;
      if (lower < 0.0) lower = 0.0;
    }
    return std::make_pair(lower, upper);
  };

  // store val, and error
  double alpha_val, alpha_err, n_val, n_err, xA_val, xA_err, f_val, f_err, sigma_val, sigma_err;
  auto alpha_range = getValErr("alpha_1_A", alpha_val, alpha_err);
  auto n_range = getValErr("n_1_A", n_val, n_err);
  auto xA_range = getValErr("x_A", xA_val, xA_err);
  auto f_range = getValErr("f", f_val, f_err);
  auto sigma_range = getValErr("sigma_1_A", sigma_val, sigma_err);

  // --- buld mass fit model ---
  // signla parameters
  // mean and sigma is always free parameter - user can change lower, upper bounds in config file
  RooRealVar mean("mean", "mean of the signal", pdgMass.JPsi, pdgMass.JPsi - 0.05, pdgMass.JPsi + 0.05);
  RooRealVar sigma_1_A("sigma_1_A", "width of the signal", sigma_val, sigma_range.first, sigma_range.second);

  // fixed parameters
  RooRealVar x_A("x_A", "sigma ratio", xA_val);
  RooRealVar alpha_1_A("alpha_1_A", "tail shift", alpha_val);
  RooRealVar n_1_A("n_1_A", "power order", n_val);
  RooRealVar f("f", "cb fraction", f_val);

  // for legacy
  // RooRealVar x_A("x_A", "sigma ratio", xA_val, xA_range.first, xA_range.second);
  // RooRealVar alpha_1_A("alpha_1_A", "tail shift", alpha_val, alpha_range.first, alpha_range.second);
  // RooRealVar n_1_A("n_1_A", "power order", n_val, n_range.first, n_range.second);
  // RooRealVar f("f", "cb fraction", f_val, f_range.first, f_range.second);

  RooFormulaVar sigma_2_A("sigma_2_A", "@0*@1", RooArgList(sigma_1_A, x_A));
  RooFormulaVar alpha_2_A("alpha_2_A", "1.0*@0", RooArgList(alpha_1_A));
  RooFormulaVar n_2_A("n_2_A", "1.0*@0", RooArgList(n_1_A));

  // signal models
  RooCBShape cb_1_A("cball_1_A", "Crystal Ball 1", *(ws->var("mass")), mean, sigma_1_A, alpha_1_A, n_1_A);
  RooCBShape cb_2_A("cball_2_A", "Crystal Ball 2", *(ws->var("mass")), mean, sigma_2_A, alpha_2_A, n_2_A);
  ws->import(cb_1_A);
  ws->import(cb_2_A);

  RooAddPdf pdfMASS_Jpsi("pdfMASS_Jpsi", "Signal Shape", RooArgList(cb_1_A, cb_2_A), RooArgList(f));
  ws->import(pdfMASS_Jpsi, RecycleConflictNodes());
}

void MassFit::buildCBG()
{
  std::cout << "===== buildCBG() =====\n\n";
  // --- get MC fit result ---
  // open MC fit file
  TString mcFilePath = Form("roots/2DFit_%s/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root",
                            DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP);
  fMcParams = new TFile(mcFilePath.Data());
  if (!fMcParams || fMcParams->IsZombie())
  {
    std::cerr << "[ERROR] Can't open MC fit result file: " << mcFilePath.Data() << "\n";
    exit(1);
  }

  // load MC dataset to get the fit result
  RooDataSet *dataset_fit = (RooDataSet *)fMcParams->Get("datasetMass");
  RooWorkspace *ws_fit = new RooWorkspace("ws_mc_fit");
  ws_fit->import(*dataset_fit);

  // extract fit results
  double sigma_index = 5; // +- 5 sigma

  // use lambda function
  auto getValErr = [&](const char *name, double &val, double &err)
  {
    val = ws_fit->var(name)->getVal();
    err = ws_fit->var(name)->getError();
    double lower = val - sigma_index * err;
    double upper = val + sigma_index * err;
    if (std::string(name) == "f")
    {
      if (upper > 1.0)
        upper = 1.0;
      if (lower < 0.0)
        lower = 0.0;
    }
    return std::make_pair(lower, upper);
  };

  // store val, and error
  double alpha_val, alpha_err, n_val, n_err, xA_val, xA_err, f_val, f_err, sigma_val, sigma_err;
  auto alpha_range = getValErr("alpha_cb", alpha_val, alpha_err);
  auto n_range = getValErr("n_cb", n_val, n_err);
  auto xA_range = getValErr("x_A", xA_val, xA_err);
  auto f_range = getValErr("f", f_val, f_err);
  auto sigma_range = getValErr("sigma_cb", sigma_val, sigma_err);

  // --- buld mass fit model ---
  // signal parameters
  RooRealVar mean("mean", "cb mean", pdgMass.JPsi, pdgMass.JPsi - 0.05, pdgMass.JPsi + 0.05);
  RooRealVar sigma_cb("sigma_cb", "CB sigma", sigma_val, sigma_range.first, sigma_range.second);

  // fixed parameters
  RooRealVar x_A("x_A", "sigma ratio", xA_val);
  RooRealVar alpha_cb("alpha_cb", "tail shift", alpha_val);
  RooRealVar n_cb("n_cb", "power order", n_val);
  RooRealVar f("f", "cb fraction", f_val);

  RooFormulaVar sigma_gauss("sigma_gauss", "sigma_cb * x_A", "@0*@1", RooArgList(sigma_cb, x_A));

  RooGaussian gauss("gauss", "gaussian PDF", *(ws->var("mass")), mean, sigma_gauss);
  RooCBShape cb("cb", "crystal ball", *(ws->var("mass")), mean, sigma_cb, alpha_cb, n_cb);
  
  // To call individual componets
  ws->import(gauss);
  ws->import(cb);

  RooAddPdf pdfMASS_Jpsi("pdfMASS_Jpsi", "Signal Shape", RooArgList(gauss, cb), RooArgList(f));

  ws->import(pdfMASS_Jpsi, RecycleConflictNodes());
}


void MassFit::buildExpo()
{
  std::cout << "===== buildExpo() =====\n\n";
  // bkg parameter
  RooRealVar lambda("lambda", "mass bkg expo lambda", -0.01, -1., 0.); // expo decay constant

  // bkg model
  RooExponential pdfMASS_bkg("pdfMASS_bkg", "Background", *(ws->var("mass")), lambda);
  ws->import(pdfMASS_bkg);
}

void MassFit::buildCheby(int order)
{
  std::cout << "===== buildCheby() =====\n\n";
  // bkg parameters
  RooRealVar sl1("sl1", "mass bkg sl1", 0.01, -1., 1.);
  RooRealVar sl2("sl2", "mass bkg sl2", 0.01, -1., 1.);
  RooRealVar sl3("sl3", "mass bkg sl3", 0.01, -1., 1.);
  RooRealVar sl4("sl4", "mass bkg sl4", 0.01, -1., 1.);
  RooRealVar sl5("sl5", "mass bkg sl5", 0.01, -1., 1.);
  RooRealVar sl6("sl6", "mass bkg sl6", 0.01, -1., 1.);

  RooArgList coeffs;
  coeffs.add(sl1);
  if (order >= 2) coeffs.add(sl2);
  if (order >= 3) coeffs.add(sl3);
  if (order >= 4) coeffs.add(sl4);
  if (order >= 5) coeffs.add(sl5);
  if (order == 6) coeffs.add(sl6);

  // bkg models
  RooChebychev pdfMASS_bkg("pdfMASS_bkg", "Background", *(ws->var("mass")), coeffs);
  ws->import(pdfMASS_bkg);
}

void MassFit::defineModel()
{
  std::cout << "===== defineModel() =====\n\n";
  // signal
  if (pdfTypeSig == "doubleCB")
    buildDoubleCB();
  else if (pdfTypeSig == "CBG")
    buildCBG();
  else
  {
    std::cout << "[ERROR] failed to build sig model: " << pdfTypeSig << "\n";
    std::cout << "possible PDFs: doulbeCB, CBG\n";
    exit(1);
  }

  // bkg
  if (pdfTypeBkg == "expo")
    buildExpo();
  else if (pdfTypeBkg == "cheby1")
    buildCheby(1);
  else if (pdfTypeBkg == "cheby2")
    buildCheby(2);
  else if (pdfTypeBkg == "cheby3")
    buildCheby(3);
  else if (pdfTypeBkg == "cheby4")
    buildCheby(4);
  else if (pdfTypeBkg == "cheby5")
    buildCheby(5);
  else if (pdfTypeBkg == "cheby6")
    buildCheby(6);
  else
  {
    std::cout << "[ERROR] failed to build bkg model: " << pdfTypeBkg << "\n";
    std::cout << "possible PDFs: expo, cheby1, cheby2 ... cheby6\n";
    exit(1);
  }

  // sig + bkg model
  RooRealVar N_Jpsi("N_Jpsi", "", 100000, 1, 1e7);
  RooRealVar N_Bkg("N_Bkg", "", 100000, 1, 1e7);

  RooAddPdf pdfMASS_Tot("pdfMASS_Tot", "Jpsi + Bkg", RooArgList(*ws->pdf("pdfMASS_Jpsi"), *ws->pdf("pdfMASS_bkg")), RooArgList(N_Jpsi, N_Bkg));

  ws->import(pdfMASS_Tot); // , RecycleConflictNodes()
}

void MassFit::initVar(const std::string &varName, double init, double low, double high)
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

void MassFit::performFit()
{
  std::cout << "===== performFit() =====\n\n";
  RooAbsPdf *pdf = ws->pdf("pdfMASS_Tot");
  RooAbsData *data = ws->data("dsAB");

  bool isWeighted = ws->data("dsAB")->isWeighted();

  fitResult = pdf->fitTo(*data, Save(), Hesse(kTRUE), Range("massFitRange"), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(nCPU), Strategy(2));

  fitResult->Print("V");

  std::cout << "fixing floating parameters to constant (after fitting)\n";
  const RooArgList &floatted_params = fitResult->floatParsFinal();
  for (auto arg : floatted_params)
  {
    RooRealVar *param = dynamic_cast<RooRealVar *>(arg);
    if (param)
    {
      std::cout << "Fixing parameter: " << param->GetName() << "\n";
      ws->var(param->GetName())->setConstant(kTRUE);
    }
  }

  // Todo: turn off setConstantn - N_Jpsi, N_Bkg
}

void MassFit::makePlot()
{
  std::cout << "===== makePlot() =====\n\n";

  TCanvas *myCanvas = new TCanvas("myCanvas", "", 800, 800);
  myCanvas->cd();

  TPad *pad1 = new TPad("pad1", "", 0, 0.25, 1.0, 1.0);
  pad1->SetTicks(1, 1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  gPad->SetLogy();

  RooPlot *massFrame = ws->var("mass")->frame(nMassBin); // bins
  massFrame->SetTitle("");

  ws->data("dsAB")->plotOn(massFrame, Name("dataOS"), MarkerSize(.8));
  ws->pdf("pdfMASS_Tot")->plotOn(massFrame, VisualizeError(*fitResult, 1), FillColor(kOrange), Range("massFitRange"), NormRange("massFitRange"));
  ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("pdfMASS_Tot"), LineColor(kBlack), Range("massFitRange"), NormRange("massFitRange"));
  ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("pdfMASS_bkg"), Components(*ws->pdf("pdfMASS_bkg")), LineColor(kBlue), LineStyle(kDashed), Range("massFitRange"), NormRange("massFitRange"));

  if (pdfTypeSig == "doubleCB")
  {
    ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("bkg_plus_cb1"), Components(RooArgSet(*ws->pdf("pdfMASS_bkg"), *ws->pdf("cball_1_A"))), LineColor(44), LineStyle(kDashDotted), Range("massFitRange"), NormRange("massFitRange"));
    ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("bkg_plus_cb2"), Components(RooArgSet(*ws->pdf("pdfMASS_bkg"), *ws->pdf("cball_2_A"))), LineColor(8), LineStyle(kDashDotted), Range("massFitRange"), NormRange("massFitRange"));
  }
  else if (pdfTypeSig == "CBG")
  {
    ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("bkg_plus_cb"), Components(RooArgSet(*ws->pdf("pdfMASS_bkg"), *ws->pdf("cb"))), LineColor(44), LineStyle(kDashDotted), Range("massFitRange"), NormRange("massFitRange"));
    ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("bkg_plus_gauss"), Components(RooArgSet(*ws->pdf("pdfMASS_bkg"), *ws->pdf("gauss"))), LineColor(8), LineStyle(kDashDotted), Range("massFitRange"), NormRange("massFitRange"));
  }




  massFrame->SetFillStyle(4000);
  massFrame->GetYaxis()->SetTitleOffset(1.43);
  // massFrame->GetYaxis()->CenterTitle();
  // massFrame->GetYaxis()->SetTitleSize(0.058);
  // massFrame->GetYaxis()->SetLabelSize(0.054);
  // massFrame->GetYaxis()->SetRangeUser(ws->var("N_Jpsi")->getVal()/100, ws->var("N_Jpsi")->getVal());

  TH1 *hTmp = ws->data("dsAB")->createHistogram("hist", *ws->var("mass"), Binning(massFrame->GetNbinsX(), massFrame->GetXaxis()->GetXmin(), massFrame->GetXaxis()->GetXmax()));
  Double_t YMax = hTmp->GetBinContent(hTmp->GetMaximumBin());
  // Double_t YMin = min( hTmp->GetBinContent(hTmp->FindFirstBinAbove(0.0)), hTmp->GetBinContent(hTmp->FindLastBinAbove(0.0)) );
  Double_t YMin = 1e99;
  for (int i = 1; i <= hTmp->GetNbinsX(); i++)
    if (hTmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, hTmp->GetBinContent(i));
  double Ydown;
  double Yup;
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
  Yup = 30 * YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));
  massFrame->GetYaxis()->SetRangeUser(Ydown, Yup);
  // massFrame->SetMinimum(2*10);
  massFrame->GetXaxis()->SetLabelSize(0);
  massFrame->GetXaxis()->SetTitleSize(0);
  massFrame->GetXaxis()->CenterTitle();
  massFrame->GetXaxis()->SetRangeUser(massLow, massHigh);
  massFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massFrame->Draw();

  // --- draw data, pdfs legends ---
  TLegend *leg_B = new TLegend(text_x + 0.53, text_y - 0.2, text_x + 0.73, text_y);
  leg_B->SetTextSize(text_size);
  leg_B->SetTextFont(43);
  leg_B->SetBorderSize(0);

  // Need to check
  leg_B->AddEntry(massFrame->findObject("dataOS"), "Data", "pe");
  leg_B->AddEntry(massFrame->findObject("pdfMASS_Tot"), "Total", "l");
  leg_B->AddEntry(massFrame->findObject("pdfMASS_bkg"), "Background", "l");
  if (pdfTypeSig == "doubleCB")
  {
    leg_B->AddEntry(massFrame->findObject("bkg_plus_cb1"), "Bkg + CB1", "l");
    leg_B->AddEntry(massFrame->findObject("bkg_plus_cb2"), "Bkg + CB2", "l");
  }
  else if (pdfTypeSig == "CBG")
  {
    leg_B->AddEntry(massFrame->findObject("bkg_plus_cb"), "Bkg + CB", "l");
    leg_B->AddEntry(massFrame->findObject("bkg_plus_gauss"), "Bkg + gauss", "l");
  }
  leg_B->Draw("same");

  // --- print fit parameters ---
  if (yLow == 0)
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c, |y^{#mu#mu}| < %.1f, Cent. %d - %d%s ", ptLow, ptHigh, yHigh, cLow / 2, cHigh / 2, "%"), text_x, text_y, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c; %.1f < |y^{#mu#mu}| < %.1f; Cent. %d - %d%s", ptLow, ptHigh, yLow, yHigh, cLow / 2, cHigh / 2, "%"), text_x, text_y, text_color, text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f,  N_{Bkg} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError(), ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError()), text_x, text_y - y_diff * 1, text_color, text_size);
  drawText(Form("mean = %.4f #pm %.4f", ws->var("mean")->getVal(), ws->var("mean")->getError()), text_x, text_y - y_diff * 2, text_color, text_size);

  int y_count = 3; // After: y, N_Jpsi, N_Bkg

  if (pdfTypeSig == "doubleCB")
  {
    drawText(Form("#alpha_{J/#psi} = %.4f (fixed)", ws->var("alpha_1_A")->getVal()), text_x, text_y - y_diff * y_count++, text_color, text_size);
    drawText(Form("f_{J/#psi} = %.4f (fixed)", ws->var("f")->getVal()), text_x, text_y - y_diff * y_count++, text_color, text_size);
    drawText(Form("n_{J/#psi} = %.4f (fixed)", ws->var("n_1_A")->getVal()), text_x, text_y - y_diff * y_count++, text_color, text_size);
    drawText(Form("#sigma1_{J/#psi} = %.2f #pm %.2f MeV/c^{2}, (#sigma2/#sigma1)_{J/#psi} = %.3f (fixed)", (ws->var("sigma_1_A")->getVal()) * 1000, (ws->var("sigma_1_A")->getError()) * 1000, ws->var("x_A")->getVal()), text_x, text_y - y_diff * y_count++, text_color, text_size);
  }
  else if (pdfTypeSig == "CBG")
  {
    drawText(Form("#alpha_{J/#psi} = %.4f (fixed)", ws->var("alpha_cb")->getVal()), text_x, text_y - y_diff * y_count++, text_color, text_size);
    drawText(Form("f_{J/#psi} = %.4f (fixed)", ws->var("f")->getVal()), text_x, text_y - y_diff * y_count++, text_color, text_size);
    drawText(Form("n_{J/#psi} = %.4f (fixed)", ws->var("n_cb")->getVal()), text_x, text_y - y_diff * y_count++, text_color, text_size);
    drawText(Form("#sigma1_{J/#psi} = %.2f #pm %.2f MeV/c^{2}, (#sigma2/#sigma1)_{J/#psi} = %.3f (fixed)", (ws->var("sigma_cb")->getVal()) * 1000, (ws->var("sigma_cb")->getError()) * 1000, ws->var("x_A")->getVal()), text_x, text_y - y_diff * y_count++, text_color, text_size);
  }

  // bkgs
  if (pdfTypeBkg == "expo")
  {
    drawText(Form("#lambda_{Expo} = %.4f #pm %.4f", ws->var("lambda")->getVal(), ws->var("lambda")->getError()),
             text_x, text_y - y_diff * y_count++, text_color, text_size);
  }
  else if (pdfTypeBkg.find("cheby") == 0)
  {
    for (int i = 1; i <= 6; ++i)
    {
      std::string parName = Form("sl%d", i);
      if (ws->var(parName.c_str()))
      {
        drawText(Form("sl%d = %.4f #pm %.4f", i, ws->var(parName.c_str())->getVal(), ws->var(parName.c_str())->getError()),
                 text_x, text_y - y_diff * y_count++, text_color, text_size);
      }
    }
  }

  TPad *pad2 = new TPad("pad2", "", 0, 0.0, 1.0, 0.25);
  myCanvas->cd();
  pad2->Draw();
  pad2->cd();
  pad2->SetTopMargin(0); // Upper and lower plot are joined
  pad2->SetBottomMargin(0.67);
  pad2->SetBottomMargin(0.4);
  pad2->SetFillStyle(4000);
  pad2->SetFrameFillStyle(4000);
  pad2->SetTicks(1, 1);
  RooPlot *frameTMP = (RooPlot *)massFrame->Clone("frameTMP");
  RooHist *pullHist = frameTMP->pullHist("dataOS", "pdfMASS_Tot", true);
  pullHist->SetMarkerSize(0.8);
  RooPlot *pullFrame = ws->var("mass")->frame(Title("Pull Distribution"));
  pullFrame->addPlotable(pullHist, "P");
  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.08);
  pullFrame->GetYaxis()->SetLabelSize(0.08);
  pullFrame->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame->GetXaxis()->SetTitleOffset(1.05);
  pullFrame->GetXaxis()->SetLabelOffset(0.04);
  pullFrame->GetXaxis()->SetLabelSize(0.08);
  pullFrame->GetXaxis()->SetTitleSize(0.08);
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw();
  

  TLine *l1 = new TLine(massLow, 0, massHigh, 0);
  l1->SetLineStyle(1);
  l1->Draw("same");

  printChi2(ws, pad2, frameTMP, fitResult, "mass", "dataOS", "pdfMASS_Tot", nMassBin, false);

  Double_t theNLL = fitResult->minNll();
  std::cout << " *** NLL : " << std::setprecision(15) << theNLL << "\n";

  
  myCanvas->Update();
  myCanvas->SaveAs(Form("../figs/2DFit_%s/Mass/Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  myCanvas->SaveAs(Form("figs/2DFit_%s/Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  delete myCanvas;
  delete massFrame;
  delete pullFrame;
  delete hTmp;
  
}

void MassFit::saveResults()
{
  std::cout << "===== saveResults() =====\n\n";
  TFile *outputFile = new TFile(Form("roots/2DFit_%s/MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP), "recreate");
  if (!outputFile || outputFile->IsZombie())
  {
    std::cerr << "[ERROR] Could not create output file.\n";
    return;
  }

  ws->pdf("pdfMASS_Tot")->Write();
  ws->pdf("pdfMASS_bkg")->Write();
  fitResult->Write();

  RooArgSet *fitargs = new RooArgSet();
  fitargs->add(fitResult->floatParsFinal());
  RooDataSet *datasetMass = new RooDataSet("datasetMass", "dataset with Mass Fit result", *fitargs);
  datasetMass->add(*fitargs);
  datasetMass->Write();
  
  outputFile->Close();

  delete fitargs;
  delete datasetMass;
  delete outputFile;
}