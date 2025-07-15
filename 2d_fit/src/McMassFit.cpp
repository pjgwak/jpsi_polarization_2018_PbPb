#include "McMassFit.h"
#include <iostream>
#include <algorithm>
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLine.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooMsgService.h"
#include "../headers/polarizationUtilities.h"

using namespace RooFit;

McMassFit::McMassFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
                     float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP)
    : ptLow(ptLow), ptHigh(ptHigh), yLow(yLow), yHigh(yHigh), cLow(cLow), cHigh(cHigh),
      cosLow(cosLow), cosHigh(cosHigh), PR(PR), PRw(PRw),
      fEffW(fEffW), fAccW(fAccW), isPtW(isPtW), isTnP(isTnP)
{
  gStyle->SetEndErrorSize(0);

  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
}

McMassFit::~McMassFit()
{
  RooMsgService::instance().getStream(1).addTopic(RooFit::ObjectHandling);
  
  delete fInputData;  
  delete ws;
  delete fitResult;
}

void McMassFit::init()
{
  std::cout << "===== Start init() =====\n\n";
  setLabels();
  openInputFile();
  setupWorkspaceAndData();
  setVariableRanges();
  defineModel();
}

void McMassFit::run()
{
  std::cout << "===== Start run() =====\n\n";
  
  performFit();
  makePlot();
  saveResults();
}

void McMassFit::setLabels()
{
  std::cout << "===== Start setLabels() =====\n\n";

  // DATE = "No_Weight_2";
  gSystem->mkdir(Form("roots/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("../figs/2DFit_%s/mc_Mass", DATE.c_str()), kTRUE); // summary plots

  // TString bCont;
  // if (PR == 0)
  //   bCont = "Prompt";
  // else if (PR == 1)
  //   bCont = "NonPrompt";
  // else if (PR == 2)
  //   bCont = "Inclusive";

  kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh).Data();

  fname = (PRw == 1) ? "PR" : "NP";
}

void McMassFit::openInputFile()
{
  std::cout << "===== Start openInputFile() =====\n\n";
  fInputData = new TFile(inputFilePath.c_str());

  if (!fInputData || fInputData->IsZombie())
  {
    std::cerr << "[ERROR]: Failed to open input file: " << inputFilePath << "\n";
    exit(1);
  }
}

void McMassFit::setupWorkspaceAndData()
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

  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut
  TString osCut = "recoQQsign==0";
  TString angleCut = Form("cos_ep>%.2f && cos_ep<%.2f", cosLow, cosHigh);
  TString finalCut = TString::Format("%s && %s && %s && %s", osCut.Data(), accCut.Data(), kineCut.Data(), angleCut.Data());

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

void McMassFit::setVariableRanges()
{
  std::cout << "===== Start setVariableRanges() =====\n\n";
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->setRange("mcPlotRange", massLow, massHigh);
  ws->var("mass")->setRange("mcFitRange", massMin, massMax);
}

void McMassFit::defineModel()
{
  std::cout << "===== Start defineModel() =====\n\n";
  if (pdfType == "doubleCB")
    buildDoubleCB();
  else if (pdfType == "CBG")
    buildCBG();
  else {
    std::cout << "[ERROR] failed to build model: " << pdfType << "\n";
    std::cout << "possible PDFs: doulbeCB, CBG\n";
    exit(1);
  }
}

void McMassFit::buildDoubleCB() {
  // signal parameters
  RooRealVar mean("mean", "mean of the signal", 3.096, 3.086, 3.106);
  RooRealVar sigma_1_A("sigma_1_A", "width of the signal", 0.01, 0.001, 0.1);
  RooRealVar x_A("x_A", "sigma ratio", 1.1, 1, 3);
  RooRealVar alpha_1_A("alpha_1_A", "tail shift", 1.5, 0.8, 5);
  RooRealVar n_1_A("n_1_A", "power order", 1.5, 0.8, 5);
  RooRealVar f("f", "fraction of first CB", 0.6, 0.05, 0.95);

  RooFormulaVar sigma_2_A("sigma_2_A", "sigma_1_A * x_A", "@0*@1", RooArgList(sigma_1_A, x_A));
  RooFormulaVar alpha_2_A("alpha_2_A", "alpha_1_A", "1.0*@0", RooArgList(alpha_1_A));
  RooFormulaVar n_2_A("n_2_A", "n_1_A", "1.0*@0", RooArgList(n_1_A));

  // signal models
  RooCBShape cb_1_A("cball_1_A", "Crystal Ball 1", *(ws->var("mass")), mean, sigma_1_A, alpha_1_A, n_1_A);
  RooCBShape cb_2_A("cball_2_A", "Crystal Ball 2", *(ws->var("mass")), mean, sigma_2_A, alpha_2_A, n_2_A);

  RooAddPdf pdfMASS_Jpsi("pdfMASS_Jpsi", "Signal Shape", RooArgList(cb_1_A, cb_2_A), RooArgList(f));

  // fit model
  RooRealVar N_Jpsi("N_Jpsi", "Number of J/psi signals", 20000, 1, 800000);
  RooAddPdf pdfMASS_Tot("pdfMASS_Tot", "Total MC PDF", RooArgList(pdfMASS_Jpsi), RooArgList(N_Jpsi));

  // To call individual componets
  ws->import(cb_1_A);
  ws->import(cb_2_A);

  // RecycleConflictNodes -> Use componetns insdie worksapce to build complicated model
  ws->import(pdfMASS_Tot, RooFit::RecycleConflictNodes());
}

void McMassFit::buildCBG()
{
  // signal parameters
  RooRealVar mean("mean", "cb mean", 3.096, 3.086, 3.106);
  RooRealVar sigma_cb("sigma_cb", "CB sigma", 0.01, 0.001, 0.1);
  RooRealVar x_A("x_A", "sigma ratio", 1.1, 1, 3);
  RooRealVar alpha_cb("alpha_cb", "tail shift", 1.5, 0.8, 5);
  RooRealVar n_cb("n_cb", "power order", 1.5, 0.8, 5);
  RooRealVar f("f", "fraction of first CB", 0.6, 0.05, 0.95);

  RooFormulaVar sigma_gauss("sigma_gauss", "sigma_cb * x_A", "@0*@1", RooArgList(sigma_cb, x_A));

  RooGaussian gauss("gauss", "gaussian PDF", *(ws->var("mass")), mean, sigma_gauss);
  RooCBShape cb("cb", "crystal ball", *(ws->var("mass")), mean, sigma_cb, alpha_cb, n_cb);

  RooAddPdf pdfMASS_Jpsi("pdfMASS_Jpsi", "Signal Shape", RooArgList(gauss, cb), RooArgList(f));

  // fit model
  RooRealVar N_Jpsi("N_Jpsi", "Number of J/psi signals", 20000, 1, 800000);
  RooAddPdf pdfMASS_Tot("pdfMASS_Tot", "Total MC PDF", RooArgList(pdfMASS_Jpsi), RooArgList(N_Jpsi));

  // To call individual componets
  ws->import(gauss);
  ws->import(cb);

  // RecycleConflictNodes -> Use componetns insdie worksapce to build complicated model
  ws->import(pdfMASS_Tot, RooFit::RecycleConflictNodes());
}

void McMassFit::initVar(const std::string &varName, double init, double low, double high)
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

void McMassFit::performFit()
{
  std::cout << "===== Start performFit() =====\n\n";
  RooAbsPdf *pdf = ws->pdf("pdfMASS_Tot");
  RooAbsData *data = ws->data("dsAB");

  bool isWeighted = data->isWeighted();
  fitResult = pdf->fitTo(*data,
                         Save(),
                         Extended(kTRUE),
                         Range("mcFitRange"),
                         NumCPU(nCPU),
                         SumW2Error(isWeighted),
                         Strategy(2),
                         PrintLevel(0));

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
}

void McMassFit::makePlot()
{
  std::cout << "===== Start makePlot() =====\n\n";
  TCanvas *myCanvas = new TCanvas("myCanvas", "", 800, 800);
  myCanvas->cd();

  TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetTicks(1, 1);
  pad1->Draw();
  pad1->cd();
  gPad->SetLogy();

  RooPlot *massFrame = ws->var("mass")->frame(nMassBin); // bins
  massFrame->SetTitle("");

  ws->data("dsAB")->plotOn(massFrame, Name("dataOS"), MarkerSize(.8));
  ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("pdfMASS_Tot"), LineColor(kBlack), Range("mcFitRange"), NormRange("mcFitRange"));

  if (pdfType=="doubleCB") {
    ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("CB1"), Components("cball_1_A"), LineWidth(3), LineColor(kGreen + 2), Range("mcFitRange"), NormRange("mcFitRange"), LineStyle(kDashed));
    ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("CB2"), Components("cball_2_A"), LineWidth(3), LineColor(kBlue + 2), Range("mcFitRange"), NormRange("mcFitRange"), LineStyle(kDashed));
  } else if (pdfType=="CBG") {
    ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("CB"), Components("cb"), LineWidth(3), LineColor(kGreen + 2), Range("mcFitRange"), NormRange("mcFitRange"), LineStyle(kDashed));
    ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("gauss"), Components("gauss"), LineWidth(3), LineColor(kBlue + 2), Range("mcFitRange"), NormRange("mcFitRange"), LineStyle(kDashed));
  }
  
  // legend
  TLegend *legend = new TLegend(0.6, 0.65, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.03);
  legend->AddEntry(massFrame->findObject("dataOS"), "Data", "pe");
  legend->AddEntry(massFrame->findObject("pdfMASS_Tot"), "Total Fit", "l");

  if (pdfType == "doubleCB")
  {
    legend->AddEntry(massFrame->findObject("CB1"), "CB1", "l");
    legend->AddEntry(massFrame->findObject("CB2"), "CB2", "l");
  }
  else if (pdfType == "CBG")
  {
    legend->AddEntry(massFrame->findObject("CB"), "CB", "l");
    legend->AddEntry(massFrame->findObject("gauss"), "gauss", "l");
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
  Yup = 10 * YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));

  massFrame->GetYaxis()->SetRangeUser(Ydown, Yup);
  massFrame->GetXaxis()->SetLabelSize(0);
  massFrame->GetXaxis()->SetTitleSize(0);
  massFrame->GetXaxis()->CenterTitle();
  massFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massFrame->Draw();
  legend->Draw();

  // draw bin information
  // drawText(Form("%.1f < p_{T} < %.1f GeV/c", ptLow, ptHigh), text_x+0.3, text_y, kBlack, text_size);
  // drawText(Form("%.1f < |y| < %.1f", yLow, yHigh), text_x+0.3, text_y - y_diff, kBlack, text_size);
  // drawText(Form("Cent. %d-%d %%", cLow / 2, cHigh / 2), text_x+0.3, text_y - 2 * y_diff, kBlack, text_size);

  // draw fit parameters
  if (yLow == 0)
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c, |y^{#mu#mu}| < %.1f, Cent. %d - %d%s ", ptLow, ptHigh, yHigh, cLow / 2, cHigh / 2, "%"), text_x, text_y, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c; %.1f < |y^{#mu#mu}| < %.1f; Cent. %d - %d%s", ptLow, ptHigh, yLow, yHigh, cLow / 2, cHigh / 2, "%"), text_x, text_y, text_color, text_size);

  int y_count = 1;

  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()), text_x, text_y - y_diff * y_count++, text_color, text_size);
  drawText(Form("mean = %.4f #pm %.4f", ws->var("mean")->getVal(), ws->var("mean")->getError()), text_x, text_y - y_diff * y_count++, text_color, text_size);
  // drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
  if (pdfType == "doubleCB")
  {
    drawText(Form("#alpha = %.4f #pm %.4f", ws->var("alpha_1_A")->getVal(), ws->var("alpha_1_A")->getError()), text_x, text_y - y_diff * y_count++, text_color, text_size);
    drawText(Form("f = %.4f #pm %.4f", ws->var("f")->getVal(), ws->var("f")->getError()), text_x, text_y - y_diff * y_count++, text_color, text_size);
    drawText(Form("n_{1} = %.4f #pm %.4f", ws->var("n_1_A")->getVal(), ws->var("n_1_A")->getError()), text_x, text_y - y_diff * y_count++, text_color, text_size);
    drawText(Form("#sigma_{1} = %.4f #pm %.4f", ws->var("sigma_1_A")->getVal(), ws->var("sigma_1_A")->getError()), text_x, text_y - y_diff * y_count++, text_color, text_size);
    drawText(Form("#sigma_{2} / #sigma_{1} = %.4f #pm %.4f", ws->var("x_A")->getVal(), ws->var("x_A")->getError()), text_x, text_y - y_diff * y_count++, text_color, text_size);
  }
  else if (pdfType == "CBG")
  {
    drawText(Form("#alpha = %.4f #pm %.4f", ws->var("alpha_cb")->getVal(), ws->var("alpha_cb")->getError()), text_x, text_y - y_diff * y_diff * y_count++, text_color, text_size);
    drawText(Form("f = %.4f #pm %.4f", ws->var("f")->getVal(), ws->var("f")->getError()), text_x, text_y - y_diff * y_diff * y_count++, text_color, text_size);
    drawText(Form("n_{cb} = %.4f #pm %.4f", ws->var("n_cb")->getVal(), ws->var("n_cb")->getError()), text_x, text_y - y_diff * y_diff * y_count++, text_color, text_size);
    drawText(Form("#sigma_{cb} = %.4f #pm %.4f", ws->var("sigma_cb")->getVal(), ws->var("sigma_cb")->getError()), text_x, text_y - y_diff * y_diff * y_count++, text_color, text_size);
    drawText(Form("#sigma_{gauss} / #sigma_{cb} = %.4f #pm %.4f", ws->var("x_A")->getVal(), ws->var("x_A")->getError()), text_x, text_y - y_diff * y_diff * y_count++, text_color, text_size);
  }
  

  TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.3);
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
  RooPlot *pullFrame = ws->var("mass")->frame(Title(" "));
  pullFrame->addPlotable(pullHist, "P");
  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.15);
  pullFrame->GetYaxis()->SetLabelSize(0.15);
  pullFrame->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame->GetXaxis()->SetTitleOffset(1.05);
  pullFrame->GetXaxis()->SetLabelOffset(0.04);
  pullFrame->GetXaxis()->SetLabelSize(0.15);
  pullFrame->GetXaxis()->SetTitleSize(0.15);
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
  std::cout << " *** NLL : " << theNLL << "\n";

  myCanvas->Update();
  myCanvas->SaveAs(Form("../figs/2DFit_%s/mc_Mass/mc_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  myCanvas->SaveAs(Form("figs/2DFit_%s/mc_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  delete legend;
  delete hTmp;
  delete myCanvas;
  delete massFrame;
  delete frameTMP;
  delete pullFrame;
  delete l1;
}

void McMassFit::saveResults()
{
  std::cout << "===== Start saveResults() =====\n\n";
  TFile *outputFile = new TFile(Form("roots/2DFit_%s/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP), "recreate");

  ws->pdf("pdfMASS_Tot")->Write();

  RooArgSet *fitargs = new RooArgSet();
  fitargs->add(fitResult->floatParsFinal());
  RooDataSet *datasetMass = new RooDataSet("datasetMass", "dataset with Mass Fit result", *fitargs);
  datasetMass->add(*fitargs);
  datasetMass->Write();

  // fitResult->Print("v");

  outputFile->Close();

  delete fitargs;
  delete datasetMass;
  delete outputFile;
}