#include "MassFit.h"
#include <iostream>
#include <algorithm>
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

void MassFit::run()
{
  std::cout << "===== run() =====\n\n";
  setLabels();
  openInputFile();
  setupWorkspaceAndData();
  setVariableRanges();
  initializeParameters();
  defineModel();
  performFit();
  makePlot();
  saveResults();
}
void MassFit::setLabels()
{
  std::cout << "===== setLabels() =====\n\n";
  DATE = "No_Weight_2";
  gSystem->mkdir(Form("roots/2DFit_%s/Mass", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/Mass", DATE.c_str()), kTRUE);

  if (PRw == 1)
    fname = "PR";
  else
    fname = "NP";

  kineLabel = std::string(getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh).Data());
}

void MassFit::openInputFile()
{
  std::cout << "===== openInputFile() =====\n\n";
  fInputData = new TFile("../files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root");

  if (!fInputData || fInputData->IsZombie())
  {
    std::cerr << "CRITICAL: Input data file could not be opened. Aborting.\n";
    exit(1);
  }
}

void MassFit::setupWorkspaceAndData()
{
  std::cout << "===== setupWorkspaceAndData() =====\n\n";
  ws = new RooWorkspace("workspace");
  RooDataSet *dataset = (RooDataSet *)fInputData->Get("dataset");
  ws->import(*dataset);

  TString kineCutStr = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);
  TString accCutStr = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018
  TString osCutStr = "recoQQsign==0";
  TString angleCutStr = Form("cos_ep>%.2f && cos_ep<%.2f", cosLow, cosHigh);
  TString finalCut = TString::Format("%s && %s && %s && %s", osCutStr.Data(), accCutStr.Data(), kineCutStr.Data(), angleCutStr.Data());

  RooDataSet *datasetW = new RooDataSet("datasetW", "", *dataset->get(), Import(*dataset), WeightVar(*ws->var("weight")));
  RooDataSet *dsAB = (RooDataSet *)datasetW->reduce(RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), finalCut.Data());
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


void MassFit::initializeParameters()
{
  std::cout << "===== initializeParameters() =====\n\n";
  fMcParams = new TFile(Form("roots/2DFit_%s/mc_Mass/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root",
                             DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  if (!fMcParams || fMcParams->IsZombie())
  {
    std::cerr << "CRITICAL: MC parameters file could not be opened. Aborting.\n";
    exit(1);
  }

  RooDataSet *dataset_fit = (RooDataSet *)fMcParams->Get("datasetMass");
  RooWorkspace *ws_fit = new RooWorkspace("workspace_fit");
  ws_fit->import(*dataset_fit);
  
  double alpha_MC_value = ws_fit->var("alpha_1_A")->getVal();
  double alpha_MC_value_err = ws_fit->var("alpha_1_A")->getError();
  double n_MC_value = ws_fit->var("n_1_A")->getVal();
  double n_MC_value_err = ws_fit->var("n_1_A")->getError();
  double xA_MC_value = ws_fit->var("x_A")->getVal();
  double xA_MC_value_err = ws_fit->var("x_A")->getError();
  double f_MC_value = ws_fit->var("f")->getVal();
  double f_MC_value_err = ws_fit->var("f")->getError();
  double sigma_MC_value = ws_fit->var("sigma_1_A")->getVal();
  double sigma_MC_value_err = ws_fit->var("sigma_1_A")->getError();

  // legacy
  double sigma_index = 5;
  double alpha_lower = alpha_MC_value - (sigma_index * alpha_MC_value_err);
  double alpha_higher = alpha_MC_value + (sigma_index * alpha_MC_value_err);
  double xA_lower = xA_MC_value - (sigma_index * xA_MC_value_err);
  double xA_higher = xA_MC_value + (sigma_index * xA_MC_value_err);
  double n_lower = n_MC_value - (sigma_index * n_MC_value_err);
  double n_higher = n_MC_value + (sigma_index * n_MC_value_err);
  double f_lower = f_MC_value - (sigma_index * f_MC_value_err);
  double f_higher = f_MC_value + (sigma_index * f_MC_value_err);
  double sigma_lower = sigma_MC_value - (sigma_index * sigma_MC_value_err);
  double sigma_higher = sigma_MC_value + (sigma_index * sigma_MC_value_err);

  if (f_higher > 1.0)
    f_higher = 1.0;
  if (f_lower < 0.0)
    f_lower = 0.0;

  mc_alpha = alpha_MC_value;
  mc_n = n_MC_value;
  mc_sigma = sigma_MC_value;
  mc_x = xA_MC_value;
  mc_f = f_MC_value;

  delete ws_fit;
} 

void MassFit::defineModel()
{
  std::cout << "===== defineModel() =====\n\n";

  // signla parameters
  RooRealVar mean("m_{J/#Psi}", "mean of the signal", pdgMass.JPsi, pdgMass.JPsi - 0.05, pdgMass.JPsi + 0.05);
  RooRealVar sigma_1_A("sigma_1_A", "width of the signal", mc_sigma, mc_sigma - 0.01, mc_sigma + 0.01);
  
  RooRealVar x_A("x_A", "sigma ratio", mc_x);
  RooRealVar alpha_1_A("alpha_1_A", "tail shift", mc_alpha);
  RooRealVar n_1_A("n_1_A", "power order", mc_n);
  RooRealVar f("f", "cb fraction", mc_f);

  RooFormulaVar sigma_2_A("sigma_2_A", "@0*@1", RooArgList(sigma_1_A, x_A));
  RooFormulaVar alpha_2_A("alpha_2_A", "1.0*@0", RooArgList(alpha_1_A));
  RooFormulaVar n_2_A("n_2_A", "1.0*@0", RooArgList(n_1_A));

  // signal models
  RooCBShape cb_1_A("cball_1_A", "Crystal Ball 1", *(ws->var("mass")), mean, sigma_1_A, alpha_1_A, n_1_A);
  RooCBShape cb_2_A("cball_2_A", "Crystal Ball 2", *(ws->var("mass")), mean, sigma_2_A, alpha_2_A, n_2_A);
  ws->import(cb_1_A);
  ws->import(cb_2_A);

  RooAddPdf pdfMASS_Jpsi("pdfMASS_Jpsi", "Signal Shape", RooArgList(cb_1_A, cb_2_A), RooArgList(f));

  // bkg parameters
  RooRealVar sl1("sl1", "sl1", 0.01, -1., 1.);
  RooRealVar sl2("sl2", "sl2", 0.01, -1., 1.);
  // RooRealVar sl3("sl2", "sl2", 0.01, -1., 1.);
  // RooRealVar sl4("sl2", "sl2", 0.01, -1., 1.);
  // RooRealVar sl5("sl2", "sl2", 0.01, -1., 1.);

  // bkg models
  RooChebychev *pdfMASS_bkg;
  if (ptLow < 6.5)
    pdfMASS_bkg = new RooChebychev("pdfMASS_bkg", "Background", *(ws->var("mass")), RooArgList(sl1, sl2));
  else
    pdfMASS_bkg = new RooChebychev("pdfMASS_bkg", "Background", *(ws->var("mass")), RooArgList(sl1));
  ws->import(*pdfMASS_bkg);

  // fit model
  RooRealVar N_Jpsi("N_Jpsi", "", 100000, 1, 1e7);
  RooRealVar N_Bkg("N_Bkg", "", 100000, 1, 1e7);

  RooAddPdf pdfMASS_Tot("pdfMASS_Tot", "Jpsi + Bkg", RooArgList(pdfMASS_Jpsi, *ws->pdf("pdfMASS_bkg")), RooArgList(N_Jpsi, N_Bkg));

  ws->import(pdfMASS_Tot, RecycleConflictNodes());

  delete pdfMASS_bkg;
}

void MassFit::performFit()
{
  std::cout << "===== performFit() =====\n\n";
  RooAbsPdf *pdf = ws->pdf("pdfMASS_Tot");
  RooAbsData *data = ws->data("dsAB");

  bool isWeighted = data->isWeighted();

  // bool isWeighted = true; // test
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

  ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("bkg_plus_cb1"), Components(RooArgSet(*ws->pdf("pdfMASS_bkg"), *ws->pdf("cball_1_A"))), LineColor(44), LineStyle(kDashDotted));
  ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("bkg_plus_cb2"), Components(RooArgSet(*ws->pdf("pdfMASS_bkg"), *ws->pdf("cball_2_A"))), LineColor(8), LineStyle(kDashDotted));


  // ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Components(RooArgSet(*pdfMASS_bkg, *cb_1_A)), LineColor(44));
  // ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Components(RooArgSet(*pdfMASS_bkg, *cb_2_A)), LineColor(8));


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
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));
  massFrame->GetYaxis()->SetRangeUser(Ydown, Yup);
  // massFrame->SetMinimum(2*10);
  massFrame->GetXaxis()->SetLabelSize(0);
  massFrame->GetXaxis()->SetTitleSize(0);
  massFrame->GetXaxis()->CenterTitle();
  massFrame->GetXaxis()->SetRangeUser(massLow, massHigh);
  massFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massFrame->Draw();

  TLegend *leg_B = new TLegend(text_x + 0.53, text_y - 0.2, text_x + 0.73, text_y);
  leg_B->SetTextSize(text_size);
  leg_B->SetTextFont(43);
  leg_B->SetBorderSize(0);
  leg_B->AddEntry(massFrame->findObject("dataOS"), "Data", "pe");
  leg_B->AddEntry(massFrame->findObject("pdfMASS_Tot"), "Total", "l");
  leg_B->AddEntry(massFrame->findObject("pdfMASS_bkg"), "Background", "l");
  leg_B->AddEntry(massFrame->findObject("bkg_plus_cb1"), "Bkg + CB1", "l");
  leg_B->AddEntry(massFrame->findObject("bkg_plus_cb2"), "Bkg + CB2 Component", "l");
  leg_B->Draw("same");

  if (yLow == 0)
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c, |y^{#mu#mu}| < %.1f, Cent. %d - %d%s ", ptLow, ptHigh, yHigh, cLow / 2, cHigh / 2, "%"), text_x, text_y, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c; %.1f < |y^{#mu#mu}| < %.1f; Cent. %d - %d%s", ptLow, ptHigh, yLow, yHigh, cLow / 2, cHigh / 2, "%"), text_x, text_y, text_color, text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f,  N_{Bkg} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError(), ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError()), text_x, text_y - y_diff * 1, text_color, text_size);
  drawText(Form("m_{J/#psi} = %.4f #pm %.4f", ws->var("m_{J/#Psi}")->getVal(), ws->var("m_{J/#Psi}")->getError()), text_x, text_y - y_diff * 2, text_color, text_size);
  drawText(Form("#alpha_{J/#psi} = %.4f (fixed)", ws->var("alpha_1_A")->getVal()), text_x, text_y - y_diff * 3, text_color, text_size);
  drawText(Form("f_{J/#psi} = %.4f (fixed)", ws->var("f")->getVal()), text_x, text_y - y_diff * 4, text_color, text_size);
  drawText(Form("n_{J/#psi} = %.4f (fixed)", ws->var("n_1_A")->getVal()), text_x, text_y - y_diff * 5, text_color, text_size);
  drawText(Form("#sigma1_{J/#psi} = %.2f #pm %.2f MeV/c^{2}, (#sigma2/#sigma1)_{J/#psi} = %.3f (fixed)", (ws->var("sigma_1_A")->getVal()) * 1000, (ws->var("sigma_1_A")->getError()) * 1000, ws->var("x_A")->getVal()), text_x, text_y - y_diff * 6, text_color, text_size);
  // drawText(Form("(#sigma2/#sigma1)_{J/#psi} = %.3f (fixed)",ws->var("x_A")->getVal()),text_x,text_y-y_diff*7,text_color,text_size);

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
  myCanvas->SaveAs(Form("figs/2DFit_%s/Mass/Mass_Fixed_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  
  delete myCanvas;
  delete massFrame;
  delete pullFrame;
  delete hTmp;
}

void MassFit::saveResults()
{
  std::cout << "===== saveResults() =====\n\n";
  TFile *outputFile;
  outputFile = new TFile(Form("roots/2DFit_%s/Mass//Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP), "recreate");

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