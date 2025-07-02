#include "MassFit.h"
#include <iostream>
#include <algorithm> // min()
#include "TSystem.h" // gSystem - to control OS
#include "TStyle.h"  // gStyle - to control ROOT cosmetic
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
#include "../../../headers/polarizationUtilities.h"

using namespace std;
using namespace RooFit;

MassFit::MassFit(float ptLow, float ptHigh,
                     float yLow, float yHigh,
                     int cLow, int cHigh,
                     float cosLow, float cosHigh,
                     int PR,
                     int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP)
    : ptLow(ptLow), ptHigh(ptHigh),
      yLow(yLow), yHigh(yHigh),
      cLow(cLow), cHigh(cHigh),
      cosLow(cosLow), cosHigh(cosHigh),
      PR(PR), PRw(PRw),
      fEffW(fEffW), fAccW(fAccW), isPtW(isPtW), isTnP(isTnP)
{
  gStyle->SetEndErrorSize(0);
  kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

  // ignore import success
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
  // RooMsgService::instance().getStream(1).removeTopic(RooFit::InputArguments);

  // ignore data handling
  // RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);

  // ignore empty bin warning of residual pull. - comment them out when you want to test the codes.
  // RooMsgService::instance().getStream(0).removeTopic(Plotting);
  // RooMsgService::instance().getStream(1).removeTopic(Plotting);
}

MassFit::~MassFit()
{
  RooMsgService::instance().getStream(1).addTopic(RooFit::ObjectHandling);
  delete ws;
  delete fitMass;
}

void MassFit::run()
{
  std::cout << "===== Start run() =====\n";
  setLabels();
  reduceDataset();
  buildModel();
  fitModel();
  drawResult();
  saveResult();
}

void MassFit::setLabels()
{
  std::cout << "===== Start setLabels() =====\n";
  DATE = "No_Weight_2";
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("roots/2DFit_%s/Mass", DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/Mass", DATE.Data()), kTRUE);

  if (PRw == 1)
    fname = "PR";
  else if (PRw == 2)
    fname = "NP";
}

void MassFit::reduceDataset()
{
  std::cout << "===== Start reduceDataset() =====\n";
  TFile *f1 = new TFile("../../../files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root");

  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
  TString OS = "recoQQsign==0 &&";

  TString angle_cut = Form("&& cos_ep>%.2f && cos_ep<%.2f", cosLow, cosHigh);

  kineCut = OS + accCut + kineCut + angle_cut;

  RooDataSet *dataset = (RooDataSet *)f1->Get("dataset");
  ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  ws->data("dataset")->Print();
  cout << "pt: " << ptLow << "-" << ptHigh << ", y: " << yLow << "-" << yHigh << ", Cent: " << cLow / 2 << "-" << cHigh / 2 << "%" << endl;
  cout << "####################################" << endl;
  // Todo: user custom weight flag
  RooDataSet *datasetW = new RooDataSet("datasetW", "A sample", *dataset->get(), Import(*dataset), WeightVar(*ws->var("weight")));
  // RooDataSet *datasetW = new RooDataSet("datasetW","A sample",*dataset->get(),Import(*dataset));
  RooDataSet *dsAB = (RooDataSet *)datasetW->reduce(RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data());
  cout << "******** New Combined Dataset ***********" << endl;
  dsAB->SetName("dsAB");
  // ws->import(*dsAB);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();
  ws->import(*dsAB);

  delete dsAB;
  delete datasetW;
  delete f1;
}

void MassFit::buildModel()
{
  std::cout << "===== Start buildModel() =====\n";

  // ===== get parameters from MC =====
  TFile *f_fit = new TFile(Form("roots/2DFit_%s/mc_Mass/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));

  auto ws_mc = (RooWorkspace *)f_fit->Get("ws_mc");
  Double_t alpha_MC_value = ws_mc->var("alpha_1_A")->getVal();
  // Double_t alpha_MC_value_err = ws_mc->var("alpha_1_A")->getError();
  Double_t n_MC_value = ws_mc->var("n_1_A")->getVal();
  // Double_t n_MC_value_err = ws_mc->var("n_1_A")->getError();
  Double_t xA_MC_value = ws_mc->var("x_A")->getVal();
  // Double_t xA_MC_value_err = ws_mc->var("x_A")->getError();
  Double_t f_MC_value = ws_mc->var("f")->getVal();
  Double_t f_MC_value_err = ws_mc->var("f")->getError();
  Double_t sigma_MC_value = ws_mc->var("sigma_1_A")->getVal();
  // Double_t sigma_MC_value_err = ws_mc->var("sigma_1_A")->getError();

  double sigma_index = 5;

  // Double_t alpha_lower = alpha_MC_value - (sigma_index * alpha_MC_value_err);
  // Double_t alpha_higher = alpha_MC_value + (sigma_index * alpha_MC_value_err);
  // Double_t xA_lower = xA_MC_value - (sigma_index * xA_MC_value_err);
  // Double_t xA_higher = xA_MC_value + (sigma_index * xA_MC_value_err);
  // Double_t n_lower = n_MC_value - (sigma_index * n_MC_value_err);
  // Double_t n_higher = n_MC_value + (sigma_index * n_MC_value_err);
  Double_t f_lower = f_MC_value - (sigma_index * f_MC_value_err);
  Double_t f_higher = f_MC_value + (sigma_index * f_MC_value_err);


  // Double_t sigma_lower = sigma_MC_value - (sigma_index * sigma_MC_value_err);
  // Double_t sigma_higher = sigma_MC_value + (sigma_index * sigma_MC_value_err);

  // if (n_lower<0.0)n_lower==0.0;
  if (f_higher > 1.0)
    f_higher = 1.0;
  if (f_lower < 0.0)
    f_lower = 0.0;

  // ===== build mass model =====
  // signla parameters
  components["mean"] = std::make_unique<RooRealVar>("m_{J/#Psi}", "mean of the signal", pdgMass.JPsi, pdgMass.JPsi - 0.05, pdgMass.JPsi + 0.05);
  components["x_A"] = std::make_unique<RooRealVar>("x_A", "sigma ratio", xA_MC_value);
  components["sigma_1_A"] = std::make_unique<RooRealVar>("sigma_1_A", "", sigma_MC_value, sigma_MC_value - 0.01, sigma_MC_value + 0.01);
  components["alpha_1_A"] = std::make_unique<RooRealVar>("alpha_1_A", "tail shift", alpha_MC_value);
  components["n_1_A"] = std::make_unique<RooRealVar>("n_1_A", "power order", n_MC_value);
  components["f"] = std::make_unique<RooRealVar>("f", "cb fraction", f_MC_value);

  // RooFormular variables
  auto *sigma_1_A_ptr = dynamic_cast<RooRealVar *>(components["sigma_1_A"].get());
  auto *x_A_ptr = dynamic_cast<RooRealVar *>(components["x_A"].get());
  components["sigma_2_A"] = std::make_unique<RooFormulaVar>("sigma_2_A", "@0*@1", RooArgList(*sigma_1_A_ptr, *x_A_ptr));

  auto *alpha_1_A_ptr = dynamic_cast<RooRealVar *>(components["alpha_1_A"].get());
  components["alpha_2_A"] = std::make_unique<RooFormulaVar>("alpha_2_A", "1.0*@0", RooArgList(*alpha_1_A_ptr));

  auto *n_1_A_ptr = dynamic_cast<RooRealVar *>(components["n_1_A"].get());
  components["n_2_A"] = std::make_unique<RooFormulaVar>("n_2_A", "1.0*@0", RooArgList(*n_1_A_ptr));

  // signal models
  auto *mean_ptr = dynamic_cast<RooRealVar *>(components["mean"].get());
  auto *sigma_2_A_ptr = dynamic_cast<RooFormulaVar *>(components["sigma_2_A"].get());
  auto *alpha_2_A_ptr = dynamic_cast<RooFormulaVar *>(components["alpha_2_A"].get());
  auto *n_2_A_ptr = dynamic_cast<RooFormulaVar *>(components["n_2_A"].get());
  components["cb_1_A"] = std::make_unique<RooCBShape>("cb_1_A", "Crystal Ball 1", *(ws->var("mass")), *mean_ptr, *sigma_1_A_ptr, *alpha_1_A_ptr, *n_1_A_ptr);
  components["cb_2_A"] = std::make_unique<RooCBShape>("cb_2_A", "Crystal Ball 2", *(ws->var("mass")), *mean_ptr, *sigma_2_A_ptr, *alpha_2_A_ptr, *n_2_A_ptr);

  auto *cb_1_A_ptr = dynamic_cast<RooCBShape *>(components["cb_1_A"].get());
  auto *cb_2_A_ptr = dynamic_cast<RooCBShape *>(components["cb_2_A"].get());
  auto *f_ptr = dynamic_cast<RooRealVar *>(components["f"].get());
  components["pdfMASS_Jpsi"] = std::make_unique<RooAddPdf>("pdfMASS_Jpsi", "Signal", RooArgList(*cb_1_A_ptr, *cb_2_A_ptr), RooArgList(*f_ptr));

  // bkg parameters and model
  components["sl1"] = std::make_unique<RooRealVar>("sl1", "sl1", 0.0, -1., 1.);
  components["sl2"] = std::make_unique<RooRealVar>("sl2", "sl2", 0.0, -1., 1.);

  auto *sl1_ptr = dynamic_cast<RooRealVar *>(components["sl1"].get());
  auto *sl2_ptr = dynamic_cast<RooRealVar *>(components["sl2"].get());
  components["pdfMASS_bkg"] = std::make_unique<RooChebychev>("pdfMASS_bkg", "Background", *(ws->var("mass")), RooArgList(*sl1_ptr, *sl2_ptr));

  // total model
  components["N_Jpsi"] = std::make_unique<RooRealVar>("N_Jpsi", "inclusive Jpsi signals", 100000, 1, 1000000);
  components["N_Bkg"] = std::make_unique<RooRealVar>("N_Bkg", "background signals", 100000, 1, 1000000);

  auto *pdfMASS_Jpsi_ptr = dynamic_cast<RooAddPdf *>(components["pdfMASS_Jpsi"].get());
  auto *pdfMASS_bkg_ptr = dynamic_cast<RooChebychev *>(components["pdfMASS_bkg"].get());
  auto *N_Jpsi_ptr = dynamic_cast<RooRealVar *>(components["N_Jpsi"].get());
  auto *N_Bkg_ptr = dynamic_cast<RooRealVar *>(components["N_Bkg"].get());
  components["pdfMASS_Tot"] = std::make_unique<RooAddPdf>("pdfMASS_Tot", "Jpsi + Bkg", RooArgList(*pdfMASS_Jpsi_ptr, *pdfMASS_bkg_ptr), RooArgList(*N_Jpsi_ptr, *N_Bkg_ptr));

  ws->import(*components["pdfMASS_Tot"].get());

  delete f_fit;
}

void MassFit::fitModel()
{
  std::cout << "===== Start fitModel() =====\n";
  bool isWeighted = ws->data("dsAB")->isWeighted();
  // bool isWeighted = false; // test
  fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*ws->data("dsAB"), Save(), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(nCPU), Strategy(2));
  fitMass->Print("V");

  // fix parameters for next steps
  std::cout << "fixing floating parameters to constant (after fitting)\n";
  const RooArgList &floatted_params = fitMass->floatParsFinal();
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

void MassFit::drawResult()
{
  // std::cout << "===== Start drawResult() =====\n";
  TCanvas *c_A = new TCanvas("canvas_A", "My plots", 4, 4, 550, 520);
  c_A->cd();
  TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.25, 0.98, 1.0);
  pad_A_1->SetTicks(1, 1);
  pad_A_1->Draw();
  pad_A_1->cd();
  RooPlot *myPlot_A = ws->var("mass")->frame(nMassBin); // bins
  myPlot_A->SetTitle("");
  ws->data("dsAB")->plotOn(myPlot_A, Name("dataHist_A"));
  pad_A_1->cd();
  gPad->SetLogy();
  RooPlot *myPlot2_A = (RooPlot *)myPlot_A->Clone();
  ws->data("dsAB")->plotOn(myPlot2_A, Name("dataOS"), MarkerSize(.8));

  auto *bkg_ptr = components["pdfMASS_bkg"].get();
  auto *cb1_ptr = components["cb_1_A"].get();
  auto *cb2_ptr = components["cb_2_A"].get();

  RooArgSet bkg_plus_cb1(*bkg_ptr, *cb1_ptr);
  RooArgSet bkg_plus_cb2(*bkg_ptr, *cb2_ptr);

  // ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A, VisualizeError(*fitMass, 1), FillColor(kOrange));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A, Name("pdfMASS_Tot"), LineColor(kBlack));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A, Components(bkg_plus_cb1), LineColor(44), LineWidth(1), LineStyle(kDashed));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A, Components(bkg_plus_cb2), LineColor(8), LineWidth(1), LineStyle(kDashed));
  // model.plotOn(myPlot2_A,Name("model"), LineColor(kBlack));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A, Name("pdfMASS_bkg"), Components("pdfMASS_bkg"), LineColor(kBlue), LineStyle(kDashed), LineWidth(2));
  ws->data("dsAB")->plotOn(myPlot2_A, Name("dataOS"), MarkerSize(.8));
  // model.plotOn(myPlot2_A,Name("model_Bkg"),Components(RooArgSet(*pdfMASS_bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
  // make a pretty plot
  myPlot2_A->SetFillStyle(4000);
  myPlot2_A->GetYaxis()->SetTitleOffset(1.43);
  TH1 *h = ws->data("dsAB")->createHistogram("hist", *ws->var("mass"), Binning(myPlot_A->GetNbinsX(), myPlot_A->GetXaxis()->GetXmin(), myPlot_A->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  // Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
  Double_t YMin = 1e99;
  for (int i = 1; i <= h->GetNbinsX(); i++)
    if (h->GetBinContent(i) > 0)
      YMin = min(YMin, h->GetBinContent(i));
  double Ydown;
  double Yup;
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));
  myPlot2_A->GetYaxis()->SetRangeUser(Ydown, Yup);
  // myPlot2_A->SetMinimum(2*10);
  myPlot2_A->GetXaxis()->SetLabelSize(0);
  myPlot2_A->GetXaxis()->SetTitleSize(0);
  myPlot2_A->GetXaxis()->CenterTitle();
  myPlot2_A->GetXaxis()->SetRangeUser(massLow, massHigh);
  myPlot2_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  myPlot2_A->Draw();

  TLegend *leg_B = new TLegend(text_x + 0.5, text_y - 0.2, text_x + 0.7, text_y);
  leg_B->SetTextSize(text_size);
  leg_B->SetTextFont(43);
  leg_B->SetBorderSize(0);
  leg_B->AddEntry(myPlot2_A->findObject("dataOS"), "Data", "pe");
  leg_B->AddEntry(myPlot2_A->findObject("pdfMASS_Tot"), "Total", "l");
  leg_B->AddEntry(myPlot2_A->findObject("pdfMASS_bkg"), "Background", "l");
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

  TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.001, 0.98, 0.32);
  c_A->cd();
  pad_A_2->Draw();
  pad_A_2->cd();
  pad_A_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_A_2->SetBottomMargin(0.67);
  pad_A_2->SetBottomMargin(0.4);
  pad_A_2->SetFillStyle(4000);
  pad_A_2->SetFrameFillStyle(4000);
  pad_A_2->SetTicks(1, 1);

  RooPlot *frameTMP = (RooPlot *)myPlot2_A->Clone("TMP");
  RooHist *hpull_A = frameTMP->pullHist("dataOS", "pdfMASS_Tot", true);
  hpull_A->SetMarkerSize(0.8);
  RooPlot *pullFrame_A = ws->var("mass")->frame(Title("Pull Distribution"));
  pullFrame_A->addPlotable(hpull_A, "P");
  pullFrame_A->SetTitle("");
  pullFrame_A->SetTitleSize(0);
  pullFrame_A->GetYaxis()->SetTitleOffset(0.3);
  pullFrame_A->GetYaxis()->SetTitle("Pull");
  pullFrame_A->GetYaxis()->SetTitleSize(0.08);
  pullFrame_A->GetYaxis()->SetLabelSize(0.08);
  pullFrame_A->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame_A->GetYaxis()->CenterTitle();

  pullFrame_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame_A->GetXaxis()->SetTitleOffset(1.05);
  pullFrame_A->GetXaxis()->SetLabelOffset(0.04);
  pullFrame_A->GetXaxis()->SetLabelSize(0.08);
  pullFrame_A->GetXaxis()->SetTitleSize(0.08);
  pullFrame_A->GetXaxis()->CenterTitle();

  pullFrame_A->GetYaxis()->SetTickSize(0.04);
  pullFrame_A->GetYaxis()->SetNdivisions(404);
  pullFrame_A->GetXaxis()->SetTickSize(0.03);
  pullFrame_A->Draw();

  TLine *l1 = new TLine(massLow, 0, massHigh, 0);
  l1->SetLineStyle(1);
  l1->Draw("same");

  printChi2(ws, pad_A_2, frameTMP, fitMass, "mass", "dataOS", "pdfMASS_Tot", nMassBin, false);

  TH1D *outh = new TH1D("fitResults", "fit result", 20, 0, 20);
  outh->GetXaxis()->SetBinLabel(1, "Jpsi");

  float temp1 = ws->var("N_Jpsi")->getVal();
  float temp1err = ws->var("N_Jpsi")->getError();

  outh->SetBinContent(1, temp1);
  outh->SetBinError(1, temp1err);

  fitMass->Print();
  Double_t theNLL = fitMass->minNll();
  cout << " *** NLL : " << std::setprecision(15) << theNLL << endl;

  // ws->Print();
  RooArgSet *fitargs = new RooArgSet();
  fitargs->add(fitMass->floatParsFinal());
  RooDataSet *datasetMass = new RooDataSet("datasetMass", "dataset with Mass Fit result", *fitargs);
  datasetMass->add(*fitargs);

  c_A->Update();
  c_A->SaveAs(Form("figs/2DFit_%s/Mass/Mass_Fixed_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));

  delete datasetMass;
  delete fitargs;
  delete outh;
  delete c_A;
  delete myPlot_A;
  delete frameTMP;
}

void MassFit::saveResult()
{
  std::cout << "===== start saveResult() =====\n";
  TFile *outfile;
  outfile = new TFile(Form("roots/2DFit_%s/Mass//Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP), "recreate");

  fitMass->Write("fitMass");
  // pdfMASS_Tot->Write();
  // pdfMASS_bkg->Write();
  // datasetMass->Write();
  // outh->Write();

  // add new workspace
  ws->SetName("old_ws");

  RooWorkspace newWs("ws_mass", "");
  newWs.import(*ws->pdf("pdfMASS_Tot"));
  // newWs.import(*ws->data("dsAB"), Rename("ds_mass"));
  newWs.Write();

  outfile->Close();
  delete outfile;
}