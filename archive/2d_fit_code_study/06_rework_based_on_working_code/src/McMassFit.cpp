#include "McMassFit.h"
#include <iostream>
#include <algorithm> // min()
#include "TSystem.h" // gSystem - to control OS
#include "TStyle.h" // gStyle - to control ROOT cosmetic
#include "TFile.h"
#include "TH1D.h"
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


McMassFit::McMassFit(float ptLow, float ptHigh,
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

  // ignore import success
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
  // RooMsgService::instance().getStream(1).removeTopic(RooFit::InputArguments);

  // ignore data handling
  // RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);

  // ignore empty bin warning of residual pull. - comment them out when you want to test the codes.
  // RooMsgService::instance().getStream(0).removeTopic(Plotting);
  // RooMsgService::instance().getStream(1).removeTopic(Plotting);
}

McMassFit::~McMassFit()
{
  RooMsgService::instance().getStream(1).addTopic(RooFit::ObjectHandling);
  delete ws;
}

void McMassFit::run()
{
  std::cout << "===== Start run() =====\n";
  setLabels();
  reduceDataset();
  buildModel();
  drawResult();
  saveResult();
}

void McMassFit::setLabels()
{
  std::cout << "===== Start setLabels() =====\n";
  DATE = "2DFit_No_Weight_2";
  gSystem->mkdir(Form("roots/%s/mc_Mass", DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/%s/mc_Mass", DATE.Data()), kTRUE);

  TString bCont; // not used
  if (PR == 0)
    bCont = "Prompt";
  else if (PR == 1)
    bCont = "NonPrompt";
  else if (PR == 2)
    bCont = "Inclusive";

  TString fname;
  if (PRw == 1)
    fname = "PR";
  else if (PRw == 2)
    fname = "NP";
}

void McMassFit::reduceDataset()
{
  std::cout << "===== Start reduceDataset() =====\n";

  TFile *f1 = new TFile("../../../files_roodata/RooDataSet_miniAOD_isMC1_PR_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root");
  // massLow, massHigh -> defiend in the common header

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
  // RooDataSet *datasetW = new RooDataSet("datasetW", "A sample", *dataset->get(), Import(*dataset), WeightVar(*ws->var("weight")));
  RooDataSet *datasetW = new RooDataSet("datasetW", "A sample", *dataset->get(), Import(*dataset));
  RooDataSet *dsAB = (RooDataSet *)datasetW->reduce(RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data());
  cout << "******** New Combined Dataset ***********" << endl;
  dsAB->SetName("dsAB");
  ws->import(*dsAB);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->setRange("fitRange", massLow, fit_limit);
  ws->var("mass")->Print();
}

void McMassFit::buildModel()
{
  std::cout << "===== Start buildModel() =====\n";
  // basic parameters
  components["mean"] = std::make_unique<RooRealVar>("m_{J/#Psi}", "mean...", 3.096, 3.086, 3.106);
  components["sigma_1_A"] = std::make_unique<RooRealVar>("sigma_1_A", "width...", 0.01, 0.001, 1);
  components["x_A"] = std::make_unique<RooRealVar>("x_A", "sigma ratio", 0.1, 0.03, 5);
  components["alpha_1_A"] = std::make_unique<RooRealVar>("alpha_1_A", "tail shift", 2, 0.1, 5);
  components["n_1_A"] = std::make_unique<RooRealVar>("n_1_A", "power order", 1.5, 0.1, 5);
  components["f"] = std::make_unique<RooRealVar>("f", "cb fraction", 0.6, 0.01, 0.99);

  // FormulaVar
  auto *sigma_1_A_ptr = dynamic_cast<RooRealVar *>(components["sigma_1_A"].get());
  auto *x_A_ptr = dynamic_cast<RooRealVar *>(components["x_A"].get());
  components["sigma_2_A"] = std::make_unique<RooFormulaVar>("sigma_2_A", "@0*@1", RooArgList(*sigma_1_A_ptr, *x_A_ptr)); // 버그 수정!

  auto *alpha_1_A_ptr = dynamic_cast<RooRealVar *>(components["alpha_1_A"].get());
  components["alpha_2_A"] = std::make_unique<RooFormulaVar>("alpha_2_A", "1.0*@0", RooArgList(*alpha_1_A_ptr)); // 버그 수정!

  auto *n_1_A_ptr = dynamic_cast<RooRealVar *>(components["n_1_A"].get());
  components["n_2_A"] = std::make_unique<RooFormulaVar>("n_2_A", "1.0*@0", RooArgList(*n_1_A_ptr)); // 버그 수정!

  // signal PDFs
  auto *mean_ptr = dynamic_cast<RooRealVar *>(components["mean"].get());
  auto *sigma_2_A_ptr = dynamic_cast<RooFormulaVar *>(components["sigma_2_A"].get());
  auto *alpha_2_A_ptr = dynamic_cast<RooFormulaVar *>(components["alpha_2_A"].get());
  auto *n_2_A_ptr = dynamic_cast<RooFormulaVar *>(components["n_2_A"].get());

  components["cb_1_A"] = std::make_unique<RooCBShape>("cball_1_A", "Crystal Ball 1", *(ws->var("mass")), *mean_ptr, *sigma_1_A_ptr, *alpha_1_A_ptr, *n_1_A_ptr);
  components["cb_2_A"] = std::make_unique<RooCBShape>("cball_2_A", "Crystal Ball 2", *(ws->var("mass")), *mean_ptr, *sigma_2_A_ptr, *alpha_2_A_ptr, *n_2_A_ptr);

  // total pdf
  auto *cb_1_A_ptr = dynamic_cast<RooCBShape *>(components["cb_1_A"].get());
  auto *cb_2_A_ptr = dynamic_cast<RooCBShape *>(components["cb_2_A"].get());
  auto *f_ptr = dynamic_cast<RooRealVar *>(components["f"].get());
  components["pdfMASS_Jpsi"] = std::make_unique<RooAddPdf>("pdfMASS_Jpsi", "Signal", RooArgList(*cb_1_A_ptr, *cb_2_A_ptr), RooArgList(*f_ptr));

  components["N_Jpsi"] = std::make_unique<RooRealVar>("N_Jpsi", "inclusive Jpsi signals", 1600000, 1, 2000000);

  auto *pdfMASS_Jpsi_ptr = dynamic_cast<RooAddPdf *>(components["pdfMASS_Jpsi"].get());
  auto *N_Jpsi_ptr = dynamic_cast<RooRealVar *>(components["N_Jpsi"].get());
  components["pdfMASS_Tot"] = std::make_unique<RooAddPdf>("pdfMASS_Tot", "Jpsi + Bkg", RooArgList(*pdfMASS_Jpsi_ptr), RooArgList(*N_Jpsi_ptr));

  ws->import(*components["pdfMASS_Tot"].get());
}

void McMassFit::drawResult()
{
  std::cout << "===== Start drawResult() =====\n";
  // nMassBin = 25;
  TCanvas *c_A = new TCanvas("canvas_A", "My plots", 4, 4, 550, 520);
  c_A->cd();
  TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.16, 0.98, 1.0);
  pad_A_1->SetTicks(1, 1);
  pad_A_1->Draw();
  pad_A_1->cd();
  RooPlot *myPlot_A = ws->var("mass")->frame(nMassBin); // bins
  myPlot_A->SetTitle("");
  ws->data("dsAB")->plotOn(myPlot_A, Name("dataHist_A"));

  pad_A_1->cd();
  bool logY_flag = true;
  if (logY_flag == true)
  {
    gPad->SetLogy();
  }
  RooPlot *myPlot2_A = (RooPlot *)myPlot_A->Clone();
  ws->data("dsAB")->plotOn(myPlot2_A, Name("dataOS"), MarkerSize(.8));
  // bool isWeighted = ws->data("dsAB")->isWeighted();
  bool isWeighted = false;
  cout << endl
       << "********* Starting Mass Dist. Fit **************" << endl
       << endl;
  fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*ws->data("dsAB"), Save(), Hesse(kTRUE), Range(massLow, fit_limit), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(24), Strategy(2));
  cout << endl
       << "********* Finished Mass Dist. Fit **************" << endl
       << endl;
  fitMass->Print("V");

  // Check and get fitted parameters
  // cout << fitN_Jpsi.getVal() << endl;

  /*
  for ( int i = 0; i < fitParams.getSize(); ++i)
  {
    auto & fitPar = (RooRealVar &) fitParams[i];
    std::cout << fitPar.GetName() << " " << fitPar.getVal() << " +- " << fitPar.getError() << std::endl;
  }
  */

  // double f_factor = (double)fitFraction.getVal();
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A, Name("pdfMASS_tot"), LineColor(kBlack), Range(2.6, fit_limit), NormRange("fitRange"));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A, Name("cball_1_A"), Components("cball_1_A"), LineWidth(1), LineColor(kGreen + 2), Range(2.6, fit_limit), NormRange("fitRange"), LineStyle(kDashed));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A, Name("cball_2_A"), Components("cball_2_A"), LineWidth(1), LineColor(kBlue + 2), Range(2.6, fit_limit), NormRange("fitRange"), LineStyle(kDashed));

  // make a pretty plot
  myPlot2_A->SetFillStyle(4000);
  myPlot2_A->GetYaxis()->SetTitleOffset(1.43);
  // myPlot2_A->GetYaxis()->CenterTitle();
  // myPlot2_A->GetYaxis()->SetTitleSize(0.058);
  // myPlot2_A->GetYaxis()->SetLabelSize(0.054);
  // myPlot2_A->GetYaxis()->SetRangeUser(ws->var("N_Jpsi")->getVal()/100, ws->var("N_Jpsi")->getVal());
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

  myPlot2_A->GetXaxis()->SetLabelSize(0);
  myPlot2_A->GetXaxis()->SetTitleSize(0);
  myPlot2_A->GetXaxis()->CenterTitle();
  // if (logY_flag == true) {
  //     myPlot2_A->SetMinimum(2*10);
  //     myPlot2_A->SetMaximum(10000000000);
  // }
  // else {
  //     myPlot2_A->GetYaxis()->SetRangeUser(0,280000);
  // }
  myPlot2_A->GetYaxis()->SetRangeUser(Ydown, Yup);

  myPlot2_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  myPlot2_A->Draw();

  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff, text_color, text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()), text_x, text_y - y_diff * 2, text_color, text_size);
  // drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
  drawText(Form("#alpha = %.4f #pm %.4f", ws->var("alpha_1_A")->getVal(), ws->var("alpha_1_A")->getError()), text_x, text_y - y_diff * 3, text_color, text_size);
  drawText(Form("f = %.4f #pm %.4f", ws->var("f")->getVal(), ws->var("f")->getError()), text_x, text_y - y_diff * 4, text_color, text_size);
  drawText(Form("n_{1} = %.4f #pm %.4f", ws->var("n_1_A")->getVal(), ws->var("n_1_A")->getError()), text_x, text_y - y_diff * 5, text_color, text_size);
  drawText(Form("#sigma_{1} = %.4f #pm %.4f", ws->var("sigma_1_A")->getVal(), ws->var("sigma_1_A")->getError()), text_x, text_y - y_diff * 6, text_color, text_size);
  drawText(Form("#sigma_{2} / #sigma_{1} = %.4f #pm %.4f", ws->var("x_A")->getVal(), ws->var("x_A")->getError()), text_x, text_y - y_diff * 7, text_color, text_size);

  TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.006, 0.98, 0.227);
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
  RooHist *hpull_A = frameTMP->pullHist("dataOS", "pdfMASS_tot", true);
  hpull_A->SetMarkerSize(0.8);
  RooPlot *pullFrame_A = ws->var("mass")->frame(Title("Pull Distribution"));
  pullFrame_A->addPlotable(hpull_A, "P");
  pullFrame_A->SetTitle("");
  pullFrame_A->SetTitleSize(0);
  pullFrame_A->GetYaxis()->SetTitleOffset(0.3);
  pullFrame_A->GetYaxis()->SetTitle("Pull");
  pullFrame_A->GetYaxis()->SetTitleSize(0.15);
  pullFrame_A->GetYaxis()->SetLabelSize(0.15);
  pullFrame_A->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame_A->GetYaxis()->CenterTitle();

  pullFrame_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame_A->GetXaxis()->SetTitleOffset(1.05);
  pullFrame_A->GetXaxis()->SetLabelOffset(0.04);
  pullFrame_A->GetXaxis()->SetLabelSize(0.15);
  pullFrame_A->GetXaxis()->SetTitleSize(0.15);
  pullFrame_A->GetXaxis()->CenterTitle();

  pullFrame_A->GetYaxis()->SetTickSize(0.04);
  pullFrame_A->GetYaxis()->SetNdivisions(404);
  pullFrame_A->GetXaxis()->SetTickSize(0.03);
  pullFrame_A->Draw();

  TLine *l1 = new TLine(massLow, 0, massHigh, 0);
  l1->SetLineStyle(1);
  l1->Draw("same");

  printChi2(ws, pad_A_2, frameTMP, fitMass, "mass", "dataOS", "pdfMASS_tot", nMassBin, false);

  TH1D *outh = new TH1D("fitResults", "fit result", 20, 0, 20);
  outh->GetXaxis()->SetBinLabel(1, "Jpsi");

  float temp1 = ws->var("N_Jpsi")->getVal();
  float temp1err = ws->var("N_Jpsi")->getError();

  outh->SetBinContent(1, temp1);
  outh->SetBinError(1, temp1err);

  fitMass->Print();
  Double_t theNLL = fitMass->minNll();
  cout << " *** NLL : " << theNLL << endl;
  // RooRealVar *Par1 = ws->var("sl1");
  // RooRealVar *Par2 = ws->var("sl2");
  // cout << "Chebychev Par1 : " << Par1->getVal() << " +/- " << Par1->getError() << endl;
  // cout << "Chebychev Par2 : " << Par2->getVal() << " +/- " << Par2->getError() << endl;

  c_A->Update();
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
  c_A->SaveAs(Form("figs/%s/mc_Mass/mc_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
}

void McMassFit::saveResult()
{
  std::cout << "===== start saveResult() =====\n";
  TFile *outfile;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
  outfile = new TFile(Form("roots/%s/mc_Mass/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP), "recreate");
 
  // add fit result
  fitMass->Print("v");

  // add new workspace
  ws->SetName("old_ws_true");

  RooWorkspace newWs("ws_mass", "");
  newWs.import(*ws->pdf("pdfMASS_Tot"));
  newWs.import(*ws->data("dsAB"), Rename("ds_mass"));
  newWs.Write();

  outfile->Close();
}