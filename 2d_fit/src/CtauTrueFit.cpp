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

CtauTrueFit::CtauTrueFit(float ptLow, float ptHigh,
                       float yLow, float yHigh,
                       int cLow, int cHigh,
                       float cosLow, float cosHigh,
                       int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP)
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

void CtauTrueFit::run()
{
  std::cout << "===== Start run() =====\n\n";
  setLabels();
  openInputFile();
  createKinematicCut();
  setupWorksapceAndData();
  setVariableRanges();
  defineModel();
  performFit();
  makePlot();
  saveResults();
}

void CtauTrueFit::setLabels()
{
  std::cout << "===== Start setLabels() =====\n\n";

  DATE = "No_Weight_2";
  gSystem->mkdir(Form("roots/2DFit_%s/CtauTrue", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/CtauTrue", DATE.c_str()), kTRUE);

  TString kineTmp = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
  kineLabel = std::string(kineTmp.Data());
}


void CtauTrueFit::openInputFile()
{
  std::cout << "===== Start openInputFile() =====\n\n";
  // Todo: user custom input path
  inputFile = new TFile("../files_roodata/RooDataSet_miniAOD_isMC1_NP_Jpsi_GenOnly_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250307.root");

  if (!inputFile || inputFile->IsZombie())
  {
    std::cerr << "CRITICAL: Input file could not be opened. Aborting.\n";
    exit(1);
  }
}

void CtauTrueFit::createKinematicCut()
{
  std::cout << "===== Start createKinematicCut() =====\n\n";
  TString cut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5", ptLow, ptHigh, yLow, yHigh);

  // It doesn't use the accCut
  // TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
  // cut += accCut;

  kineCutMC = std::string(cut.Data());
}

void CtauTrueFit::setupWorksapceAndData()
{
  std::cout << "===== Start setupWorksapceAndData() =====\n\n";
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
  std::cout << "===== Start setVariableRanges() =====\n\n";
  RooRealVar *ctauVar = ws->var("ctau3D");

  ctauVar->setRange(ctau3DMin, ctau3DMax);
  ctauVar->setRange("truePlotRange", ctau3DMin, ctauHigh);
  ctauVar->setRange("trueFitRange", ctau3DMin, ctau3DMax);
}

void CtauTrueFit::defineModel()
{
  std::cout << "===== Start defineModel() =====\n\n";
  // resolution model
  RooTruthModel delta_fcn("RooTruthModel", "", *ws->var("ctau3D"));

  // decay parameters
  RooRealVar lambdaDSS("lambdaDSS", "", 0.1, 0.01, 5);
  RooRealVar lambdaDSS2("lambdaDSS2", "", 0.1, 0.01, 5);
  RooRealVar fDSS("fDSS", "fracion of decay1", 0.6, 0.01, 1);

  // decay models
  RooDecay pdfCTAUTRUEDSS1("pdfCTAUTRUEDSS1", "", *ws->var("ctau3D"), lambdaDSS, delta_fcn, RooDecay::SingleSided);
  RooDecay pdfCTAUTRUEDSS2("pdfCTAUTRUEDSS2", "", *ws->var("ctau3D"), lambdaDSS2, delta_fcn, RooDecay::SingleSided);

  // combine decay models
  RooAddPdf pdfCTAUTRUE("pdfCTAUTRUE", "", RooArgList(pdfCTAUTRUEDSS1, pdfCTAUTRUEDSS2), RooArgList(fDSS));

  // extend decay model for fitting
  RooRealVar N_Jpsi_MC("N_Jpsi_MC", "", 500000, 100000, 1000000);
  auto ctauTrueModel = new RooExtendPdf("TrueModel_Tot", "", pdfCTAUTRUE, N_Jpsi_MC);
  ws->import(*ctauTrueModel);

  delete ctauTrueModel;
}

void CtauTrueFit::performFit()
{
  std::cout << "===== Start performFit() =====\n\n";
  RooDataSet *dataToFit = (RooDataSet *)ws->data("reducedDS_MC")->reduce(Form("ctau3D>=%.f&&ctau3D<%.f", ctau3DMin, ctau3DMax));
  ws->import(*dataToFit, Rename("dataToFit"));

  // bool isWeighted = ws->data("dataToFit")->isWeighted();
  fitResult = ws->pdf("TrueModel_Tot")->fitTo(*ws->data("dataToFit"), Save(), SumW2Error(false), Extended(kTRUE), NumCPU(nCPU), Strategy(2), Range("trueFitRange"), PrintLevel(-1));
  // Range("trueFitRange") makes fit be better

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

void CtauTrueFit::makePlot()
{
  std::cout << "===== Start makePlot() =====\n\n";
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
  ws->data("dataToFit")->plotOn(ctauFrame, Name("dataToFit")); // DataError(RooAbsData::SumW2)
  ws->pdf("TrueModel_Tot")->plotOn(ctauFrame, Range("truePlotRange"), NormRange("trueFitRange"), Name("TrueModel_Tot"));
  if (nExp >= 2)
  {
    ws->pdf("TrueModel_Tot")->plotOn(ctauFrame, Name("comp1"), Components(*ws->pdf("pdfCTAUTRUEDSS1")), LineStyle(kDashed), LineColor(kGreen + 2), Range("truePlotRange"), NormRange("trueFitRange"));
    ws->pdf("TrueModel_Tot")->plotOn(ctauFrame, Name("comp2"), Components(*ws->pdf("pdfCTAUTRUEDSS2")), LineStyle(kDashed), LineColor(kAzure - 4), Range("truePlotRange"), NormRange("trueFitRange"));
  }
  if (nExp >= 3)
  {
    ws->pdf("TrueModel_Tot")->plotOn(ctauFrame, Name("comp3"), Components(*ws->pdf("pdfCTAUTRUEDSS3")), LineStyle(kDotted), LineColor(kOrange + 1), Range("truePlotRange"), NormRange("trueFitRange"));
  }

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
  TLegend *leg = new TLegend(0.6, 0.65, 0.93, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(ctauFrame->findObject("dataToFit"), "Gen ctau3D", "pe");
  leg->AddEntry(ctauFrame->findObject("TrueModel_Tot"), "Total fit", "l");
  if (nExp >= 2)
  {
    leg->AddEntry(ctauFrame->findObject("comp1"), "Decay Exp1", "l");
    leg->AddEntry(ctauFrame->findObject("comp2"), "Decay Exp2", "l");
  }
  if (nExp >= 3)
  {
    leg->AddEntry(ctauFrame->findObject("comp3"), "Decay Exp3", "l");
  }
  leg->Draw();

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
  myCanvas->SaveAs(Form("figs/2DFit_%s/CtauTrue/ctauTrue_%s_%s.png", DATE.c_str(), bCont.Data(), kineLabel.c_str()));

  delete pullTmp;
  delete myCanvas;
}

void CtauTrueFit::saveResults()
{
  std::cout << "===== Start saveResults() =====\n\n";

  TFile *outputFile = new TFile(Form("roots/2DFit_%s/CtauTrue/CtauTrueResult_Inclusive_%s.root", DATE.c_str(), kineLabel.c_str()), "RECREATE");

  fitResult->Write();
  ws->pdf("TrueModel_Tot")->Write();

  outputFile->Close();
  delete outputFile;
}
