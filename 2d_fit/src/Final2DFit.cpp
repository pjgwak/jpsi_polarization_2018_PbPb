#include "Final2DFit.h"
#include <iostream>

#include "TStyle.h"
#include "RooMsgService.h"
#include "TFile.h"
#include "RooWorkspace.h"
#include "TString.h" // Form()
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooCurve.h"

#include "TSystem.h"
#include "TString.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "TH1D.h"
#include "RooPlot.h"
#include "TF1.h"
#include "RooFormulaVar.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "RooProdPdf.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "RooHist.h"

#include "../headers/polarizationUtilities.h"

using std::cout; using std::endl;
using namespace RooFit;

Final2DFit::Final2DFit(float ptLow, float ptHigh,
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

Final2DFit::~Final2DFit()
{
  RooMsgService::instance().getStream(0).addTopic(Caching);
  RooMsgService::instance().getStream(1).addTopic(Caching);
  RooMsgService::instance().getStream(0).addTopic(Plotting);
  RooMsgService::instance().getStream(1).addTopic(Plotting);
  RooMsgService::instance().getStream(0).addTopic(Integration);
  RooMsgService::instance().getStream(1).addTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);

  delete ws;
  if (f1) { f1->Close(); delete f1; }
  if (fMass) { fMass->Close(); delete fMass; }
  if (fCErr) { fCErr->Close(); delete fCErr; }
  if (fCRes) { fCRes->Close(); delete fCRes; }
  if (fCBkg) { fCBkg->Close(); delete fCBkg; }
  if (fCTrue) { fCTrue->Close(); delete fCTrue; }
  delete fitMass;
  delete fit2DFit;
}

void Final2DFit::run()
{
  std::cout << "===== Start run() =====\n\n";
  setLabels();
  createKinematicCut();
  openInputFiles();
  setupWorksapceAndData();
  setVariableRanges();
  defineAndFitMassModel();
  defineCtauModel();
  define2DFitModel();
  perform2DFit();
  makePlots();
  saveResults();
}

void Final2DFit::setLabels()
{
  std::cout << "===== Start setLabels() =====\n\n";
  DATE = "No_Weight_2";
  gSystem->mkdir(Form("roots/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("../figs/2DFit_%s/Final", DATE.c_str()), kTRUE);

  if (PRw == 1)
    fname = "PR";
  else if (PRw == 2)
    fname = "NP";

  TString kineTmp = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
  kineLabel = std::string(kineTmp.Data());
}

void Final2DFit::createKinematicCut()
{
  std::cout << "===== Start createKinematicCut() =====\n\n";
  TString kinePart = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);
  
  TString accPart = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
  
  TString osPart = "recoQQsign==0 &&";

  // TString anglePart = Form("&& cos_ep>%.2f && cos_ep<%.2f", cosLow, cosHigh);
  TString anglePart = "true";

  kineCut = std::string((osPart + accPart + kinePart + anglePart).Data());
}

void Final2DFit::openInputFiles()
{
  std::cout << "===== Start openInputFiles() =====\n\n";
  f1 = new TFile("../../../input_roodataset/roots/OniaRooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw0_Accw0_PtW0_TnP0_Run3_PbPb_ptWeightFit.root");
  fMass = new TFile(Form("roots/2DFit_%s/MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  fCErr = new TFile(Form("roots/2DFit_%s/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  fCRes = new TFile(Form("roots/2DFit_%s/CtauResResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  fCBkg = new TFile(Form("roots/2DFit_%s/CtauBkgResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  fCTrue = new TFile(Form("roots/2DFit_%s/CtauTrueResult_Inclusive_%s.root", DATE.c_str(), kineLabel.c_str()));
}

void Final2DFit::setupWorksapceAndData()
{
  std::cout << "===== Start setupWorksapceAndData() =====\n\n";
  RooDataSet *dataset = (RooDataSet *)f1->Get("dataset");
  RooDataSet *datasetMass = (RooDataSet *)fMass->Get("datasetMass");
  // RooAddPdf* pdfMASS_Tot = (RooAddPdf*)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Bkg = (RooDataSet *)fCErr->Get("dataw_Bkg");
  // RooDataSet *dataw_Sig = (RooDataSet*)fCErr->Get("dataw_Sig");
  RooHistPdf *pdfCTAUERR_Tot = (RooHistPdf *)fCErr->Get("pdfCTAUERR_Tot");
  RooHistPdf *pdfCTAUERR_Jpsi = (RooHistPdf *)fCErr->Get("pdfCTAUERR_Jpsi");
  RooHistPdf *pdfCTAUERR_Bkg = (RooHistPdf *)fCErr->Get("pdfCTAUERR_Bkg");
  // Todo: Res model is in the Bkg model. Also don't need to open fCRes file?
  // RooAddPdf *GaussModel_Tot = (RooAddPdf *)fCRes->Get("GaussModel_Tot");
  RooAddPdf *TrueModel_Tot = (RooAddPdf *)fCTrue->Get("TrueModel_Tot");
  RooAddPdf *pdfCTAU_Bkg_Tot = (RooAddPdf *)fCBkg->Get("pdfCTAU_Bkg_Tot");

  ws = new RooWorkspace("workspace");

  ws->import(*dataset); // total
  ws->import(*datasetMass);
  ws->import(*dataw_Bkg);
  ws->import(*TrueModel_Tot);
  ws->import(*pdfCTAUERR_Tot);
  ws->import(*pdfCTAUERR_Jpsi);
  ws->import(*pdfCTAU_Bkg_Tot);

  // make a dataset for fit
  RooArgSet argSet(*(ws->var("ctau3D")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("ctau3DRes")), *(ws->var("ctau3DErr")));
  argSet.add(*(ws->var("pt1")));
  argSet.add(*(ws->var("pt2")));
  argSet.add(*(ws->var("eta1")));
  argSet.add(*(ws->var("eta2")));
  // argSet.add(*(ws->var("cos_ep")));
  argSet.add(*(ws->var("recoQQsign")));
  argSet.add(*(ws->var("cBin")));

  RooDataSet *datasetW = new RooDataSet("datasetW", "A sample", argSet, Import(*dataset), WeightVar(*ws->var("weight")));
  ws->import(*datasetW);

  // Todo: Do we need it?
  // RooDataSet *datasetWo = new RooDataSet("datasetWo", "A sample", *argSet, Import(*dataset));
  // ws->import(*datasetWo);

  RooDataSet *dsTot = (RooDataSet *)datasetW->reduce(argSet, kineCut.c_str());
  dsTot->SetName("dsTot");
  ws->import(*dsTot);

  delete datasetW;
  delete dsTot;
  // Don't tocuh the instances got from Get() method
}

void Final2DFit::setVariableRanges()
{
  std::cout << "===== Start setVariableRanges() =====\n\n";
  
  // ctauErr
  RooDataSet *dsTot = (RooDataSet *)ws->data("dsTot");
  RooRealVar *ctauErrVar = ws->var("ctau3DErr");
  double ctauErrMin;
  double ctauErrMax;
  dsTot->getRange(*ctauErrVar, ctauErrMin, ctauErrMax);
  ctauErrVar->setRange(ctauErrMin, ctauErrMax);
  ctauErrVar->setRange("ctauRange", ctauErrMin, ctauErrMax);

  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);  
  
  // Todo: get the range from sPlot dataset?
  ws->var("ctau3DRes")->setRange(-10, 10);
  ws->var("ctau3DRes")->setRange("ctauResRange", -10, 10);

  std::cout << "--- Variable Range Set ---\n";
  ws->var("mass")->Print();
  ws->var("ctau3D")->Print();
  ws->var("ctau3DErr")->Print();
  ws->var("ctau3DRes")->Print();
  std::cout << "--------------------------\n";

  // don't free the instances got from ws->data, var, pdf etc.
}

void Final2DFit::defineAndFitMassModel()
{
  std::cout << "===== Start defineAndFitMassModel() =====\n\n";
  
  // Todo: move to define 2DFit model, user custom
  ws->factory("b_Jpsi[0.22, 1e-8, 1.0]");

  std::cout << "*************************** START MASS FIT *****************************\n";

  // ===== use MC mass fit results =====
  TFile *f_fit = new TFile(Form("roots/2DFit_%s/mc_Mass/mc_MassFitResult_%s_PRw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fEffW, fAccW, isPtW, isTnP));
  if (!f_fit || f_fit->IsZombie())
  {
    std::cerr << "Error: Could not open MC fit result file!\n";
    delete f_fit;
    return;
  }

  RooDataSet *dataset_fit = (RooDataSet *)f_fit->Get("datasetMass");
  RooWorkspace *ws_fit = new RooWorkspace("workspace_fit");
  ws_fit->import(*dataset_fit);

  Double_t alpha_MC_value = ws_fit->var("alpha_1_A")->getVal();
  Double_t alpha_MC_value_err = ws_fit->var("alpha_1_A")->getError();
  Double_t n_MC_value = ws_fit->var("n_1_A")->getVal();
  Double_t n_MC_value_err = ws_fit->var("n_1_A")->getError();
  Double_t xA_MC_value = ws_fit->var("x_A")->getVal();
  Double_t xA_MC_value_err = ws_fit->var("x_A")->getError();
  Double_t f_MC_value = ws_fit->var("f")->getVal();
  Double_t f_MC_value_err = ws_fit->var("f")->getError();
  Double_t sigma_MC_value = ws_fit->var("sigma_1_A")->getVal();
  Double_t sigma_MC_value_err = ws_fit->var("sigma_1_A")->getError();

  double sigma_index = 5;

  Double_t alpha_lower = alpha_MC_value - (sigma_index * alpha_MC_value_err);
  Double_t alpha_higher = alpha_MC_value + (sigma_index * alpha_MC_value_err);
  Double_t xA_lower = xA_MC_value - (sigma_index * xA_MC_value_err);
  Double_t xA_higher = xA_MC_value + (sigma_index * xA_MC_value_err);
  Double_t n_lower = n_MC_value - (sigma_index * n_MC_value_err);
  Double_t n_higher = n_MC_value + (sigma_index * n_MC_value_err);
  Double_t f_lower = f_MC_value - (sigma_index * f_MC_value_err);
  Double_t f_higher = f_MC_value + (sigma_index * f_MC_value_err);
  Double_t sigma_lower = sigma_MC_value - (sigma_index * sigma_MC_value_err);
  Double_t sigma_higher = sigma_MC_value + (sigma_index * sigma_MC_value_err);

  if (f_higher > 1.0)
    f_higher = 1.0;
  if (f_lower < 0.0)
    f_lower = 0.0;
  double paramslower[6] = {alpha_lower, n_lower, 1e-6, xA_lower, 0.0, 0.0};
  double paramsupper[6] = {alpha_higher, n_higher, 0.5, xA_higher, 1.0, 25.0};

  double alpha_1_init = alpha_MC_value;
  double n_1_init = n_MC_value;
  double sigma_1_init = sigma_MC_value;
  double x_init = xA_MC_value;
  double f_init = f_MC_value;

  // ===== Sig model =====
  RooRealVar mean("mean_Jpsi", "mean of the signal", pdgMass.JPsi, pdgMass.JPsi - 0.1, pdgMass.JPsi + 0.1);
  RooRealVar sigma_1_A("sigma_1_A", "sigma of CB1", sigma_MC_value, sigma_lower, sigma_higher);
  RooRealVar x_A("x_A", "sigma ratio", xA_MC_value);
  RooFormulaVar sigma_2_A("sigma_2_A", "@0*@1", RooArgList(sigma_1_A, x_A));
  RooRealVar alpha_1_A("alpha_1_A", "tail shift", alpha_MC_value);
  RooFormulaVar alpha_2_A("alpha_2_A", "1.0*@0", RooArgList(alpha_1_A));
  RooRealVar n_1_A("n_1_A", "power order", n_MC_value);
  RooFormulaVar n_2_A("n_2_A", "1.0*@0", RooArgList(n_1_A));
  RooRealVar f("f", "cb fraction", f_MC_value);

  RooCBShape cb_1_A("cball_1_A", "Crystal Ball 1", *ws->var("mass"), mean, sigma_1_A, alpha_1_A, n_1_A);
  RooCBShape cb_2_A("cball_2_A", "Crystal Ball 2", *ws->var("mass"), mean, sigma_2_A, alpha_2_A, n_2_A);
  RooAddPdf pdfMASS_Jpsi("pdfMASS_Jpsi", "Signal Shape", RooArgList(cb_1_A, cb_2_A), f);

  // ===== Bkg model =====
  RooRealVar sl1("sl1", "sl1", 0.0, -1., 1.);
  RooRealVar sl2("sl2", "sl2", 0.0, -1., 1.);
  // RooRealVar *sl3 = new RooRealVar("sl3", "sl3", 0.0, -1., 1.);
  // RooRealVar *sl4 = new RooRealVar("sl4", "sl4", 0.0, -1., 1.);
  // RooRealVar *sl5 = new RooRealVar("sl5", "sl5", 0.0, -1., 1.);
  // RooRealVar *sl6 = new RooRealVar("sl6", "sl6", 0.0, -1., 1.);
  RooChebychev pdfMASS_bkg("pdfMASS_bkg", "Background Shape", *ws->var("mass"), RooArgList(sl1, sl2));

  // ===== Sig + Bkg model =====
  Double_t NBkg_limit = 2.0e+07;
  Double_t NJpsi_limit = 10.0e+06;

  RooRealVar N_Jpsi("N_Jpsi", "Number of Jpsi", 10000, 1, 1e7);
  RooRealVar N_Bkg("N_Bkg", "Number of Background", 50000, 1, 2e7);
  RooAddPdf pdfMASS_Tot("pdfMASS_Tot", "", RooArgList(pdfMASS_Jpsi, pdfMASS_bkg), RooArgList(N_Jpsi, N_Bkg));

  ws->import(pdfMASS_Tot);

  fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*ws->data("dsTot"), Save(), Range(massLow, massHigh), Timer(kTRUE), Extended(kTRUE), SumW2Error(true), NumCPU(nCPU), PrintLevel(0));
  fitMass->Print("V");
  std::cout << "**************************** END MASS FIT ******************************\n";

  delete ws_fit;
  f_fit->Close(); delete f_fit;
}

void Final2DFit::defineCtauModel()
{
  std::cout << "===== Start defineCtauModels() =====\n\n";
  ws->pdf("pdfCTAU_Bkg_Tot")->getParameters(ws->data("dsTot")->get())->setAttribAll("Constant", kTRUE);

  // ws->pdf("pdfCTAU_Bkg_Tot")->getParameters(
  //                               // RooArgSet(*ws->var("ctau3D"), *ws->pdf("pdfCTAU_BkgPR"), *ws->pdf("pdfCTAU_BkgNoPR"),  *ws->pdf("pdfCTAUCOND_Bkg"), *ws->pdf("pdfCTAUCOND_BkgPR"), *ws->pdf("pdfCTAUCOND_BkgNoPR")
  //                               RooArgSet(*ws->var("ctau3D"), *ws->pdf("pdfCTAU_BkgPR"), *ws->pdf("pdfCTAU_BkgNoPR"), *ws->pdf("pdfCTAUCOND_BkgPR"), *ws->pdf("pdfCTAUCOND_BkgNoPR")))
  //     ->setAttribAll("Constant", kTRUE);

  double lambda = ws->var("lambdaDSS")->getVal();
  double lambda1 = ws->var("lambdaDSS2")->getVal();
  //   double lambda2 = ws->var("lambdaDSS3")->getVal();
  double fdss = ws->var("fDSS")->getVal();
  //   double fdss1 = ws->var("fDSS1")->getVal();

  ws->factory(Form("lambdaDSS_test1[%.4f]", lambda));
  ws->factory(Form("lambdaDSS_test2[%.4f]", lambda1));
  //   ws->factory(Form("lambdaDSS_test3[%.4f]", lambda2));
  ws->factory(Form("fDSS1_test[%.4f]", fdss));
  //   ws->factory(Form("fDSS2_test[%.4f]", fdss1));

  // ws->var("lambdaDSS_test1")->setConstant();
  // ws->var("lambdaDSS_test2")->setConstant();
  // //   ws->var("lambdaDSS_test3")->setConstant();
  // ws->var("fDSS1_test")->setConstant();
  // //   ws->var("fDSS2_test")->setConstant();

  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test1", "ctau3D", "lambdaDSS_test1", "pdfCTAURES")); // NP
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test2", "ctau3D", "lambdaDSS_test2", "pdfCTAURES")); // NP
  //   ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test3", "ctau3D", "lambdaDSS_test3", "pdfCTAURES")); //NP
  RooAddPdf pdfCTAUCOND_JpsiNoPR("pdfCTAUCOND_JpsiNoPR", "", RooArgList(*ws->pdf("pdfCTAUTRUE_test1"), *ws->pdf("pdfCTAUTRUE_test2")), RooArgList(*ws->var("fDSS1_test")));
  ws->import(pdfCTAUCOND_JpsiNoPR);

  ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_JpsiPR", "pdfCTAURES"));

  RooProdPdf pdfJpsiPR("pdfCTAU_JpsiPR", "", *ws->pdf("pdfCTAUERR_Jpsi"),
                       Conditional(*ws->pdf("pdfCTAUCOND_JpsiPR"), RooArgList(*ws->var("ctau3D"))));
  ws->import(pdfJpsiPR);
  RooProdPdf pdfJpsiNoPR("pdfCTAU_JpsiNoPR", "", *ws->pdf("pdfCTAUERR_Jpsi"),
                         Conditional(*ws->pdf("pdfCTAUCOND_JpsiNoPR"), RooArgList(*ws->var("ctau3D"))));
  ws->import(pdfJpsiNoPR);
}

void Final2DFit::define2DFitModel()
{
  std::cout << "===== Start define2DFitModel() =====\n\n";

  // 2D Sig model
  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_JpsiPR",
                   "pdfCTAU_JpsiPR",
                   "pdfMASS_Jpsi"));
  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_JpsiNoPR",
                   "pdfCTAU_JpsiNoPR",
                   "pdfMASS_Jpsi"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUMASS_Jpsi",
                   "b_Jpsi",
                   "pdfCTAUMASS_JpsiNoPR",
                   "pdfCTAUMASS_JpsiPR"));
  
  // 2D Bkg model
  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgPR",
                   "pdfCTAU_BkgPR",
                   "pdfMASS_bkg"));
  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgNoPR",
                   "pdfCTAU_BkgNoPR",
                   "pdfMASS_bkg"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUMASS_Bkg",
                   "b_Bkg",
                   "pdfCTAUMASS_BkgNoPR",
                   "pdfCTAUMASS_BkgPR"));

  double njpsi = ws->var("N_Jpsi")->getVal();
  ws->factory(Form("nSig2D[%.3f, %.3f, %.3f]", njpsi, 1., 1e+8));
  double nbkg = ws->var("N_Bkg")->getVal();
  ws->factory(Form("nBkg2D[%.3f, %.3f, %.3f]", nbkg, 1., 1e+8));

  // 2D Sig + Bkg model
  RooAbsPdf *themodel = NULL;
  themodel = new RooAddPdf("pdfCTAUMASS_Tot", "",
                           RooArgList(*ws->pdf("pdfCTAUMASS_Jpsi"), *ws->pdf("pdfCTAUMASS_Bkg")),
                           RooArgList(*ws->var("nSig2D"), *ws->var("nBkg2D")));
  ws->import(*themodel);

  delete themodel;
}

void Final2DFit::perform2DFit()
{
  std::cout << "===== Start perform2DFit() =====\n\n";
  std::cout << "\n******************* Start Ctau-Mass 2D Fit *********************\n\n";

  fit2DFit = ws->pdf("pdfCTAUMASS_Tot")->fitTo(*ws->data("dsTot"), Extended(kTRUE), NumCPU(nCPU), SumW2Error(true), PrintLevel(0), Save());
  fit2DFit->Print("V");

  std::cout << "\n******************** Finish Ctau-Mass 2D Fit **********************\n";
}

void Final2DFit::makePlots()
{
  std::cout << "\n===== Start makePlots() =====\n";
  plotMass();
  plotCtau();
}

void Final2DFit::plotMass()
{
  RooDataSet *dsTot = (RooDataSet *)ws->data("dsTot");

  TCanvas *c_mass = new TCanvas("c_mass", "Mass Distribution", 800, 800);

  TPad *pad1 = new TPad("pad1", "", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->SetTicks(1, 1);
  pad1->Draw();
  pad1->cd();
  gPad->SetLogy();

  RooPlot *massFrame = ws->var("mass")->frame(RooFit::Title("Mass Distribution"), Bins(nMassBin), Range(massLow, massHigh));
  dsTot->plotOn(massFrame, RooFit::Name("data_mass"));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(massFrame, RooFit::Name("fit_mass"));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(massFrame, RooFit::Components("pdfMASS_bkg"), RooFit::LineStyle(kDashed), RooFit::LineColor(kAzure - 9));
  massFrame->GetYaxis()->SetTitleSize(0.05);
  massFrame->GetYaxis()->SetLabelSize(0.04);
  massFrame->GetYaxis()->SetTitleOffset(1.0);

  TH1 *h_tmp = dsTot->createHistogram("h_tmp", *ws->var("mass"), Binning(nMassBin, massLow, massHigh));
  Double_t yMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
  Double_t yMin = 1e99;
  for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
    if (h_tmp->GetBinContent(i) > 0)
      yMin = std::min(yMin, h_tmp->GetBinContent(i));
  Double_t yUp(0.), yDown(0.);
  yUp = yMax*TMath::Power((yMax/0.1), 0.3);
  // yUp = yMax * TMath::Power((yMax / 0.01), 0.5);
  yDown = yMin / (TMath::Power((yMax / yMin), (0.2 / (1.0 - 0.9))));
  massFrame->GetYaxis()->SetRangeUser(yDown, yUp);
  massFrame->Draw();
  delete h_tmp;

  c_mass->cd();
  TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.25);
  pad2->SetTicks(1, 1);
  pad2->Draw();
  pad2->cd();

  RooPlot *frameTMP = (RooPlot *)massFrame->Clone("frameTMP");
  RooHist *hPull = massFrame->pullHist("data_mass", "fit_mass", true);
  hPull->SetMarkerSize(0.8);
  RooPlot *pullFrame = ws->var("mass")->frame(Title("Pull Distribution"));
  pullFrame->addPlotable(hPull, "P");
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

  printChi2(ws, pad2, frameTMP, fit2DFit, "mass", "data_mass", "fit_mass", nMassBin, false);

  TLine *l1 = new TLine(massLow, 0, massHigh, 0);
  l1->SetLineStyle(1);
  l1->Draw("same");

  TLine *line_at_0 = new TLine(ws->var("mass")->getMin(), 0.0, ws->var("mass")->getMax(), 0.0);
  line_at_0->SetLineColor(kRed);
  line_at_0->SetLineStyle(kDashed);
  line_at_0->Draw("same");

  c_mass->SaveAs(Form("figs/2DFit_%s/Fit_Mass_%s.png", DATE.c_str(), kineLabel.c_str()));
  c_mass->SaveAs(Form("../figs/2DFit_%s/Final/Fit_Mass_%s.pdf", DATE.c_str(), kineLabel.c_str()));
  delete frameTMP;
  delete c_mass;
}

void Final2DFit::plotCtau()
{
  RooDataSet *dsTot = (RooDataSet *)ws->data("dsTot");

  TCanvas *c_ctau = new TCanvas("c_ctau", "Lifetime Distribution", 800, 800);
  c_ctau->cd();

  TPad *pad1 = new TPad("pad1", "", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->SetTicks(1, 1);
  pad1->Draw();
  pad1->cd();
  gPad->SetLogy();

  RooPlot *ctauFrame = ws->var("ctau3D")->frame(RooFit::Title("Lifetime Distribution"), Bins(nCtauBins), Range(ctauLow, ctauHigh));
  ctauFrame->updateNormVars(RooArgSet(*ws->var("ctau3D")));
  dsTot->plotOn(ctauFrame, RooFit::Name("data_ctau"));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(ctauFrame, RooFit::Name("fit_ctau"), RooFit::ProjWData(*dsTot));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(ctauFrame, RooFit::Components("pdfCTAUMASS_Bkg"), RooFit::LineStyle(kDashed), RooFit::LineColor(kAzure - 9), RooFit::ProjWData(*dsTot));
  // ws->pdf("pdfCTAUMASS_Tot")->plotOn(ctauFrame, RooFit::ProjWData(*dsTot), RooFit::Components("pdfCTAUMASS_Bkg"), RooFit::LineStyle(kDashed), RooFit::LineColor(kAzure - 9));
  ctauFrame->GetYaxis()->SetTitleSize(0.05);
  ctauFrame->GetYaxis()->SetRangeUser(0.01, 1000000);
  ctauFrame->GetYaxis()->SetLabelSize(0.04);
  ctauFrame->GetYaxis()->SetTitleOffset(1.0);

  TH1 *h_tmp = dsTot->createHistogram("h_tmp", *ws->var("ctau3D"), Binning(nCtauBins, ctauLow, ctauHigh));
  Double_t yMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
  Double_t yMin = 1e99;
  for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
    if (h_tmp->GetBinContent(i) > 0)
      yMin = std::min(yMin, h_tmp->GetBinContent(i));
  Double_t yUp(0.), yDown(0.);
  yUp = yMax * TMath::Power((yMax / 0.1), 0.3);
  // yUp = yMax * TMath::Power((yMax / 0.01), 0.5);
  yDown = yMin / (TMath::Power((yMax / yMin), (0.2 / (1.0 - 0.9))));
  ctauFrame->GetYaxis()->SetRangeUser(yDown, yUp);
  ctauFrame->Draw();
  delete h_tmp;

  c_ctau->cd();
  TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.25);
  pad2->SetTicks(1, 1);
  pad2->Draw();
  pad2->cd();

  RooHist *h_ctau = (RooHist *)ctauFrame->findObject("data_ctau");
  RooCurve *c_fit = (RooCurve *)ctauFrame->findObject("fit_ctau");
  RooPlot *ratioFrame = ws->var("ctau3D")->frame(RooFit::Title(" "));

  TGraphAsymmErrors *g_ratio_ctau = new TGraphAsymmErrors();
  int point_idx_ctau = 0;
  for (int i = 0; i < h_ctau->GetN(); ++i)
  {
    double x, y;
    h_ctau->GetPoint(i, x, y);
    if (y == 0)
      continue;
    double model_val = c_fit->Eval(x);
    if (model_val > 1e-9)
    {
      g_ratio_ctau->SetPoint(point_idx_ctau, x, y / model_val);
      g_ratio_ctau->SetPointError(point_idx_ctau, 0, 0, h_ctau->GetErrorYlow(i) / model_val, h_ctau->GetErrorYhigh(i) / model_val);
      point_idx_ctau++;
    }
  }

  ratioFrame->GetYaxis()->SetRangeUser(0, 2);
  ratioFrame->GetYaxis()->SetTitle("Data / Fit");
  ratioFrame->GetYaxis()->SetNdivisions(505);
  ratioFrame->GetYaxis()->SetTitleSize(0.1);
  ratioFrame->GetYaxis()->SetLabelSize(0.08);
  ratioFrame->GetYaxis()->SetTitleOffset(0.4);
  ratioFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  ratioFrame->GetXaxis()->SetTitleSize(0.1);
  ratioFrame->GetXaxis()->SetLabelSize(0.08);
  ratioFrame->Draw();

  TLine *line_at_1_ctau = new TLine(ws->var("ctau3D")->getMin(), 1.0, ws->var("ctau3D")->getMax(), 1.0);
  line_at_1_ctau->SetLineColor(kRed);
  line_at_1_ctau->SetLineStyle(kDashed);
  line_at_1_ctau->Draw("same");

  g_ratio_ctau->SetMarkerStyle(20);
  g_ratio_ctau->SetMarkerSize(0.8);
  g_ratio_ctau->Draw("P SAME");

  c_ctau->SaveAs(Form("figs/2DFit_%s/Fit_Ctau_%s.png", DATE.c_str(), kineLabel.c_str()));
  c_ctau->SaveAs(Form("../figs/2DFit_%s/Final/Fit_Ctau_%s.pdf", DATE.c_str(), kineLabel.c_str()));
  delete c_ctau;
}

void Final2DFit::saveResults()
{
  std::cout << "\n===== Start saveResults() =====\n";
  TFile *outputFile = new TFile(Form("roots/2DFit_%s/FitResult_%s.root", DATE.c_str(), kineLabel.c_str()), "RECREATE");
  outputFile->cd();

  RooWorkspace wsTmp("ws_final");
  wsTmp.import(*ws->pdf("pdfCTAUMASS_Tot"));
  wsTmp.Write();

  fit2DFit->Write("fit2DFit");

  outputFile->Close();
  delete outputFile;
}