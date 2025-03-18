#include <iostream>
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"
#include "../headers/rootFitHeaders.h"
#include "../headers/commonUtility.h"
#include "../headers/JpsiUtility.h"
#include "../headers/cutsAndBin.h"
#include "../headers/CMS_lumi_v2mass.C"
#include "../headers/tdrstyle.C"

void final()
{
    TStopwatch *t = new TStopwatch;
    t->Start();
    cout << "\n=================================\n";
    cout << "\n Start 2D final fit\n";
    cout << "\n=================================\n";

    // set macro configure
    using namespace std;
    using namespace RooFit;

    // input parameters
    float cos_low = 0.0, cos_high = 0.1;
    float ptLow = 3;
    float ptHigh = 6.5;
    float yLow = 1.6;
    float yHigh = 2.4;
    int cLow = 60;
    int cHigh = 180;

    // Usually not used
    int PR = 0; // 0=PR, 1=NP, 2=Inc.
    int PRw = 1;
    bool fEffW = true;
    bool fAccW = true;
    bool isPtW = true;
    bool isTnP = true;
    double massLow = 2.6;
    double massHigh = 3.5; // Jpsi mass range

    // be quiet please
    // RooMsgService::instance().getStream(0).removeTopic(Caching); // in RooFit:: namespace
    // RooMsgService::instance().getStream(1).removeTopic(Caching);
    // RooMsgService::instance().getStream(0).removeTopic(Plotting);
    // RooMsgService::instance().getStream(1).removeTopic(Plotting);
    // RooMsgService::instance().getStream(0).removeTopic(Integration);
    // RooMsgService::instance().getStream(1).removeTopic(Integration);
    // RooMsgService::instance().getStream(0).removeTopic(InputArguments);
    // RooMsgService::instance().getStream(1).removeTopic(InputArguments);
    // RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    // RooMsgService::instance().getStream(1).removeTopic(Minimization);
    // RooMsgService::instance().getStream(1).removeTopic(Caching);
    RooMsgService::instance().getStream(0).removeTopic(Tracing);
    RooMsgService::instance().getStream(1).removeTopic(Tracing);
    RooMsgService::instance().setGlobalKillBelow(WARNING);


    // ===== labeling ===== //
    std::ostringstream out_file_path;
    out_file_path << "final/final_pt" << ptLow << "-" << ptHigh
                  << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
                  << "_cos" << cos_low << "_" << cos_high;
    string out_ss = out_file_path.str();

    std::ostringstream input_file;
    input_file << "pt" << ptLow << "-" << ptHigh
               << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
               << "_cos" << cos_low << "_" << cos_high
               << ".root";
    string in_ss = input_file.str();

    std::ostringstream input_true_file;
    input_true_file << "pt" << ptLow << "-" << ptHigh
                    << "_absY" << yLow << "-" << yHigh
                    << ".root";
    string in_ture_ss = input_true_file.str();

    // ===== make output folders ===== //
    gSystem->mkdir(Form("roots/final/"), kTRUE);
    gSystem->mkdir(Form("figs/final"), kTRUE);


    // ===== kinematic cut ===== //
    TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);
    TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
    TString OS = "recoQQsign==0 &&";

    TString angle_cut = Form("&& cos_ep>%.2f && cos_ep<%.2f", cos_low, cos_high);

    TString final_cut = OS + accCut + kineCut + angle_cut;


    // ===== import files ===== //
    auto f_data = new TFile("../files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root");
    auto f_mass = new TFile(("roots/mass/mass_" + in_ss).c_str());
    auto f_err = new TFile(("roots/err/err_" + in_ss).c_str());
    auto f_res = new TFile(("roots/res/res_" + in_ss).c_str());
    auto f_bkg = new TFile(("roots/bkg/bkg_" + in_ss).c_str());
    auto f_true = new TFile(("roots/true/true_" + in_ture_ss).c_str());

    // import objects
    auto ws_mass = (RooWorkspace *)f_mass->Get("ws_mass");
    auto ws_err = (RooWorkspace *)f_err->Get("ws_err");
    auto ws_res = (RooWorkspace *)f_res->Get("ws_res");
    auto ws_true = (RooWorkspace *)f_true->Get("ws_true");
    auto ws_bkg = (RooWorkspace *)f_bkg->Get("ws_bkg");
    // // RooAddPdf *mass_pdf = (RooAddPdf *)((RooWorkspace *)f_mass->Get("ws_mass"))->pdf("mass_pdf");

    

    // ===== make a dataset ===== //
    RooDataSet *ds_tmp = (RooDataSet *)f_data->Get("dataset");
    ds_tmp->SetName("ds_tmp");
    RooDataSet *ds_tmp_weight = new RooDataSet("ds_tmp_weight", "", *ds_tmp->get(), Import(*ds_tmp), WeightVar("weight"));

    // tmporal local variables
    // auto mass = (RooRealVar *)ds_tmp_weight->get(0)->find("mass");
    // auto ctau3D = (RooRealVar *)ds_tmp_weight->get(0)->find("ctau3D");
    // auto ctau3DErr = (RooRealVar *)ds_tmp_weight->get(0)->find("ctau3DErr");

    double ctauErrMin = ws_err->var("ctau3DErr")->getMin();
    double ctauErrMax = ws_err->var("ctau3DErr")->getMax();

    auto ws_tmp = new RooWorkspace("ws_tmp");
    ws_tmp->import(*ds_tmp);
    ws_tmp->import(*ds_tmp_weight);

    // auto ds_final = (RooDataSet *)ds_tmp_weight->reduce(RooArgSet(*ws_tmp->var("mass"), *ws_tmp->var("ctau3D"), *ws_tmp->var("ctau3DErr")), Form("ctau3DErr>%f&&ctau3DErr<%f", ctauErrMin, ctauErrMax));
    auto ds_final = (RooDataSet *)ds_tmp->reduce(RooArgSet(*ws_tmp->var("mass"), *ws_tmp->var("ctau3D"), *ws_tmp->var("ctau3DErr")), Form("ctau3DErr>%f&&ctau3DErr<%f", ctauErrMin, ctauErrMax));
    ds_final->SetName("ds_final");

    ws_tmp->import(*ds_final);

    auto ws_final = new RooWorkspace("ws_final");
    ws_final->import(*ds_final);

    // // re-connet pointer to ds_final observables
    // mass = (RooRealVar *)ds_final->get(0)->find("mass");
    // ctau3D = (RooRealVar *)ds_final->get(0)->find("ctau3D");
    // ctau3DErr = (RooRealVar *)ds_final->get(0)->find("ctau3DErr");

    // relase memories of temporal pointers
    // delete mass;
    // delete ctau3D;
    // delete ctau3DErr;
    
    // set variable ranges
    ws_final->var("mass")->setRange(massLow, massHigh);
    ws_final->var("ctau3D")->setRange(ctauLow, ctauHigh);
    ws_final->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
    // ws_final->data("ds_final")->Print("V");


    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Build mass model *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Clone() is deep copy. It makes variable which disconnect to original objects.
    auto mass_sig = (RooAddPdf *)ws_mass->pdf("mass_sig")->Clone("mass_sig");
    auto mass_bkg = (RooAddPdf *)ws_mass->pdf("mass_bkg")->Clone("mass_bkg");
    ws_final->import(*mass_sig);
    ws_final->import(*mass_bkg);

    
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Bring Punzi terms *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    auto sig_punzi = (RooHistPdf *)ws_err->pdf("err_sig_pdf")->Clone("sig_punzi");
    auto bkg_punzi = (RooHistPdf *)ws_err->pdf("err_bkg_pdf")->Clone("bkg_punzi");
    ws_final->import(*sig_punzi);
    ws_final->import(*bkg_punzi);

    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Build ctauBkg model *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    
    // np
    auto ctau_NP = (RooAddModel *)ws_bkg->pdf("pdfCTAUCOND_BkgNoPR")->Clone("ctau_NP");
    ws_final->import(*ctau_NP);

    // pr
    ws_final->factory(Form("SUM::%s(%s)", "ctau_PR", "pdfCTAURES"));

    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Build res model *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // it was imported when we imported ctau_NP -> Use pdfCTAURES
    // just for test
    // auto res_gauss = (RooAddModel *)ws_final->pdf("pdfCTAURES")->Clone("res_gauss");
    // res_gauss->Print("V");

    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Build NP model *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    auto true_np = (RooAddPdf *)ws_true->pdf("true_decay")->Clone("ctau_np");
    ws_final->import(*true_np);

    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Product Puniz term x ctau PDFs *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // RooProdPdf ctau_PR("ctau_PR", "CtPDF with PEE", *ws_final->pdf("sig_punzi"), Conditional(*ws_final->pdf("res_gauss"), RooArgList(*ws_final->var("ctau3D"))));
    // RooProdPdf ctau_sig_NP("ctau_sig_NP", "CtPDF with PEE", *ws_final->pdf("sig_punzi"), Conditional(*ws_final->pdf("ctau_NP"), RooArgList(*ws_final->var("ctau3D"))));

    // RooRealVar f_b_bkg("f_b_bkg", "bfraction in bkg - its not the b_frac we measure", 0.8, 0.3, 1);
    // RooProdPdf CtBkg_PR("CtBkg_PR", "CtPDF with PEE", *err_pdf_bkg, Conditional(res_gauss, RooArgList(*ctau3D)));                
    // RooProdPdf CtBkg_NP("CtBkg_NP", "CtPDF with PEE", *err_pdf_bkg, Conditional(decay_bkg, RooArgList(*ctau3D)));

    RooProdPdf ctau_PR_pee("ctau_PR_pee", "", *ws_final->pdf("sig_punzi"),
                       Conditional(*ws_final->pdf("ctau_PR"), RooArgList(*ws_final->var("ctau3D"))));
    ws_final->import(ctau_PR_pee);

    RooProdPdf ctau_NP_pee("ctau_NP_pee", "", *ws_final->pdf("bkg_punzi"),
                       Conditional(*ws_final->pdf("ctau_NP"), RooArgList(*ws_final->var("ctau3D"))));
    ws_final->import(ctau_NP_pee);


    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Product mass x ctau - sig and bkg individually *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // RooProdPdf mass_ct_sig_PR("mass_ct_sig_PR", "mass x ctau3D", RooArgList(mass_sig, ctau_PR));
    // RooProdPdf mass_ct_sig_NP("mass_ct_sig_NP", "mass x ctau3D", RooArgList(mass_sig, ctau_sig_NP));
    // RooAddPdf mass_ct_sig("mass_ct_sig", "", RooArgSet(mass_ct_sig_NP, mass_ct_sig_PR), b_frac);

    // RooProdPdf mass_ct_bkg_PR("mass_ct_bkg_PR", "mass x ctau3D", RooArgList(mass_bkg, CtBkg_PR));
    // RooProdPdf mass_ct_bkg_NP("mass_ct_bkg_NP", "mass x ctau3D", RooArgList(mass_bkg, CtBkg_NP));
    // RooAddPdf mass_ct_bkg("mass_ct_bkg", "", RooArgSet(mass_ct_bkg_NP, mass_ct_bkg_PR), f_b_bkg);


    // ===== The b_fraction ===== //
    ws_final->factory("b_frac[0.25, 0.2, 1.0]");

    // bkg
    ws_final->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgPR",
                           "ctau_PR_pee",
                           "mass_bkg"));
    ws_final->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgNoPR",
                           "ctau_NP_pee",
                           "mass_bkg"));
    ws_final->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUMASS_Bkg",
                     "b_Bkg[0.5,0.001,1]",
                     "pdfCTAUMASS_BkgNoPR",
                     "pdfCTAUMASS_BkgPR")); // b_Bkg is bkg PR, NP ratio

    // sig
    ws_final->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_JpsiPR",
                     "ctau_PR_pee",
                     "mass_sig"));
    ws_final->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_JpsiNoPR",
                     "ctau_NP_pee",
                     "mass_sig"));
    ws_final->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUMASS_Jpsi",
                     "b_frac",
                     "pdfCTAUMASS_JpsiNoPR",
                     "pdfCTAUMASS_JpsiPR"));

    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Final fit model *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    auto n_sig = new RooRealVar("n_sig", "", 5000, 1, 25000);
    auto n_bkg = new RooRealVar("n_bkg", "", 20000, 1, 100000);
    ws_final->import(*n_sig);
    ws_final->import(*n_bkg);

    RooAbsPdf *themodel = NULL;
    themodel = new RooAddPdf("themodel", "themodel",
                             RooArgList(*ws_final->pdf("pdfCTAUMASS_Bkg"), *ws_final->pdf("pdfCTAUMASS_Jpsi")),
                             RooArgList(*ws_final->var("n_bkg"), *ws_final->var("n_sig")));
    ws_final->import(*themodel);

    // RooRealVar N_JPSI("N_JPSI", "bfraction in bkg - it's not the b_frac we measure", 100000, 1, 1e7);
    // RooExtendPdf pdfTot_sig("pdfTot_sig", "", mass_ct_sig, N_JPSI);

    // RooRealVar N_BKG("N_BKG", "bfraction in bkg - it's not the b_frac we measure", 1e5, 1, 1e6);
    // RooExtendPdf pdfTot_bkg("pdfTot_bkg", "", mass_ct_bkg, N_BKG);

    // RooFormulaVar nn_frac("nn_frac", "@0/(@0+@1)", RooArgSet(N_JPSI, N_BKG));
    // RooAddPdf final_2d_model("final_2d_model", "Use it for 2D fit", RooArgSet(pdfTot_sig, pdfTot_bkg), RooArgSet(nn_frac));
    // // RooAddPdf final_2d_model("final_2d_model", "Use it for 2D fit", RooArgSet(mass_ct_sig_PEE, mass_ct_bkg_PEE), RooArgSet(*n_sig, *n_bkg));


    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Fit here *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    auto fit_2d = ws_final->pdf("themodel")->fitTo(*ds_final, Save(1), SumW2Error(true), NumCPU(12), Extended(1), Strategy(2), RecoverFromUndefinedRegions(1.), PrintEvalErrors(-1));
    fit_2d->Print("V");


    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Plotting - very long but simple *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    // two parts: mass and ctau3D

    // ===== main canvas - mass ===== //
    auto c_mass = new TCanvas("c_mass", "", 800, 800);
    
    // ===== mass distribution ===== //
    auto mass_frame = ws_tmp->var("mass")->frame(Range(massLow, massHigh));
    ws_tmp->data("ds_final")->plotOn(mass_frame, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), MarkerColor(kRed));
    ws_tmp->data("ds_tmp")->plotOn(mass_frame, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), MarkerColor(kBlack));
    themodel->plotOn(mass_frame);
    mass_frame->Draw();
    c_mass->Draw();

    // ===== main canvas - ctau3D ===== //
    auto c_ctau = new TCanvas("c_ctau", "", 800, 800);
    c_ctau->cd();

    // ===== ctau3D dist ===== //
    auto ctau_frame = ws_tmp->var("ctau3D")->frame(Range(ctauLow, ctauHigh));
    ctau_frame->updateNormVars(RooArgSet(*ws_tmp->var("ctau3D")));
    ws_tmp->data("ds_final")->plotOn(ctau_frame, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), MarkerColor(kRed));
    ws_tmp->data("ds_tmp")->plotOn(ctau_frame, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), MarkerColor(kBlack));
    themodel->plotOn(ctau_frame);
    ctau_frame->Draw();
    c_ctau->Draw();

    // ===== ===== //
    // ===== ===== //
    // ===== ===== //
    // ===== ===== //


    cout << "\n=================================\n";
    cout << "\n Finish 2D final fit\n";
    cout << "\n=================================\n";
    t->Stop();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
    ds_tmp->Print("V");
}

void to_store_test_codes_please_ignore()
{
    // ===== ===== //
}