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

void true_fit()
{
    TStopwatch *t = new TStopwatch;
    t->Start();
    cout << "\n=================================\n";
    cout << "\n Start ctau true fit\n";
    cout << "\n=================================\n";

    double nCtauBins = 1000; // for test

    // set macro configure
    using namespace std;
    using namespace RooFit;

    // input parameters
    float cos_low = 0.0, cos_high = 0.1;
    float ptLow = 3;
    float ptHigh = 6.5;
    float yLow = 1.6;
    float yHigh = 2.4;
    int cLow = 0;
    int cHigh = 180;

    // Usually not used
    int PR = 0; // 0=PR, 1=NP, 2=Inc.
    int PRw = 1;
    bool fEffW = false;
    bool fAccW = false;
    bool isPtW = false;
    bool isTnP = false;
    double massLow = 2.6;
    double massHigh = 3.5; // Jpsi mass range

    // be quiet please
    // RooMsgService::instance().getStream(0).removeTopic(Caching); // in RooFit:: namespace
    // RooMsgService::instance().getStream(1).removeTopic(Caching);
    // RooMsgService::instance().getStream(0).removeTopic(Plotting);
    // RooMsgService::instance().getStream(1).removeTopic(Plotting);
    // RooMsgService::instance().getStream(0).removeTopic(InputArguments);
    // RooMsgService::instance().getStream(1).removeTopic(InputArguments);
    // RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    // RooMsgService::instance().getStream(1).removeTopic(Minimization);
    // RooMsgService::instance().getStream(1).removeTopic(Caching);
    RooMsgService::instance().setGlobalKillBelow(WARNING);


    // ===== kinematic cut ===== //
    TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f", ptLow, ptHigh, yLow, yHigh, massLow, massHigh);
    TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut

    TString final_cut = accCut + kineCut;


    // ===== labeling ===== //
    std::ostringstream out_file_path;
    out_file_path << "true/true_pt" << ptLow << "-" << ptHigh
                  << "_absY" << yLow << "-" << yHigh;
    // << "_cos" << cos_low << "_" << cos_high;
    string out_ss = out_file_path.str();


    // ===== make output folders ===== //
    gSystem->mkdir(Form("roots/true/"), kTRUE);
    gSystem->mkdir(Form("figs/true"), kTRUE);


    // ===== import files ===== //
    auto f_np_mc = new TFile("../files_roodata/RooDataSet_miniAOD_isMC1_NP_Jpsi_GenOnly_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250307.root");
    f_np_mc->SetName("f_np_mc");
    RooDataSet *dataset = (RooDataSet *)f_np_mc->Get("dataset"); // full NP MC dataset -> Only used in reducing step

    auto ds_temp_weight = new RooDataSet("ds_temp_weight", "", *dataset->get(), Import(*dataset)); // set weight if exist
    auto ds_temp_kinematics = (RooDataSet *)ds_temp_weight->reduce(final_cut.Data());

    auto ws_true = new RooWorkspace("ws_true");
    

    // ===== set ranges ===== //
    // set parameter range
    double trueMin = 1e-6,
           trueMax = 10;

    // make a reduced dataset to use for fitting and plotting
    RooDataSet *ds_red_np_mc = (RooDataSet *)ds_temp_kinematics->reduce(Form("ctau3D>=%.4f&&ctau3D<%.4f", trueMin, trueMax));
    ds_red_np_mc->SetName("ds_red_np_mc");

    ws_true->import(*ds_red_np_mc);

    // auto ctau3D = (RooRealVar *)ds_red_np_mc->get(0)->find("ctau3D");
    // ctau3D->Print("V");

    double ctau3DMin = 1e-6,
           ctau3DMax = 10;
    ws_true->var("ctau3D")->setRange(ctau3DMin, ctau3DMax);
    ws_true->var("ctau3D")->setRange("trueFitWindow", ctau3DMin, 8);
    ws_true->var("ctau3D")->Print();

    //***********************************************************************
    //**************************** CTAU TRUE FIT ****************************
    //***********************************************************************
    cout << "\n************ Start MC Ctau True Fit ***************\n\n";

    // mc resolution model
    RooRealVar true_bias1("true_bias1", "bias of detector resolution 1", 1e-4);
    RooRealVar true_sigma1("true_sigma1", "", 0.002, 0.0003, 1);
    RooGaussModel true_gm1("true_gm1", "gauss model 1", *ws_true->var("ctau3D"), RooConst(0), true_sigma1);

    RooRealVar true_sigma2("true_sigma2", "", 0.002, 0.0003, 1);
    RooGaussModel true_gm2("true_gm2", "gauss model 2", *ws_true->var("ctau3D"), true_bias1, true_sigma2);

    RooRealVar true_gw1("true_gw1", "Weight for decay1", 0.5, 0.0, 1.0); // 가중치 w1
    RooFormulaVar true_gw2("true_gw2", "1 - true_gw1", RooArgList(true_gw1));
    RooAddModel true_gmsum("true_gmsum", "sum of true_gm1 and true_gm2", RooArgList(true_gm1, true_gm1), RooArgList(true_gw1));

    RooTruthModel delta_fcn("RooTruthModel", "", *ws_true->var("ctau3D"));


    // ===== Build model ===== //
    // Decay (x) res model
    RooRealVar true_tau1("true_tau1", "lifetime", 0.02, 0.001, 1);
    RooRealVar true_tau2("true_tau2", "lifetime", 0.002, 0.0001, 1);
    RooRealVar true_tau3("true_tau3", "lifetime", 0.002, 0.0001, 1);
    RooDecay true_decay_r1("true_decay_r1", "decay model", *ws_true->var("ctau3D"), true_tau1, delta_fcn, RooDecay::SingleSided);
    RooDecay true_decay_r2("true_decay_r2", "decay model", *ws_true->var("ctau3D"), true_tau2, delta_fcn, RooDecay::SingleSided);
    //      RooDecay decay_r3("decay_r3", "decay model", *ws_true->var("ctau3D"), tau3, true_gm1, RooDecay::SingleSided);

    RooRealVar true_frac1("true_frac1", "Weight for decay1", 0.5, 0.01, 1.0); // 가중치 true_frac1

    RooAddPdf true_decay("true_decay", "", RooArgList(true_decay_r1, true_decay_r2), RooArgList(true_frac1));
    //      RooAddPdf decay("decay", "", RooArgList(decay_r1, true_decay_r2, decay_r3), RooArgList(true_frac1, w2, w3));

    // extend decay model
    RooRealVar true_n_dimu("true_n_dimu", "", 500000, 200000, 1000000);
    auto fit_model = new RooExtendPdf("fit_model", "", true_decay, true_n_dimu);
    ws_true->import(*fit_model);

    // ===== Fit here ===== //
    bool isWeighted = ds_red_np_mc->isWeighted();

    // auto fit_true = decay.fitTo(*ds_red_np_mc, Save(1), SumW2Error(isWeighted), NumCPU(6), Extended(1), Range("ctauTrueRange"), PrintLevel(-1));
    auto fit_true = ws_true->pdf("fit_model")->fitTo(*ds_red_np_mc, Save(1), Extended(1), SumW2Error(false), NumCPU(12), Range("trueFitWindow"), Strategy(2), RecoverFromUndefinedRegions(1.), PrintEvalErrors(-1));
    fit_true->Print("v");

    // set constant
    true_bias1.setConstant(true);
    true_sigma1.setConstant(true);
    true_sigma2.setConstant(true);
    true_gw1.setConstant(true);
    true_tau1.setConstant(true);
    true_tau2.setConstant(true);
    true_tau3.setConstant(true);
    true_frac1.setConstant(true);


    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Plotting - very long but simple *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    // build main canvas
    TCanvas *c_true = new TCanvas("c_true", "", 800, 800);
    c_true->cd();
    TPad *pad_true = new TPad("pad_true", "pad_true", 0, 0.25, 0.98, 1.0);
    pad_true->SetTicks(1, 1);
    pad_true->Draw();
    pad_true->cd();
    gPad->SetLogy();


    // ===== ctau3DTrue frame ==== //
    RooPlot *true_frame = ws_true->var("ctau3D")->frame(Bins(nCtauBins),Range(trueMin, trueMax));
    // fit_model.setNormRange("ctauTrueRange");
    ds_red_np_mc->plotOn(true_frame, DataError(RooAbsData::SumW2), Name("ds_red_np_mc"));
    ws_true->pdf("fit_model")->plotOn(true_frame, Range("trueFitWindow"), NormRange("trueFitWindow"), Name("fit_model"));

    // set y plot range
    // RooPlot *true_frame = mass->frame(nMassBin); // bins
    TH1 *h_tmp = ds_red_np_mc->createHistogram("hist", *ws_true->var("ctau3D"), Binning(true_frame->GetNbinsX(), true_frame->GetXaxis()->GetXmin(), true_frame->GetXaxis()->GetXmax()));
    Double_t YMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
    Double_t YMin = 1e99;
    for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
        if (h_tmp->GetBinContent(i) > 0)
            YMin = min(YMin, h_tmp->GetBinContent(i));
    double Ydown;
    double Yup;
    Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
    Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));
    true_frame->GetYaxis()->SetRangeUser(Ydown, Yup);

    // some cosmetics
    true_frame->SetFillStyle(4000);
    true_frame->GetYaxis()->SetTitleOffset(1.43);
    true_frame->GetXaxis()->SetLabelSize(0); // to hide x-axis label
    true_frame->GetXaxis()->SetTitleSize(0);
    // true_frame->SetMaximum(1e7);
    // true_frame->SetMinimum(1e-1);
    true_frame->Draw();

    
    // ===== draw pull dist. ===== //
    TPad *pull_pad = new TPad("pull_pad", "pull_pad", 0, 0.001, 0.98, 0.32);
    c_true->cd();
    pull_pad->Draw();
    pull_pad->cd();
    pull_pad->SetTopMargin(0); // Upper and lower plot are joined
    pull_pad->SetBottomMargin(0.67);
    pull_pad->SetBottomMargin(0.4);
    pull_pad->SetFillStyle(4000);
    pull_pad->SetFrameFillStyle(4000);
    pull_pad->SetTicks(1, 1);

    // to ignore empty residual bin warning - it's not a matter.
    RooMsgService::instance().getStream(0).removeTopic(Plotting);
    RooMsgService::instance().getStream(1).removeTopic(Plotting);

    RooPlot *tmp_pull_frame = (RooPlot *)true_frame->Clone("tmp_pull_frame");
    RooHist *pull_hist = tmp_pull_frame->pullHist("ds_red_np_mc", "fit_model", true);
    // RooHist *pull_hist = tmp_pull_frame->pullHist("data_ctauBkg", "ctauBkg_Tot", true); // flow code
    pull_hist->SetMarkerSize(0.8);
    RooPlot *pull_frame = ws_true->var("ctau3D")->frame(Title("Pull Distribution"), Bins(nCtauBins), Range(trueMin, trueMax));
    pull_frame->addPlotable(pull_hist, "PX");
    
    pull_frame->SetTitle("");
    pull_frame->SetTitleSize(0);
    pull_frame->GetYaxis()->SetTitleOffset(0.3);
    pull_frame->GetYaxis()->SetTitle("Pull");
    pull_frame->GetYaxis()->SetTitleSize(0.08);
    pull_frame->GetYaxis()->SetLabelSize(0.08);
    pull_frame->GetYaxis()->SetRangeUser(-9.5, 9.5);
    pull_frame->GetYaxis()->CenterTitle();
    pull_frame->GetYaxis()->SetTickSize(0.04);
    pull_frame->GetYaxis()->SetNdivisions(404);
    pull_frame->GetXaxis()->SetTickSize(0.03);

    pull_frame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
    pull_frame->GetXaxis()->SetTitleOffset(1.55);
    pull_frame->GetXaxis()->SetLabelOffset(0.04);
    pull_frame->GetXaxis()->SetLabelSize(0.08);
    pull_frame->GetXaxis()->SetTitleSize(0.08);
    pull_frame->GetXaxis()->CenterTitle();
    pull_frame->Draw();
    

    // ===== horizontoal lines on pull hist ===== //
    // maybe we can make function for them...
    TLine *line_at_2 = new TLine(pull_frame->GetXaxis()->GetXmin(), 2, pull_frame->GetXaxis()->GetXmax(), 2);
    line_at_2->SetLineColor(kBlack);
    line_at_2->SetLineStyle(3);
    line_at_2->SetLineWidth(1);
    line_at_2->Draw();

    TLine *line_at_4 = new TLine(pull_frame->GetXaxis()->GetXmin(), 4, pull_frame->GetXaxis()->GetXmax(), 4);
    line_at_4->SetLineColor(kBlack);
    line_at_4->SetLineStyle(3);
    line_at_4->SetLineWidth(2);
    line_at_4->Draw();

    TLine *line_at_8 = new TLine(pull_frame->GetXaxis()->GetXmin(), 8, pull_frame->GetXaxis()->GetXmax(), 8);
    line_at_8->SetLineColor(kBlack);
    line_at_8->SetLineStyle(3);
    line_at_8->SetLineWidth(1);
    line_at_8->Draw();

    // negative
    TLine *line_at_n2 = new TLine(pull_frame->GetXaxis()->GetXmin(), -2, pull_frame->GetXaxis()->GetXmax(), -2);
    line_at_n2->SetLineColor(kBlack);
    line_at_n2->SetLineStyle(3);
    line_at_n2->SetLineWidth(1);
    line_at_n2->Draw();

    TLine *line_at_n4 = new TLine(pull_frame->GetXaxis()->GetXmin(), -4, pull_frame->GetXaxis()->GetXmax(), -4);
    line_at_n4->SetLineColor(kBlack);
    line_at_n4->SetLineStyle(3);
    line_at_n4->SetLineWidth(2);
    line_at_n4->Draw();

    TLine *line_at_n8 = new TLine(pull_frame->GetXaxis()->GetXmin(), -8, pull_frame->GetXaxis()->GetXmax(), -8);
    line_at_n8->SetLineColor(kBlack);
    line_at_n8->SetLineStyle(3);
    line_at_n8->SetLineWidth(1);
    line_at_n8->Draw();


    // ===== chi2 ===== //
    printChi2(ws_true, pull_pad, tmp_pull_frame, fit_true, "ctau3D", "ds_red_np_mc", "fit_model", nCtauBins, false);
    pull_pad->Update();


    // ===== draw main canvas ===== //
    c_true->Draw();
    c_true->SaveAs(("figs/" + out_ss + ".png").c_str());


    // ===== Export results ===== //
    auto out_file = new TFile(("roots/" + out_ss + ".root").c_str(), "recreate");

    ws_true->Write();
    out_file->Close();

    cout << "\n=================================\n";
    cout << "\n Finish ctau bkg fit\n";
    cout << "\n=================================\n";
    t->Stop();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());

    ws_true->Print("V");
    ws_true->var("true_tau1")->Print("V");
    cout << ws_true->var("true_tau1")->getVal() << endl;
}


void not_used_test_code() // please ignore this part! - pjgwak
{
    // ctau3D->setRange(0, 0.2);
    // auto c2 = new TCanvas("c2", "", 800, 800);
    // RooPlot *ctau_frame_temp = ctau3D->frame();
    // ds_red_np_mc->plotOn(ctau_frame_temp, DataError(RooAbsData::SumW2));
    // gPad->SetLogy();
    // ctau_frame_temp->SetMinimum(1e-6);
    // ctau_frame_temp->Draw();
    // c2->Draw();
}