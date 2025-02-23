#include <iostream>
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include "RooMsgService.h"
#include "RooWorkspace.h"
#include "RooFit.h"
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../headers/rootFitHeaders.h"
#include "../headers/commonUtility.h"
#include "../headers/JpsiUtility.h"
#include "../headers/cutsAndBin.h"
#include "../headers/CMS_lumi_v2mass.C"
#include "../headers/tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

// prompt MC fit
void mc_mass()
{
    TStopwatch *t = new TStopwatch;
    t->Start();
    cout << "\n=================================\n";
    cout << "\n Start MC mass fit\n";
    cout << "\n=================================\n";

    using namespace RooFit;

    // mc mass fit range - You can change here! (MC mass only)
    float mass_min = 2.6;
    float mass_max = 3.25;
    // mass range max:
    // 3.34 for forward
    // 3.28 for mid , 3.19837 (3.3까지 그리고 그림 확대해서 데이터랑 pdf 엇나가는 부분을 눈으로 보고 잘랐다.) -> AN21-003 Fig.21을 봐도 Pull이 6.8까지 튀어도 일단 넘어갈 수 있다. 풀 전체가 요동치지만 않으면

    // input parameters
    float cos_low = 0.5, cos_high = 0.6;
    float ptLow = 3;
    float ptHigh = 6.5;
    float yLow = 1.6;
    float yHigh = 2.4;
    int cLow = 60;
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
    RooMsgService::instance().getStream(0).removeTopic(Caching); // in RooFit::
    RooMsgService::instance().getStream(1).removeTopic(Caching);
    RooMsgService::instance().getStream(0).removeTopic(Plotting);
    RooMsgService::instance().getStream(1).removeTopic(Plotting);
    RooMsgService::instance().getStream(0).removeTopic(Integration);
    RooMsgService::instance().getStream(1).removeTopic(Integration);
    RooMsgService::instance().getStream(0).removeTopic(InputArguments);
    RooMsgService::instance().getStream(1).removeTopic(InputArguments);
    RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    RooMsgService::instance().getStream(1).removeTopic(Minimization);
    RooMsgService::instance().setGlobalKillBelow(WARNING);


    // ===== labeling ===== //
    std::ostringstream out_file_path;
    out_file_path << "mc_mass/mc_mass_pt" << ptLow << "-" << ptHigh
                  << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
                  << "_cos" << cos_low << "_" << cos_high;
    string out_ss = out_file_path.str();


    // ===== make output folders ===== //
    gSystem->mkdir(Form("roots/mc_mass/"), kTRUE);
    gSystem->mkdir(Form("figs/mc_mass"), kTRUE);


    // ===== kinematic cut ===== //
    TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);
    TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
    TString OS = "recoQQsign==0 &&";

    TString angle_cut = Form("&& cos_ep>%.2f && cos_ep<%.2f", cos_low, cos_high);

    TString final_cut = OS + accCut + kineCut + angle_cut;


    // ===== import inputs ===== //
    auto f_pr_mc = new TFile("../files_roodata/RooDataSet_miniAOD_isMC1_PR_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root");
    f_pr_mc->SetName("f_pr_mc");
    // f_pr_mc.Print("V");
    RooDataSet *ds_mc_temp = (RooDataSet *)f_pr_mc->Get("dataset");
    ds_mc_temp->SetName("ds_mc_temp");

    // apply kinematic cuts
    auto ds_red_pr_mc = (RooDataSet *)ds_mc_temp->reduce(final_cut.Data());
    ds_red_pr_mc->SetName("ds_red_pr_mc");
    // ds_red_pr_mc->Print("V");

    // set mass range
    RooRealVar *mc_mass = (RooRealVar *)ds_red_pr_mc->get(0)->find("mass");
    mc_mass->setRange(massLow, massHigh);
    mc_mass->setRange("mcMassFit", mass_min, mass_max);

    // // make a workspace
    // auto ws = new RooWorkspace("ws");
    // ws->import(*ds_red_pr_mc);


    // ===== define mc fit model ===== //
    // signal
    auto mass_mean = new RooRealVar("mass_mean", "", 3.096, 3.086, 3.106);
    auto mass_sigma1 = new RooRealVar("mass_sigma1", "", 0.01, 0.001, 0.1);
    auto mass_sigma2 = new RooRealVar("mass_sigma2", "", 0.03, 0.001, 0.1);
    auto mass_alpha1 = new RooRealVar("mass_alpha1", "", 2, 1, 5);
    auto mass_power1 = new RooRealVar("mass_power1", "", 2, 1, 5);
    auto mass_fracG1 = new RooRealVar("mass_fracG1", "", 0.6, 0.05, 0.95);

    auto mc_G1Sig = new RooGaussian("mc_G1Sig", "", *mc_mass, *mass_mean, *mass_sigma1);
    auto mc_CB1Sig = new RooCBShape("mc_CB1Sig", "", *mc_mass, *mass_mean, *mass_sigma2, *mass_alpha1, *mass_power1);
    auto mc_G1CB1Sig = new RooAddPdf("mc_G1CB1Sig", "", RooArgList(*mc_G1Sig, *mc_CB1Sig), RooArgList(*mass_fracG1));

    // mc doesn't use bkg
    // ws->factory("Exponential::expBkg(mass,coefExp[-1.,-3.,1.])");

    auto mc_n_sig = new RooRealVar("mc_n_sig", "", 20000, 1, 50000);
    auto mc_mass_pdf = new RooAddPdf("mc_mass_pdf", "", *mc_G1CB1Sig, *mc_n_sig);


    // ======= Fit here ===== //
    bool is_weighted = ds_red_pr_mc->isWeighted();
    cout << "Is weighted: " << is_weighted << endl;

    auto fit_mc_mass = mc_mass_pdf->fitTo(*ds_red_pr_mc, Extended(1), Save(1), SumW2Error(false), NumCPU(6), Strategy(2), Range("mcMassFit"), Timer(1));
    // PrintLevel(-1)
    fit_mc_mass->Print("V");


    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Plotting - very long but simple *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    TCanvas *c_mc_mass = new TCanvas("c_mc_mass", "", 800, 800);
    c_mc_mass->cd();
    TPad *pad_mass = new TPad("pad_mass", "pad_mass", 0, 0.25, 0.98, 1.0);
    pad_mass->SetTicks(1, 1);
    pad_mass->Draw();
    pad_mass->cd();
    gPad->SetLogy();


    // ===== draw mass dist. ===== //
    auto mc_mass_frame = (RooPlot *)mc_mass->frame(Range(massLow, massHigh), Binning(nMassBin));
    ds_red_pr_mc->plotOn(mc_mass_frame, Name("ds_red_pr_mc"), DataError(RooAbsData::SumW2), MarkerSize(.7));
    mc_mass_pdf->plotOn(mc_mass_frame, Name("mc_mass_pdf"), NormRange("mcMassFit"), LineColor(kBlack));
    mc_mass_pdf->plotOn(mc_mass_frame, LineStyle(kDashed), Components(RooArgSet(*mc_CB1Sig)), Name("mc_CB1Sig"), LineColor(44), LineWidth(2));
    mc_mass_pdf->plotOn(mc_mass_frame, LineStyle(kDashed), Components(RooArgSet(*mc_G1Sig)), Name("mc_G1Sig"), LineColor(8), LineWidth(2));

    // set y plot range
    // RooPlot *mc_mass_frame = mc_mass->frame(nMassBin); // bins
    TH1 *h_tmp = ds_red_pr_mc->createHistogram("hist", *mc_mass, Binning(mc_mass_frame->GetNbinsX(), mc_mass_frame->GetXaxis()->GetXmin(), mc_mass_frame->GetXaxis()->GetXmax()));
    Double_t YMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
    Double_t YMin = 1e99;
    for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
        if (h_tmp->GetBinContent(i) > 0)
            YMin = min(YMin, h_tmp->GetBinContent(i));
    double Ydown;
    double Yup;
    Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
    Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));
    // mc_mass_frame->SetMaximum(1e6);
    // mc_mass_frame->SetMinimum(20); // yLog can't get 0 as minimum
    mc_mass_frame->GetYaxis()->SetRangeUser(Ydown, Yup);

    mc_mass_frame->SetFillStyle(4000);
    mc_mass_frame->GetYaxis()->SetTitleOffset(1.43);
    mc_mass_frame->GetXaxis()->SetLabelSize(0); // to hide x-axis label
    mc_mass_frame->GetXaxis()->SetTitleSize(0);
    // mc_mass_frame->GetXaxis()->CenterTitle();
    // mc_mass_frame->GetXaxis()->SetRangeUser(massLow, massHigh);
    // mc_mass_frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    mc_mass_frame->Draw();


    // ===== draw legends for PDFs and dataset ===== //
    // mc_mass_pdf->plotOn(mc_mass_frame, Name("pdfMASS_Tot"), LineColor(kBlack));
    TLegend *leg_pdfs = new TLegend(text_x + 0.45, text_y - 0.2, text_x + 0.65, text_y);
    leg_pdfs->SetTextSize(text_size);
    leg_pdfs->SetTextFont(43);
    leg_pdfs->SetBorderSize(0);
    leg_pdfs->AddEntry(mc_mass_frame->findObject("ds_red_pr_mc"), "MC", "pe");
    leg_pdfs->AddEntry(mc_mass_frame->findObject("mc_mass_pdf"), "Total", "l");
    leg_pdfs->AddEntry(mc_mass_frame->findObject("mc_CB1Sig"), "CB1", "l");
    leg_pdfs->AddEntry(mc_mass_frame->findObject("mc_G1Sig"), "Gauss1", "l");
    leg_pdfs->Draw("same");
    

    // ===== draw legends ===== //
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c ; Cent. %d - %d%s; %.2f < cos#theta_{EP} < %.2f", ptLow, ptHigh, cLow / 2, cHigh / 2, "%", cos_low, cos_high), text_x, text_y, text_color, text_size);
    if (yLow == 0)
        drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff, text_color, text_size);
    else if (yLow != 0)
        drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff, text_color, text_size);
    // drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
    drawText(Form("N_{J/#psi} = %.f #pm %.f", mc_n_sig->getVal(), mc_n_sig->getError()), text_x, text_y - y_diff * 2, text_color, text_size);
    drawText(Form("m_{J/#psi} = %.4f #pm %.4f", mass_mean->getVal(), mass_mean->getError()), text_x, text_y - y_diff * 3, text_color, text_size);
    drawText(Form("#alpha_{J/#psi} = %.4f #pm %.4f", mass_alpha1->getVal(), mass_alpha1->getError()), text_x, text_y - y_diff * 4, text_color, text_size);
    drawText(Form("f_{G1} = %.4f #pm %.4f", mass_fracG1->getVal(), mass_fracG1->getError()), text_x, text_y - y_diff * 5, text_color, text_size);
    drawText(Form("n_{J/#psi} = %.4f #pm %.4f", mass_power1->getVal(), mass_power1->getError()), text_x, text_y - y_diff * 6, text_color, text_size);
    drawText(Form("#sigma1_{J/#psi} = %.2f #pm %.2f MeV/c^{2}", (mass_sigma1->getVal()) * 1000, (mass_sigma1->getError()) * 1000), text_x, text_y - y_diff * 7, text_color, text_size);
    drawText(Form("#sigma2_{J/#psi} = %.2f #pm %.2f MeV/c^{2}", (mass_sigma2->getVal()) * 1000, (mass_sigma2->getError()) * 1000), text_x, text_y - y_diff * 8, text_color, text_size);

    // ===== draw pull dist. ===== //
    TPad *pad_mass_pull = new TPad("pad_mass_pull", "pad_mass_pull", 0, 0.001, 0.98, 0.32);
    c_mc_mass->cd();
    pad_mass_pull->Draw();
    pad_mass_pull->cd();
    pad_mass_pull->SetTopMargin(0); // Upper and lower plot are joined
    pad_mass_pull->SetBottomMargin(0.67);
    pad_mass_pull->SetBottomMargin(0.4);
    pad_mass_pull->SetFillStyle(4000);
    pad_mass_pull->SetFrameFillStyle(4000);
    pad_mass_pull->SetTicks(1, 1);

    auto mc_mass_pull = mc_mass_frame->pullHist("ds_red_pr_mc", "mc_mass_pdf", true);
    auto mc_mass_pull_frame = mc_mass->frame(Title(";;Pull Distribution"), Range(2.6, 3.5));
    mc_mass_pull->SetMarkerSize(0.8);
    mc_mass_pull_frame->addPlotable(mc_mass_pull, "P");

    // cosmetics
    mc_mass_pull_frame->SetTitle("");
    mc_mass_pull_frame->SetTitleSize(0);
    mc_mass_pull_frame->GetYaxis()->SetTitleOffset(0.3);
    mc_mass_pull_frame->GetYaxis()->SetTitle("Pull");
    mc_mass_pull_frame->GetYaxis()->SetTitleSize(0.08);
    mc_mass_pull_frame->GetYaxis()->SetLabelSize(0.08);
    mc_mass_pull_frame->GetYaxis()->SetRangeUser(-9, 9);
    mc_mass_pull_frame->GetYaxis()->CenterTitle();

    mc_mass_pull_frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    mc_mass_pull_frame->GetXaxis()->SetTitleOffset(1.55);
    mc_mass_pull_frame->GetXaxis()->SetLabelOffset(0.04);
    mc_mass_pull_frame->GetXaxis()->SetLabelSize(0.08);
    mc_mass_pull_frame->GetXaxis()->SetTitleSize(0.08);
    mc_mass_pull_frame->GetXaxis()->CenterTitle();

    mc_mass_pull_frame->GetYaxis()->SetTickSize(0.04);
    mc_mass_pull_frame->GetYaxis()->SetNdivisions(404);
    mc_mass_pull_frame->GetXaxis()->SetTickSize(0.03);
    mc_mass_pull_frame->Draw();


    // ===== horizontoal lines on pull hist ===== //
    // maybe we can make function for them...
    TLine *line_at_2 = new TLine(mc_mass_pull_frame->GetXaxis()->GetXmin(), 2, mc_mass_pull_frame->GetXaxis()->GetXmax(), 2);
    line_at_2->SetLineColor(kBlack);
    line_at_2->SetLineStyle(3);
    line_at_2->SetLineWidth(1);
    line_at_2->Draw();

    TLine *line_at_4 = new TLine(mc_mass_pull_frame->GetXaxis()->GetXmin(), 4, mc_mass_pull_frame->GetXaxis()->GetXmax(), 4);
    line_at_4->SetLineColor(kBlack);
    line_at_4->SetLineStyle(3);
    line_at_4->SetLineWidth(2);
    line_at_4->Draw();

    TLine *line_at_8 = new TLine(mc_mass_pull_frame->GetXaxis()->GetXmin(), 8, mc_mass_pull_frame->GetXaxis()->GetXmax(), 8);
    line_at_8->SetLineColor(kBlack);
    line_at_8->SetLineStyle(3);
    line_at_8->SetLineWidth(1);
    line_at_8->Draw();

    // negative
    TLine *line_at_n2 = new TLine(mc_mass_pull_frame->GetXaxis()->GetXmin(), -2, mc_mass_pull_frame->GetXaxis()->GetXmax(), -2);
    line_at_n2->SetLineColor(kBlack);
    line_at_n2->SetLineStyle(3);
    line_at_n2->SetLineWidth(1);
    line_at_n2->Draw();

    TLine *line_at_n4 = new TLine(mc_mass_pull_frame->GetXaxis()->GetXmin(), -4, mc_mass_pull_frame->GetXaxis()->GetXmax(), -4);
    line_at_n4->SetLineColor(kBlack);
    line_at_n4->SetLineStyle(3);
    line_at_n4->SetLineWidth(2);
    line_at_n4->Draw();

    TLine *line_at_n8 = new TLine(mc_mass_pull_frame->GetXaxis()->GetXmin(), -8, mc_mass_pull_frame->GetXaxis()->GetXmax(), -8);
    line_at_n8->SetLineColor(kBlack);
    line_at_n8->SetLineStyle(3);
    line_at_n8->SetLineWidth(1);
    line_at_n8->Draw();

    c_mc_mass->Draw();
    c_mc_mass->SaveAs(("figs/"+out_ss+".png").c_str());


    // ===== Export results ===== //
    auto out_file = new TFile(("roots/"+out_ss+".root").c_str(), "recreate");

    auto ws_mc_mass = new RooWorkspace("ws_mc_mass");
    ws_mc_mass->import(*mc_G1CB1Sig);

    c_mc_mass->Write();
    fit_mc_mass->Write();
    ws_mc_mass->Write();
    out_file->Close();

    cout << "\n=================================\n";
    cout << "\n Finish MC mass fit\n";
    cout << "\n=================================\n";
    t->Stop();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
}