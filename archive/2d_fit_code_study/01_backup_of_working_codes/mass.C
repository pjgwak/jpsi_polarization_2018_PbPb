#include "RooFit.h"
#include "../headers/rootFitHeaders.h"
#include "../headers/commonUtility.h"
#include "../headers/JpsiUtility.h"
#include "../headers/cutsAndBin.h"
#include "../headers/CMS_lumi_v2mass.C"
#include "../headers/tdrstyle.C"

void mass()
{
    TStopwatch *t = new TStopwatch;
    t->Start();
    cout << "\n=================================\n";
    cout << "\n Start data mass fit\n";
    cout << "\n=================================\n";

    using namespace RooFit;

    // input parameters
    float cos_low = 0, cos_high = 0.1;
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
    RooMsgService::instance().getStream(0).removeTopic(Caching); // in RooFit:: namespace
    RooMsgService::instance().getStream(1).removeTopic(Caching);
    RooMsgService::instance().getStream(0).removeTopic(Plotting);
    RooMsgService::instance().getStream(1).removeTopic(Plotting);
    RooMsgService::instance().getStream(0).removeTopic(InputArguments);
    RooMsgService::instance().getStream(1).removeTopic(InputArguments);
    RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    RooMsgService::instance().getStream(1).removeTopic(Minimization);
    RooMsgService::instance().setGlobalKillBelow(WARNING);


    // ===== labeling ===== //
    std::ostringstream out_file_path;
    out_file_path << "mass/mass_pt" << ptLow << "-" << ptHigh
                  << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
                  << "_cos" << cos_low << "_" << cos_high;
    string out_ss = out_file_path.str();


    // ===== make output folders ===== //
    gSystem->mkdir(Form("roots/mass/"), kTRUE);
    gSystem->mkdir(Form("figs/mass"), kTRUE);


    // ===== kinematic cut ===== //
    TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);
    TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
    TString OS = "recoQQsign==0 &&";

    TString angle_cut = Form("&& cos_ep>%.2f && cos_ep<%.2f", cos_low, cos_high);

    TString final_cut = OS + accCut + kineCut + angle_cut;


    // ===== import inputs ===== //
    // mass data
    auto f_data = new TFile("../files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root");
    f_data->SetName("f_data");
    // f_pr_mc.Print("V");

    // make a reduced dataset
    RooDataSet *ds_tmp = (RooDataSet *)f_data->Get("dataset");
    ds_tmp->SetName("ds_tmp");
    RooDataSet *ds_tmp_weight = new RooDataSet("ds_tmp_weight", "", *ds_tmp->get(), Import(*ds_tmp), WeightVar("weight"));

    auto ds_red_mass = (RooDataSet *)ds_tmp_weight->reduce(final_cut.Data());
    ds_red_mass->SetName("ds_red_mass");
    // ds_red_mass->Print("V");
    // ws_mass->import(*ds_red_mass);

    auto ws_mass = new RooWorkspace("ws_mass");

    // set mass range
    RooRealVar *mass = (RooRealVar *)ds_red_mass->get(0)->find("mass");
    mass->setRange(massLow, massHigh);

    // mc mass fit result
    std::ostringstream in_mc_fit;
    in_mc_fit << "roots/mc_mass/mc_mass_pt" << ptLow << "-" << ptHigh
              << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
              << "_cos" << cos_low << "_" << cos_high
              << ".root";
    string ss = in_mc_fit.str();
    auto f_mc_fit = new TFile(ss.c_str());
    auto ws_mc_mass = (RooWorkspace *)f_mc_fit->Get("ws_mc_mass");


    // ===== build model ===== //
    // signal
    auto mass_mean = (RooRealVar *)ws_mc_mass->var("mass_mean")->Clone("mass_mean");
    auto mass_sigma1 = (RooRealVar *)ws_mc_mass->var("mass_sigma1")->Clone("mass_sigma1");
    auto mass_sigma2 = (RooRealVar *)ws_mc_mass->var("mass_sigma2")->Clone("mass_sigma2");
    auto mass_alpha1 = (RooRealVar *)ws_mc_mass->var("mass_alpha1")->Clone("mass_alpha1");
    auto mass_power1 = (RooRealVar *)ws_mc_mass->var("mass_power1")->Clone("mass_power1");
    auto mass_fracG1 = (RooRealVar *)ws_mc_mass->var("mass_fracG1")->Clone("mass_fracG1");

    mass_alpha1->setConstant(kTRUE);
    mass_power1->setConstant(kTRUE);
    mass_fracG1->setConstant(kTRUE);

    // manual adjustment - when you need
    // mass_mean->setRange(3.07, 3.12);
    // mass_sigma1->setRange(0.001, 0.2);

    auto G1Sig = new RooGaussian("G1Sig", "", *mass, *mass_mean, *mass_sigma1);
    auto CB1Sig = new RooCBShape("CB1Sig", "", *mass, *mass_mean, *mass_sigma2, *mass_alpha1, *mass_power1);
    auto mass_sig = new RooAddPdf("mass_sig", "", RooArgList(*G1Sig, *CB1Sig), RooArgList(*mass_fracG1));

    // bkg
    RooRealVar *mass_sl1 = new RooRealVar("mass_sl1", "sl1", 0.1, -1., 1.);
    RooRealVar *mass_sl2 = new RooRealVar("mass_sl2", "sl2", 0.01, -1., 1.);
    auto mass_bkg = new RooChebychev("mass_bkg", "Background", *mass, RooArgList(*mass_sl1, *mass_sl2));

    // total
    auto n_sig = new RooRealVar("n_sig", "", 5000, 1, 25000);
    auto n_bkg = new RooRealVar("n_bkg", "", 20000, 1, 100000);
    auto mass_pdf = new RooAddPdf("mass_pdf", "", RooArgList(*mass_sig, *mass_bkg), RooArgList(*n_sig, *n_bkg));
    ws_mass->import(*mass_pdf);

    // ===== fit here ===== //
    bool is_weighted = ds_red_mass->isWeighted();
    cout << "Is weighted: " << is_weighted << endl;

    auto fit_mass = ws_mass->pdf("mass_pdf")->fitTo(*ds_red_mass, Extended(1), Save(1), AsymptoticError(is_weighted), NumCPU(12), Strategy(2));
    // PrintLevel(-1), AsymptoticError(), SumW2Error();
    fit_mass->Print("v");

    // fix parameters for next steps (err, res ...)
    ws_mass->var("mass_mean")->setConstant(kTRUE);
    ws_mass->var("mass_sigma1")->setConstant(kTRUE);
    ws_mass->var("mass_sigma2")->setConstant(kTRUE);
    ws_mass->var("mass_sl1")->setConstant(kTRUE);
    ws_mass->var("mass_sl2")->setConstant(kTRUE);


    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Plotting - very long but simple *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    TCanvas *c_mass = new TCanvas("c_mass", "", 800, 800);
    c_mass->cd();
    TPad *pad_mass = new TPad("pad_mass", "pad_mass", 0, 0.25, 0.98, 1.0);
    pad_mass->SetTicks(1, 1);
    pad_mass->Draw();
    pad_mass->cd();
    gPad->SetLogy();

    // mass dist.
    auto mass_frame = (RooPlot *)mass->frame(Range(2.6, 3.5));
    ds_red_mass->plotOn(mass_frame, Name("ds_red_mass"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7));
    ws_mass->pdf("mass_pdf")->plotOn(mass_frame, Name("mass_pdf"), LineColor(kBlack));
    ws_mass->pdf("mass_pdf")->plotOn(mass_frame, LineStyle(kDashed), Components(RooArgSet(*mass_bkg, *CB1Sig)), Name("CB1Sig"), LineColor(44), LineWidth(2));
    ws_mass->pdf("mass_pdf")->plotOn(mass_frame, LineStyle(kDashed), Components(RooArgSet(*mass_bkg, *G1Sig)), Name("G1Sig"), LineColor(8), LineWidth(2));
    ws_mass->pdf("mass_pdf")->plotOn(mass_frame, LineStyle(kDashed), Components(RooArgSet(*mass_bkg)), Name("mass_bkg"), LineColor(kBlue+2), LineWidth(2));

    // set y plot range
    // RooPlot *mass_frame = mass->frame(nMassBin); // bins
    TH1 *h_tmp = ds_red_mass->createHistogram("hist", *mass, Binning(mass_frame->GetNbinsX(), mass_frame->GetXaxis()->GetXmin(), mass_frame->GetXaxis()->GetXmax()));
    Double_t YMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
    Double_t YMin = 1e99;
    for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
        if (h_tmp->GetBinContent(i) > 0)
            YMin = min(YMin, h_tmp->GetBinContent(i));
    double Ydown;
    double Yup;
    Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
    Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));
    // mass_frame->SetMaximum(1e6);
    // mass_frame->SetMinimum(20); // yLog can't get 0 as minimum
    mass_frame->GetYaxis()->SetRangeUser(Ydown, Yup);

    mass_frame->SetFillStyle(4000);
    mass_frame->GetYaxis()->SetTitleOffset(1.43);
    mass_frame->GetXaxis()->SetLabelSize(0); // to hide x-axis label
    mass_frame->GetXaxis()->SetTitleSize(0);
    // mass_frame->GetXaxis()->CenterTitle();
    // mass_frame->GetXaxis()->SetRangeUser(massLow, massHigh);
    // mass_frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    mass_frame->Draw();


    // ===== draw legends for PDFs and dataset ===== //
    TLegend *leg_pdfs = new TLegend(text_x + 0.45, text_y - 0.2, text_x + 0.65, text_y);
    leg_pdfs->SetTextSize(text_size);
    leg_pdfs->SetTextFont(43);
    leg_pdfs->SetBorderSize(0);
    leg_pdfs->AddEntry(mass_frame->findObject("ds_red_mass"), "Data", "pe");
    leg_pdfs->AddEntry(mass_frame->findObject("mass_pdf"), "Total", "l");
    leg_pdfs->AddEntry(mass_frame->findObject("CB1Sig"), "CB1+Bkg", "l");
    leg_pdfs->AddEntry(mass_frame->findObject("G1Sig"), "Gauss1+Bkg", "l");
    leg_pdfs->AddEntry(mass_frame->findObject("mass_bkg"), "2nd Chebychev", "l");
    leg_pdfs->Draw("same");

    // ws_mass->var(;
    
    // ===== draw legends ===== //
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c ; Cent. %d - %d%s; %.2f < cos#theta_{EP} < %.2f", ptLow, ptHigh, cLow / 2, cHigh / 2, "%", cos_low, cos_high), text_x, text_y, text_color, text_size);
    if (yLow == 0)
        drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff, text_color, text_size);
    else if (yLow != 0)
        drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff, text_color, text_size);
    // drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
    drawText(Form("N_{J/#psi} = %.f #pm %.f; N_{Bkg} = %.f #pm %.f", ws_mass->var("n_sig")->getVal(), ws_mass->var("n_sig")->getError(), ws_mass->var("n_bkg")->getVal(), ws_mass->var("n_bkg")->getError()), text_x, text_y - y_diff * 2, text_color, text_size);
    drawText(Form("m_{J/#psi} = %.4f #pm %.4f", ws_mass->var("mass_mean")->getVal(), ws_mass->var("mass_mean")->getError()), text_x, text_y - y_diff * 3, text_color, text_size);
    drawText(Form("#alpha_{J/#psi} = %.4f (fixed); n_{J/#psi} = %.4f (fixed)", ws_mass->var("mass_alpha1")->getVal(), ws_mass->var("mass_power1")->getVal()), text_x, text_y - y_diff * 4, text_color, text_size);
    drawText(Form("f_{G1} = %.4f (fixed)", mass_fracG1->getVal()), text_x, text_y - y_diff * 5, text_color, text_size);
    drawText(Form("#sigma1_{J/#psi} = %.2f #pm %.2f MeV/c^{2}", ws_mass->var("mass_sigma1")->getVal() * 1000, ws_mass->var("mass_sigma1")->getError() * 1000), text_x, text_y - y_diff * 6, text_color, text_size);
    drawText(Form("#sigma2_{J/#psi} = %.2f #pm %.2f MeV/c^{2}", ws_mass->var("mass_sigma2")->getVal() * 1000, ws_mass->var("mass_sigma2")->getError() * 1000), text_x, text_y - y_diff * 7, text_color, text_size);
    drawText(Form("slope1_{Bkg} = %.2f #pm %.2f MeV/c^{2}", (ws_mass->var("mass_sl1")->getVal()), (ws_mass->var("mass_sl1")->getError())), text_x, text_y - y_diff * 8, text_color, text_size);
    drawText(Form("slope2_{Bkg} = %.2f #pm %.2f MeV/c^{2}", (ws_mass->var("mass_sl2")->getVal()), (ws_mass->var("mass_sl2")->getError())), text_x, text_y - y_diff * 9, text_color, text_size);


    // ===== draw pull dist. ===== //
    TPad *pad_mass_pull = new TPad("pad_mass_pull", "pad_mass_pull", 0, 0.001, 0.98, 0.32);
    c_mass->cd();
    pad_mass_pull->Draw();
    pad_mass_pull->cd();
    pad_mass_pull->SetTopMargin(0); // Upper and lower plot are joined
    pad_mass_pull->SetBottomMargin(0.67);
    pad_mass_pull->SetBottomMargin(0.4);
    pad_mass_pull->SetFillStyle(4000);
    pad_mass_pull->SetFrameFillStyle(4000);
    pad_mass_pull->SetTicks(1, 1);

    auto mass_pull = mass_frame->pullHist("ds_red_mass", "mass_pdf", true); // true: isExtended?
    auto mass_pull_frame = mass->frame(Title(";;Pull Distribution"), Range(2.6, 3.5));
    mass_pull->SetMarkerSize(0.8);
    mass_pull_frame->addPlotable(mass_pull, "P");

    // cosmetics
    mass_pull_frame->SetTitle("");
    mass_pull_frame->SetTitleSize(0);
    mass_pull_frame->GetYaxis()->SetTitleOffset(0.3);
    mass_pull_frame->GetYaxis()->SetTitle("Pull");
    mass_pull_frame->GetYaxis()->SetTitleSize(0.08);
    mass_pull_frame->GetYaxis()->SetLabelSize(0.08);
    mass_pull_frame->GetYaxis()->SetRangeUser(-9, 9);
    mass_pull_frame->GetYaxis()->CenterTitle();

    mass_pull_frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    mass_pull_frame->GetXaxis()->SetTitleOffset(1.55);
    mass_pull_frame->GetXaxis()->SetLabelOffset(0.04);
    mass_pull_frame->GetXaxis()->SetLabelSize(0.08);
    mass_pull_frame->GetXaxis()->SetTitleSize(0.08);
    mass_pull_frame->GetXaxis()->CenterTitle();

    mass_pull_frame->GetYaxis()->SetTickSize(0.04);
    mass_pull_frame->GetYaxis()->SetNdivisions(404);
    mass_pull_frame->GetXaxis()->SetTickSize(0.03);
    mass_pull_frame->Draw();


    // ===== horizontoal lines on pull hist ===== //
    // maybe we can make function for them...
    TLine *line_at_2 = new TLine(mass_pull_frame->GetXaxis()->GetXmin(), 2, mass_pull_frame->GetXaxis()->GetXmax(), 2);
    line_at_2->SetLineColor(kBlack);
    line_at_2->SetLineStyle(3);
    line_at_2->SetLineWidth(1);
    line_at_2->Draw();

    TLine *line_at_4 = new TLine(mass_pull_frame->GetXaxis()->GetXmin(), 4, mass_pull_frame->GetXaxis()->GetXmax(), 4);
    line_at_4->SetLineColor(kBlack);
    line_at_4->SetLineStyle(3);
    line_at_4->SetLineWidth(2);
    line_at_4->Draw();

    TLine *line_at_8 = new TLine(mass_pull_frame->GetXaxis()->GetXmin(), 8, mass_pull_frame->GetXaxis()->GetXmax(), 8);
    line_at_8->SetLineColor(kBlack);
    line_at_8->SetLineStyle(3);
    line_at_8->SetLineWidth(1);
    line_at_8->Draw();

    // negative
    TLine *line_at_n2 = new TLine(mass_pull_frame->GetXaxis()->GetXmin(), -2, mass_pull_frame->GetXaxis()->GetXmax(), -2);
    line_at_n2->SetLineColor(kBlack);
    line_at_n2->SetLineStyle(3);
    line_at_n2->SetLineWidth(1);
    line_at_n2->Draw();

    TLine *line_at_n4 = new TLine(mass_pull_frame->GetXaxis()->GetXmin(), -4, mass_pull_frame->GetXaxis()->GetXmax(), -4);
    line_at_n4->SetLineColor(kBlack);
    line_at_n4->SetLineStyle(3);
    line_at_n4->SetLineWidth(2);
    line_at_n4->Draw();

    TLine *line_at_n8 = new TLine(mass_pull_frame->GetXaxis()->GetXmin(), -8, mass_pull_frame->GetXaxis()->GetXmax(), -8);
    line_at_n8->SetLineColor(kBlack);
    line_at_n8->SetLineStyle(3);
    line_at_n8->SetLineWidth(1);
    line_at_n8->Draw();

    c_mass->Draw();
    c_mass->SaveAs(("figs/" + out_ss + ".png").c_str());

    // ===== Export results ===== //    
    auto out_file = new TFile(("roots/" + out_ss + ".root").c_str(), "recreate");

    

    c_mass->Write();
    fit_mass->Write();
    mass_pdf->Write();
    ds_red_mass->Write();
    ws_mass->Write();
    out_file->Close();

    cout << "\n=================================\n";
    cout << "\n Finish data mass fit\n";
    cout << "\n=================================\n";
    t->Stop();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
}