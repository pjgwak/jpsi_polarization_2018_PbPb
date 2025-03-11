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

using namespace std;
using namespace RooFit;

void res()
{
    TStopwatch *t = new TStopwatch;
    t->Start();
    cout << "\n=================================\n";
    cout << "\n Start res fit\n";
    cout << "\n=================================\n";

    // set macro configure
    using namespace std;
    using namespace RooFit;
    gStyle->SetEndErrorSize(0);

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
    out_file_path << "res/res_pt" << ptLow << "-" << ptHigh
                  << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
                  << "_cos" << cos_low << "_" << cos_high;
    string out_ss = out_file_path.str();

    std::ostringstream input_file;
    input_file << "pt" << ptLow << "-" << ptHigh
              << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
              << "_cos" << cos_low << "_" << cos_high
              << ".root";
    string in_ss = input_file.str();
    
    // PR, NP name label
    // TString fname = "";
    // if (PRw == 1) fname = "PR";
    // else if (PRw == 2) fname = "NP";

    // TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

    // ===== make output folders ===== //
    gSystem->mkdir(Form("roots/res/"), kTRUE);
    gSystem->mkdir(Form("figs/res"), kTRUE);


    // ===== kinematic cut ===== //
    // cut is applied in f_err dataset 


    // ===== import inputs ===== //
    auto f_mass = new TFile(("roots/mass/mass_" + in_ss).c_str());
    auto f_err = new TFile(("roots/err/err_" + in_ss).c_str());

    RooAddPdf *mass_pdf = (RooAddPdf *)f_mass->Get("mass_pdf");
    RooDataSet *ds_sig_err = (RooDataSet *)f_err->Get("ds_sig_err");

    RooWorkspace *ws_res = new RooWorkspace("ws_res");
    ws_res->import(*mass_pdf);
    ws_res->import(*ds_sig_err);
    ws_res->var("ctau3DRes")->Print();

 
    // ===== set ranges ??? ===== //
    double ctauErrMin;
    double ctauErrMax;
    ctauErrMin = ws_res->var("ctau3DErr")->getMin();
    ctauErrMax = ws_res->var("ctau3DErr")->getMax();

    // apply Res < 0 cut. Nominal only use negative values
    RooDataSet *ctauResCutDS = (RooDataSet *)ds_sig_err->reduce(RooArgSet(*(ws_res->var("ctau3DRes")), *(ws_res->var("ctau3D")), *(ws_res->var("ctau3DErr"))), Form("ctau3DRes<0&&ctau3DErr>%f&&ctau3DErr<%f", ctauErrMin, ctauErrMax));
    ctauResCutDS->SetName("ctauResCutDS");
    ws_res->import(*ctauResCutDS);
    cout << "\nn_sig: " << ws_res->var("n_sig")->getVal() << "+/-" << ws_res->var("n_sig")->getError() << "\n\n";


    // ===== build model ===== //
    // how many gauss
    int nGauss = 2;

    // set variables
    ws_res->factory("One[1.0]");
    ws_res->factory("ctauRes_mean[0.]"); // resolution bias -> almost 0
    ws_res->factory("ctau1_CtauRes[0.]"); // resolution bias -> almost 0
    ws_res->factory("ctau2_CtauRes[0.]"); // resolution bias -> almost 0
    ws_res->factory("ctau3_CtauRes[0.]"); // resolution bias -> almost 0
    ws_res->factory("ctau4_CtauRes[0.]"); // resolution bias -> almost 0
    
    // sigma of resolutions. 
    //Ratio btw bigger gauss / smaller gaus nubmer >= 1
    ws_res->factory("s1_CtauRes[0.5, 1e-6, 1.0]");
    ws_res->factory("rS21_CtauRes[1.5, 1.0, 5.0]");
    ws_res->factory("rS32_CtauRes[2.5, 1.0, 5.0]");
    ws_res->factory("rS43_CtauRes[1.5, 1.0, 10.0]");

    ws_res->factory("RooFormulaVar::s2_CtauRes('@0*@1',{rS21_CtauRes,s1_CtauRes})");
    ws_res->factory("RooFormulaVar::s3_CtauRes('@0*@1',{rS32_CtauRes,s2_CtauRes})");
    ws_res->factory("RooFormulaVar::s4_CtauRes('@0*@1',{rS43_CtauRes,s3_CtauRes})");

    // portion of each gauss in total pdf
    ws_res->factory("f_CtauRes[0.2, 0., 1.]");
    ws_res->factory("f2_CtauRes[0.2, 0., 1.]");
    ws_res->factory("f3_CtauRes[0.5, 0., 1.]");

    // build gauss
    TString varName = "ctau3DRes";
    ws_res->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel1_ctauRes", varName.Data(),
                     "ctau1_CtauRes", //"ctau1_CtauRes",
                     "s1_CtauRes"));
    ws_res->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel2_ctauRes", varName.Data(),
                     "ctau2_CtauRes", //"ctau2_CtauRes",
                     "s2_CtauRes"));
    ws_res->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel3_ctauRes", varName.Data(),
                     "ctau3_CtauRes", //"ctau3_CtauRes",
                     "s3_CtauRes"));
    ws_res->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel4_ctauRes", varName.Data(),
                     "ctau4_CtauRes", //"ctau3_CtauRes",
                     "s4_CtauRes"));
    
    // combine gausses
    if (nGauss == 4)
    {
        ws_res->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel43_ctauRes", "GaussModel4_ctauRes", "GaussModel3_ctauRes", "f3_CtauRes"));
        ws_res->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel23_ctauRes", "GaussModel43_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));
    }
    else if (nGauss == 3)
    {
        ws_res->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel23_ctauRes", "GaussModel3_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));
        ws_res->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel23_ctauRes", "f_CtauRes"));
    }
    else if (nGauss == 2)
    {
        ws_res->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel2_ctauRes", "f_CtauRes"));
    }

    // build total fit model
    // strange initializing - I left it due to it's legacy
    RooAddPdf *gaus_res_model = new RooAddPdf("gaus_res_model");
    RooAbsPdf *ctauResModel = ctauResModel = new RooAddPdf("gaus_res_model", "gaus_res_model", *ws_res->pdf("GaussModelCOND_ctauRes"), *ws_res->var("n_sig"));
    ws_res->import(*ctauResModel);


    // ===== set new ctau3DRes range ===== //
    // nominal only use negative values
    TCanvas *c_res = new TCanvas("c_res", "My plots", 600, 600);
    c_res->cd();
    TPad *pad_res1 = new TPad("pad_res1", "pad_res1", 0, 0.25, 0.98, 1.0);
    pad_res1->SetTicks(1, 1);
    pad_res1->Draw();
    pad_res1->cd();

    // new decide new range
    RooPlot *res_frame1 = ws_res->var("ctau3DRes")->frame(Bins(nCtauResBins), Range(ctauResLow, ctauResHigh)); // bins
    res_frame1->SetTitle("");

    TH1D *h_tmp = (TH1D *)ws_res->data("ds_sig_err")->createHistogram(("h_tmp"), *ws_res->var("ctau3DRes"), Binning(nCtauResBins, ctauResLow, ctauResHigh));
    double ctauResMin = -10;
    // double ctauResMax = h_tmp->GetBinCenter(h_tmp->FindLastBinAbove(1,1));
    double ctauResMax = 0;
    cout << "NBins: " << h_tmp->GetNbinsX() << endl;
    for (int i = 0; i < h_tmp->GetNbinsX() / 2; i++)
    {
        // cout<<"Content: "<<h_tmp->GetBinContent(i)<<endl;
        if (h_tmp->GetBinContent(i) <= 0 && h_tmp->GetBinContent(i + 1) <= 0)
        {
            // cout<<"#####"<<i<<": "<<h_tmp->GetBinLowEdge(i+2)<<endl;
            ctauResMin = h_tmp->GetBinLowEdge(i + 2);
        }
        // if(h_tmp->GetBinContent(i)>1)ctauResMax = h_tmp->GetBinCenter(i)+h_tmp->GetBinWidth(i);
    }

    ws_res->var("ctau3DRes")->setRange("ctauResWindow", ctauResMin, 0);
    cout << "Fit Range: " << ctauResMin << " - 0" << endl;

    
    // ===== fit here ===== //
    // buld dataset for fit with new Res range
    RooDataSet *dataToFit = (RooDataSet *)ctauResCutDS->reduce(Form("ctau3DRes>=%.f&&ctau3DRes<=0", ctauResMin))->Clone("dataToFit");
    ws_res->import(*dataToFit);
    
    pad_res1->cd();
    gPad->SetLogy();
    RooPlot *res_frame2 = (RooPlot *)res_frame1->Clone();

    bool isWeighted = ws_res->data("dataToFit")->isWeighted();
    RooFitResult *fit_res = ws_res->pdf("gaus_res_model")->fitTo(*dataToFit, Save(), SumW2Error(isWeighted), Extended(kTRUE), NumCPU(6), PrintLevel(-1));
    fit_res->Print("V");

    // fis parameters
    ws_res->var("s1_CtauRes")->setConstant(kTRUE);
    ws_res->var("rS21_CtauRes")->setConstant(kTRUE);
    ws_res->var("rS32_CtauRes")->setConstant(kTRUE);
    ws_res->var("rS43_CtauRes")->setConstant(kTRUE);
    ws_res->var("f_CtauRes")->setConstant(kTRUE);
    ws_res->var("f2_CtauRes")->setConstant(kTRUE);
    ws_res->var("f3_CtauRes")->setConstant(kTRUE);

    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Plotting - very long but simple *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // setFixedVarsToContantVars(ws_res);
    ws_res->data("dataToFit")->plotOn(res_frame2, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), LineColor(kBlack), MarkerColor(kBlack));
    ws_res->pdf("gaus_res_model")->plotOn(res_frame2, Name("modelHist_ctauRes"), Precision(1e-4), NormRange("ctauResWindow"), Normalization(ws_res->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), LineColor(kBlack));
    ws_res->pdf("gaus_res_model")->plotOn(res_frame2, Name("modelHist_gm1"), Precision(1e-4), NormRange("ctauResWindow"), Normalization(ws_res->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws_res->pdf("GaussModel1_ctauRes")), LineColor(kGreen + 2));
    ws_res->pdf("gaus_res_model")->plotOn(res_frame2, Name("modelHist_gm2"), Precision(1e-4), NormRange("ctauResWindow"), Normalization(ws_res->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws_res->pdf("GaussModel2_ctauRes")), LineColor(kRed + 2));
    ws_res->pdf("gaus_res_model")->plotOn(res_frame2, Name("modelHist_gm3"), Precision(1e-4), NormRange("ctauResWindow"), Normalization(ws_res->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws_res->pdf("GaussModel3_ctauRes")), LineColor(kBlue + 2));
    if (nGauss == 4)
    {
        ws_res->pdf("gaus_res_model")->plotOn(res_frame2, Name("modelHist_gm4"), Precision(1e-4), NormRange("ctauResWindow"), Normalization(ws_res->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws_res->pdf("GaussModel4_ctauRes")), LineColor(kMagenta + 2));
    }
    ws_res->data("ctauResCutDS")->plotOn(res_frame2, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), LineColor(kRed + 2), MarkerColor(kRed + 2));

    // y-axis range
    TH1 *h = ws_res->data("ctauResCutDS")->createHistogram("hist", *ws_res->var("ctau3DRes"), Binning(res_frame1->GetNbinsX(), res_frame1->GetXaxis()->GetXmin(), res_frame1->GetXaxis()->GetXmax()));
    Double_t YMax = h->GetBinContent(h->GetMaximumBin());
    Double_t YMin = 1e99;
    for (int i = 1; i <= h->GetNbinsX(); i++)
        if (h->GetBinContent(i) > 0)
            YMin = min(YMin, h->GetBinContent(i));
    Double_t Yup(0.), Ydown(0.);
    Yup = YMax * TMath::Power((YMax / YMin), (0.5 / (1.0 - 0.5 - 0.2)));
    Ydown = YMin / (TMath::Power((YMax / YMin), (0.2 / (1.0 - 0.5 - 0.2))));
    res_frame2->GetYaxis()->SetRangeUser(Ydown, Yup);
    cout << Ydown << " to " << Yup << endl;
    cout << "###" << endl;

    // check lost event from new Res cut
    Double_t outTot = ws_res->data("ctauResCutDS")->numEntries();
    Double_t outRes = ws_res->data("dataToFit")->numEntries();
    cout << "ctauRes: " << ctauResMin << " ~ " << ctauResMax << endl;
    cout << "Tot evt: (" << outTot << ")" << endl;
    cout << "Res evt: (" << outRes << ")" << endl;
    cout << "lost evt: (" << (outRes * 100) / outTot << ")%, " << outRes << "evts" << "\n\n";

    // draw lines to present Res rantge
    if (outRes > 0.0)
    {
        // TLine   *minline = new TLine(rangeErr[0], 0.0, rangeErr[0], Ydown*TMath::Power((Yup/Ydown),0.4));
        TLine *minline = new TLine(ctauResMin, 0.0, ctauResMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
        minline->SetLineStyle(2);
        minline->SetLineColor(1);
        minline->SetLineWidth(3);
        res_frame2->addObject(minline);
        TLine *maxline = new TLine(ctauResMax, 0.0, ctauResMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
        maxline->SetLineStyle(2);
        maxline->SetLineColor(1);
        maxline->SetLineWidth(3);
        res_frame2->addObject(maxline);
    }
    res_frame2->GetXaxis()->CenterTitle();
    res_frame2->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
    res_frame2->SetFillStyle(4000);
    res_frame2->GetYaxis()->SetTitleOffset(2);
    res_frame2->GetXaxis()->SetLabelSize(0);
    res_frame2->GetXaxis()->SetTitleSize(0);
    res_frame2->Draw();


    // ===== legend ===== //
    TLegend *leg_C = new TLegend(text_x + 0.29, text_y + 0.03, text_x + 0.39, text_y - 0.17);
    leg_C->SetTextSize(text_size);
    leg_C->SetTextFont(43);
    leg_C->SetBorderSize(0);
    leg_C->AddEntry(res_frame2->findObject("dataHist_ctauRes"), "Data", "pe");
    leg_C->AddEntry(res_frame2->findObject("modelHist_ctauRes"), "Total PDF", "l");
    leg_C->AddEntry(res_frame2->findObject("modelHist_gm1"), "Gauss 1", "l");
    leg_C->AddEntry(res_frame2->findObject("modelHist_gm2"), "Gauss 2", "l");
    if (nGauss >= 3)
        leg_C->AddEntry(res_frame2->findObject("modelHist_gm3"), "Gauss 3", "l");
    if (nGauss == 4)
        leg_C->AddEntry(res_frame2->findObject("modelHist_gm4"), "Gauss 4", "l");
    leg_C->Draw("same");
    // cout<<"s2/s1: "<<ws_res->var("s2_CtauRes")->getVal()/ws_res->var("s1_CtauRes")->getVal()<<endl;
    // cout<<"s3/s2: "<<ws_res->var("s3_CtauRes")->getVal()/ws_res->var("s2_CtauRes")->getVal()<<endl;
    // cout << "s1: " << ws_res->var("s1_CtauRes")->getVal() << endl;


    // ===== draw latex ===== //
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y, text_color, text_size);
    if (yLow == 0)
        drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff, text_color, text_size);
    else if (yLow != 0)
        drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff, text_color, text_size);
    drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x, text_y - y_diff * 2, text_color, text_size);
    drawText(Form("Loss: (%.4f%s) %.f evts", (outTot - outRes) * 100 / outTot, "%", outTot - outRes), text_x, text_y - y_diff * 3, text_color, text_size);
    // cout<<"lost evt: ("<<(outRes*100)/outTot<<")%, "<<outRes<<"evts"<<endl;
    drawText(Form("N_{J/#psi} = %.f #pm %.f", ws_res->var("n_sig")->getVal(), ws_res->var("n_sig")->getError()), text_x + 0.5, text_y, text_color, text_size);
    drawText(Form("s1_{Res} = %.4f #pm %.3f", ws_res->var("s1_CtauRes")->getVal(), ws_res->var("s1_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 1, text_color, text_size);
    drawText(Form("(s2/s1)_{Res} = %.3f #pm %.3f", ws_res->var("rS21_CtauRes")->getVal(), ws_res->var("rS21_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 2, text_color, text_size);
    drawText(Form("f_{Res} = %.4f #pm %.4f", ws_res->var("f_CtauRes")->getVal(), ws_res->var("f_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 3, text_color, text_size);
    if (nGauss >= 3)
        drawText(Form("(s3/s2)_{Res} = %.4f #pm %.3f", ws_res->var("rS32_CtauRes")->getVal(), ws_res->var("rS32_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 4, text_color, text_size);
    
    if (nGauss >= 3)
        drawText(Form("f2_{Res} = %.4f #pm %.3f", ws_res->var("f2_CtauRes")->getVal(), ws_res->var("f2_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 5, text_color, text_size);


    // ===== draw pull ===== //
    TPad *pull_pad = new TPad("pull_pad", "pull_pad", 0, 0.001, 0.98, 0.32);
    c_res->cd();
    pull_pad->Draw();
    pull_pad->cd();
    pull_pad->SetTopMargin(0); // Upper and lower plot are joined
    pull_pad->SetBottomMargin(0.67);
    pull_pad->SetBottomMargin(0.4);
    pull_pad->SetFillStyle(4000);
    pull_pad->SetFrameFillStyle(4000);
    pull_pad->SetTicks(1, 1);

    RooPlot *pull_tmp = (RooPlot *)res_frame2->Clone("TMP");
    RooHist *pull_hist = pull_tmp->pullHist("dataHist_ctauRes", "modelHist_ctauRes", true);
    pull_hist->SetMarkerSize(0.8);
    RooPlot *pull_frame = ws_res->var("ctau3DRes")->frame(Title("Pull Distribution"), Bins(nCtauResBins), Range(ctauResLow, ctauResHigh));
    pull_frame->addPlotable(pull_hist, "PX");

    // cosmetics
    pull_frame->SetTitle("");
    pull_frame->SetTitleSize(0);
    pull_frame->GetYaxis()->SetTitleOffset(0.3);
    pull_frame->GetYaxis()->SetTitle("Pull");
    pull_frame->GetYaxis()->SetTitleSize(0.08);
    pull_frame->GetYaxis()->SetLabelSize(0.08);
    pull_frame->GetYaxis()->SetRangeUser(-9, 9);
    pull_frame->GetYaxis()->CenterTitle();

    pull_frame->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}}}");
    pull_frame->GetXaxis()->SetTitleOffset(1.55);
    pull_frame->GetXaxis()->SetLabelOffset(0.04);
    pull_frame->GetXaxis()->SetLabelSize(0.08);
    pull_frame->GetXaxis()->SetTitleSize(0.08);
    pull_frame->GetXaxis()->CenterTitle();

    pull_frame->GetYaxis()->SetTickSize(0.04);
    pull_frame->GetYaxis()->SetNdivisions(404);
    pull_frame->GetXaxis()->SetTickSize(0.03);
    pull_frame->Draw();

    // draw center line
    TLine *lC = new TLine(ctauResLow, 0, ctauResHigh, 0);
    lC->SetLineStyle(1);
    lC->Draw("same");


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


    // ===== draw chi_c, ndf ===== //
    // printChi2(ws_res, pull_pad, frameTMP_C, fit_res, "ctau3DRes", "dataHist_ctauRes", "modelHist_ctauRes", nCtauResBins, false);
    // pull_pad->Update();


    // ===== draw canvas ===== //
    c_res->Update();
    // c_res->Draw();
    c_res->SaveAs(("figs/" + out_ss + ".png").c_str());


    // ===== Export results ===== //
    auto out_file = new TFile(("roots/" + out_ss + ".root").c_str(), "recreate");

    // RooArgSet *fitargs = new RooArgSet();
    // fitargs->add(fit_res->floatParsFinal());
    // RooDataSet *datasetRes = new RooDataSet("datasetRes", "dataset with Resolution Fit result", *fitargs);
    // ws_res->import(*fit_res);

    ctauResModel->Write(); // it's content will be saved as gaus_res_model in root file
    ws_res->Write();
    // gaus_res_model->Write(); // Never do it. It's content was moved to ctauResModel.
    // gaus_res_model->Print("V");
    // //	ctauResCutDS->Write();
    // // datasetRes->Write();
    // // fit_res->Write();
    out_file->Close();

    cout << "\n=================================\n";
    cout << "\n Finish res fit\n";
    cout << "\n=================================\n";
    t->Stop();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
}