#include <iostream>
#include "TH2D.h"
#include <TLorentzVector.h>
#include <TAttMarker.h>
#include "../../headers/commonUtility.h"
#include "../../headers/JpsiUtility.h"
#include "../../headers/HiEvtPlaneList.h"
#include "../../headers/cutsAndBin.h"
#include "../../headers/Style.h"
#include "../../headers/tnp_weight_lowptPbPb_num_den_new.h"
// #include "../headers/tnp_weight_lowptPbPb.h"

void divide_2d_hist(const std::shared_ptr<TH2D> &h_num, const std::shared_ptr<TH2D> &h_den,
                    std::shared_ptr<TH2D> &h_eff);
void draw_2d_hist(const std::shared_ptr<TH2D> &h_eff, string legend, int MCtype);
bool IsAcceptable(double pt, double eta);

void acc_pp(string data_label_ = "all_event", int nevt = -1, int MCtype = 1, bool isTnP = true, bool isPtWeight = true)
{
    TStopwatch *t = new TStopwatch;
    t->Start();
    cout << "\n=================================\n";
    cout << "\n Start pp Acc Calculation\n";
    cout << "\n=================================\n\n";

    using namespace std;
    using namespace hi; // HiEvtPlaneList.h

    // ===== basic setting ===== //
    // gStyle->SetOptStat(0);

    // kinematics
    float ptLow = 0.0, ptHigh = 50.0;
    float yLow = 0.0, yHigh = 2.4;

    // jpsi mass
    float massLow = 2.6;
    float massHigh = 3.5;

    // labeling
    TString date_label = data_label_;

    // make directories for outputs
    gSystem->mkdir(Form("roots/"), kTRUE);
    gSystem->mkdir(Form("figs/"), kTRUE);

    // ===== import oniatree ===== //
    // TChain muon_chain("hionia/myTree");
    TChain muon_chain("hionia/myTree");

    if (MCtype == 1)
    {
        // pp MC prompt
        TString filename = Form("/disk1/Oniatree/Jpsi/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root");
        std::cout << "Adding pp MC prompt sample: " << filename << std::endl;
        muon_chain.Add(filename);
    }
    else if (MCtype == 2)
    {
        // not used
        // pp MC nonprompt
        TString filename = "";
        std::cout << "Adding pp MC nonprompt sample: " << filename << std::endl;

        muon_chain.Add(filename);
    }

    // ===== build histograms ===== //
    vector<double> fwd_pt_bin = {3, 6, 9, 12, 15, 20, 50};
    vector<double> mid_pt_bin = {6.5, 9.0, 12, 15.0, 20, 50};
    vector<double> cos_bin = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    vector<double> phi_bin = {-3.15, -2.835, -2.52, -2.205, -1.89, -1.575, -1.26, -0.945, -0.63, -0.315, 0.0, 0.315, 0.63, 0.945, 1.26, 1.575, 1.89, 2.205, 2.52, 2.835, 3.15};
    vector<double> y_bin = {-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4};

    double n_fwd_pt = fwd_pt_bin.size() - 1;
    double n_mid_pt = mid_pt_bin.size() - 1;
    double n_cos = cos_bin.size() - 1;
    double n_phi = phi_bin.size() - 1;
    double n_rapidity = y_bin.size() - 1;

    // 1d histograms
    auto y_lab_num = std::make_shared<TH1D>("y_lab_num", "", n_rapidity, &y_bin[0]);
    auto y_lab_den = std::make_shared<TH1D>("y_lab_den", "", n_rapidity, &y_bin[0]);

    // 3d histograms
    // ===== fwd ===== //
    auto fwd_lab_num = std::make_shared<TH3D>("fwd_lab_num", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto fwd_lab_den = std::make_shared<TH3D>("fwd_lab_den", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);

    auto fwd_hx_num = std::make_shared<TH3D>("fwd_hx_num", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto fwd_hx_den = std::make_shared<TH3D>("fwd_hx_den", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);

    auto fwd_cs_num = std::make_shared<TH3D>("fwd_cs_num", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto fwd_cs_den = std::make_shared<TH3D>("fwd_cs_den", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);

    // ===== mid ===== //
    auto mid_lab_num = std::make_shared<TH3D>("mid_lab_num", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_mid_pt, &mid_pt_bin[0]);
    auto mid_lab_den = std::make_shared<TH3D>("mid_lab_den", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_mid_pt, &mid_pt_bin[0]);

    auto mid_hx_num = std::make_shared<TH3D>("mid_hx_num", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_mid_pt, &mid_pt_bin[0]);
    auto mid_hx_den = std::make_shared<TH3D>("mid_hx_den", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_mid_pt, &mid_pt_bin[0]);

    auto mid_cs_num = std::make_shared<TH3D>("mid_cs_num", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_mid_pt, &mid_pt_bin[0]);
    auto mid_cs_den = std::make_shared<TH3D>("mid_cs_den", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_mid_pt, &mid_pt_bin[0]);

    // =====  make hists to save error ===== //
    y_lab_num->Sumw2();
    y_lab_den->Sumw2();

    // fwd
    fwd_lab_num->Sumw2();
    fwd_lab_den->Sumw2();
    fwd_hx_num->Sumw2();
    fwd_hx_den->Sumw2();
    fwd_cs_num->Sumw2();
    fwd_cs_den->Sumw2();

    // mid
    mid_lab_num->Sumw2();
    mid_lab_den->Sumw2();
    mid_hx_num->Sumw2();
    mid_hx_den->Sumw2();
    mid_cs_num->Sumw2();
    mid_cs_den->Sumw2();


    // ===== build pT weighting function ===== //
    // import pT reweighting function
    auto fPtW_mid = make_shared<TFile>("WeightedFcN_fit/ratioDataMC_AA_PromptJpsi_DATA_All_y.root", "read");
    auto fPtW_fwd = make_shared<TFile>("WeightedFcN_fit/ratioDataMC_AA_PromptJpsi_DATA_Forward_y.root", "read");

    TF1 *fptw_mid = (TF1 *)fPtW_mid->Get("dataMC_Ratio1");
    TF1 *fptw_fwd = (TF1 *)fPtW_fwd->Get("dataMC_Ratio1");

    // build local function - only for plotting????
    // 이거 왜 있지?
    TF1 *f1 = new TF1("f1", "[0]*TMath::Erf((x-[1])/[2])", 3, 50);
    f1->SetParameters(214, 22, 14);
    f1->SetParLimits(0, 0, 500);
    f1->SetParLimits(1, -1, 500);
    f1->SetParLimits(2, 0, 500);

    // ===== connect to input tree ===== //
    // muon_chain
    const int maxBranchSize = 1000;
    Int_t Gen_QQ_size;
    TClonesArray *Gen_QQ_4mom;
    TClonesArray *Gen_QQ_mupl_4mom;
    TClonesArray *Gen_QQ_mumi_4mom;
    Float_t Gen_QQ_ctau3D[1000]; //[Gen_QQ_size]
    Float_t Gen_QQ_ctau[1000];   //[Gen_QQ_size]

    TBranch *b_Gen_QQ_size;
    TBranch *b_Gen_QQ_4mom;      //!
    TBranch *b_Gen_QQ_mupl_4mom; //!
    TBranch *b_Gen_QQ_mumi_4mom; //!
    TBranch *b_Gen_QQ_ctau3D;    //!
    TBranch *b_Gen_QQ_ctau;      //!

    Gen_QQ_4mom = 0;
    Gen_QQ_mupl_4mom = 0;
    Gen_QQ_mumi_4mom = 0;

    muon_chain.SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
    muon_chain.SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
    muon_chain.SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
    muon_chain.SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
    muon_chain.SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
    muon_chain.SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);

    // ===== declare Lorentz vectors ===== //
    TLorentzVector *JP_Gen = nullptr;
    TLorentzVector *mupl_Gen = nullptr;
    TLorentzVector *mumi_Gen = nullptr;

    // ===== set beam information ===== //
    double sqrt_S_NN = 5.02; // TeV
    // Assuming p1 = -p2, two hadron beams are symmetric and and E >> m0
    // Mendelstam value reads sqrt_S_NN = 2 * E.
    double beam1_E = 5.02 / 2; // One beam have half of the energy
    beam1_E *= 1000;           // change TeV->GeV
    double p1_mag = beam1_E;   // Approximation for E >> m0.

    double beam2_E = beam1_E; // symmetry
    double p2_mag = -p1_mag;

    // build vectors
    TVector3 p1_3d_lab(0, 0, p1_mag);
    TVector3 p2_3d_lab(0, 0, p2_mag);

    auto *p1_lab = new TLorentzVector(p1_3d_lab, beam1_E);
    auto *p2_lab = new TLorentzVector(p2_3d_lab, beam2_E);

    // ===== prepare weighting variables ===== //
    double weight = 1;
    double pt_weight = 1;

    // ===== dimuon counters ===== //
    int count_den = 0;
    int count_num = 0;

    if (nevt == -1)
        nevt = muon_chain.GetEntries();

    cout << "\nTotal Events : " << nevt << endl;

    // ===== main loop start ===== //
    for (int iev = 0; iev < nevt; ++iev)
    {
        if (iev % 100000 == 0)
            cout << ">>>>> EVENT " << iev << " / " << muon_chain.GetEntries() << " (" << (int)(100. * iev / muon_chain.GetEntries()) << "%)" << endl;

        muon_chain.GetEntry(iev);

        // ===== gen loop start - denominator ===== //
        for (int igen = 0; igen < Gen_QQ_size; igen++)
        {
            JP_Gen = (TLorentzVector *)Gen_QQ_4mom->At(igen);
            mupl_Gen = (TLorentzVector *)Gen_QQ_mupl_4mom->At(igen);
            mumi_Gen = (TLorentzVector *)Gen_QQ_mumi_4mom->At(igen);

            Double_t Rapidity_gen = fabs(JP_Gen->Rapidity());

            // some kinematic varialbes in lab
            double cos_lab_ = mupl_Gen->CosTheta();
            double phi_lab_ = mupl_Gen->Phi();

            pt_weight = 1;
            if (isPtWeight)
            {
                if (Rapidity_gen > 1.6 && Rapidity_gen < 2.4)
                    pt_weight = fptw_fwd->Eval(JP_Gen->Pt());
                else
                    pt_weight = fptw_mid->Eval(JP_Gen->Pt());
            }

            // ===== Acc denominator cut ===== //
            if (!(fabs(JP_Gen->Rapidity()) < 2.4))
                continue;

            // ===== calculate polarization angles - Gen ===== //
            // polarization unit vectors in lab

            // Boost to Quarkonia rest frame
            std::unique_ptr<TLorentzVector> mupl_qrf_gen = std::make_unique<TLorentzVector>(*mupl_Gen);
            mupl_qrf_gen->Boost(-JP_Gen->BoostVector());
            std::unique_ptr<TLorentzVector> mumi_qrf_gen = std::make_unique<TLorentzVector>(*mumi_Gen);
            mumi_qrf_gen->Boost(-JP_Gen->BoostVector());
            std::unique_ptr<TLorentzVector> p1_qrf_gen = std::make_unique<TLorentzVector>(*p1_lab);
            p1_qrf_gen->Boost(-JP_Gen->BoostVector());
            std::unique_ptr<TLorentzVector> p2_qrf_gen = std::make_unique<TLorentzVector>(*p2_lab);
            p2_qrf_gen->Boost(-JP_Gen->BoostVector());

            // u_p1, u_p2
            TVector3 u_p1 = p1_qrf_gen->Vect();
            u_p1 = u_p1.Unit();
            TVector3 u_p2 = p2_qrf_gen->Vect();
            u_p2 = u_p2.Unit();

            // unit vector of mupl
            TVector3 u_mupl = mupl_qrf_gen->Vect(); // get px, py, pz
            u_mupl = u_mupl.Unit();

            // y axis for frames except EP frame
            TVector3 uy_pol = u_p1.Cross(u_p2); // p1 cross p2
            uy_pol = uy_pol.Unit();

            // CS
            TVector3 uz_cs = u_p1 - u_p2;
            uz_cs = uz_cs.Unit();
            TVector3 ux_cs = uy_pol.Cross(uz_cs);
            ux_cs = ux_cs.Unit();
            float cos_cs_ = u_mupl.Dot(uz_cs);
            float phi_cs_ = TMath::ATan2(u_mupl.Dot(uy_pol), u_mupl.Dot(ux_cs));

            // HX
            TVector3 uz_hx = JP_Gen->Vect();
            uz_hx = uz_hx.Unit();
            TVector3 ux_hx = uy_pol.Cross(uz_hx); // y cross z
            ux_hx = ux_hx.Unit();
            float cos_hx_ = u_mupl.Dot(uz_hx);
            float phi_hx_ = TMath::ATan2(u_mupl.Dot(uy_pol), u_mupl.Dot(ux_hx));

            // PX
            TVector3 uz_px = u_p1 + u_p2;
            uz_px = uz_px.Unit();
            TVector3 ux_px = uy_pol.Cross(uz_px);
            ux_px = ux_px.Unit();
            float cos_px_ = u_mupl.Dot(uz_px);
            float phi_px_ = TMath::ATan2(u_mupl.Dot(uy_pol), u_mupl.Dot(ux_px));

            // ===== fill the denominator ===== //
            ++count_den;
            // rapidity
            y_lab_den->Fill(JP_Gen->Rapidity(), weight * pt_weight);

            if (Rapidity_gen > 1.6 && Rapidity_gen < 2.4)
            {
                fwd_lab_den->Fill(cos_lab_, phi_lab_, JP_Gen->Pt(), weight * pt_weight);
                fwd_hx_den->Fill(cos_hx_, phi_hx_, JP_Gen->Pt(), weight * pt_weight);
                fwd_cs_den->Fill(cos_cs_, phi_cs_, JP_Gen->Pt(), weight * pt_weight);
            }
            else
            {
                mid_lab_den->Fill(cos_lab_, phi_lab_, JP_Gen->Pt(), weight * pt_weight);
                mid_hx_den->Fill(cos_hx_, phi_hx_, JP_Gen->Pt(), weight * pt_weight);
                mid_cs_den->Fill(cos_cs_, phi_cs_, JP_Gen->Pt(), weight * pt_weight);
            }

            // ===== numerator cuts ===== //
            bool mu1pass = IsAcceptable(mupl_Gen->Pt(), fabs(mupl_Gen->Eta()));
            bool mu2pass = IsAcceptable(mumi_Gen->Pt(), fabs(mumi_Gen->Eta()));

            if (mu1pass != true || mu2pass != true)
                continue;

            // if (!(JP_Gen->Pt() > 0 && JP_Gen->Pt() < 50 && fabs(JP_Gen->Rapidity()) < 2.4 && IsAcceptable(mupl_Gen->Pt(), fabs(mupl_Gen->Eta())) && IsAcceptable(mumi_Gen->Pt(), fabs(mumi_Gen->Eta()))))
            //     continue;

            // if (!(fabs(JP_Gen->Rapidity()) < 2.4 && fabs(mupl_Gen->Eta()) < 2.4 && fabs(mumi_Gen->Eta()) < 2.4))
            // continue;

            // if (Gen_mu_charge[Gen_QQ_mupl_idx[igen]] * Gen_mu_charge[Gen_QQ_mumi_idx[igen]] > 0)
            //     continue;

            // ===== fill the numerator ===== //
            ++count_num;
            y_lab_num->Fill(JP_Gen->Rapidity(), weight * pt_weight);

            if (Rapidity_gen > 1.6 && Rapidity_gen < 2.4)
            {
                fwd_lab_num->Fill(cos_lab_, phi_lab_, JP_Gen->Pt(), weight * pt_weight);
                fwd_hx_num->Fill(cos_hx_, phi_hx_, JP_Gen->Pt(), weight * pt_weight);
                fwd_cs_num->Fill(cos_cs_, phi_cs_, JP_Gen->Pt(), weight * pt_weight);
            }
            else
            {
                mid_lab_num->Fill(cos_lab_, phi_lab_, JP_Gen->Pt(), weight * pt_weight);
                mid_hx_num->Fill(cos_hx_, phi_hx_, JP_Gen->Pt(), weight * pt_weight);
                mid_cs_num->Fill(cos_cs_, phi_cs_, JP_Gen->Pt(), weight * pt_weight);
            }
        } // end of gen loop
    } // end of main loop

    cout << "count_den: " << count_den << endl;
    cout << "count_num: " << count_num << endl;


    // ===== save output ===== //
    TString outFileName = Form("roots/acc_pp_Jpsi_PtW%d_tnp%d_%s.root", isPtWeight, isTnP, date_label.Data());
    if (MCtype == 2)
        outFileName = Form("roots/acc_pp_Jpsi_PtW%d_tnp%d_%s.root", isPtWeight, isTnP, date_label.Data());
    TFile *outFile = new TFile(outFileName, "RECREATE");

    // histograms
    y_lab_num->Write();
    y_lab_den->Write();

    // fwd
    fwd_lab_num->Write();
    fwd_lab_den->Write();
    fwd_hx_num->Write();
    fwd_hx_den->Write();
    fwd_cs_num->Write();
    fwd_cs_den->Write();

    // mid
    mid_lab_num->Write();
    mid_lab_den->Write();
    mid_hx_num->Write();
    mid_hx_den->Write();
    mid_cs_num->Write();
    mid_cs_den->Write();

    outFile->Close();

    cout << "\n=================================\n";
    cout << "\n Finish pp Acc  Calculation\n";
    cout << "\n=================================\n\n";
    t->Stop();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
}

void divide_2d_hist(const std::shared_ptr<TH2D> &h_num, const std::shared_ptr<TH2D> &h_den,
                    std::shared_ptr<TH2D> &h_eff)
{
    if (h_num->GetNbinsX() != h_den->GetNbinsX() || h_num->GetNbinsY() != h_den->GetNbinsY())
    {
        std::cerr << "Error: The histograms have different binning!" << std::endl;
        return;
    }

    for (int i = 1; i <= h_num->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= h_num->GetNbinsY(); ++j)
        {
            double content1 = h_num->GetBinContent(i, j);
            double content2 = h_den->GetBinContent(i, j);

            if (content2 != 0)
            {
                double eff = content1 / content2;
                h_eff->SetBinContent(i, j, eff);
            }
            else
                h_eff->SetBinContent(i, j, 0);
        }
    }
}

void draw_2d_hist(const std::shared_ptr<TH2D> &h_eff, string legend, int MCtype)
{
    TCanvas cpt_eff("cpt_eff", "cpt_eff", 0, 0, 900, 800);
    cpt_eff.cd();
    h_eff->SetTitle("Acc x Eff");
    h_eff->SetMarkerStyle(24);
    h_eff->SetMarkerColor(1);
    h_eff->SetLineColor(1);
    h_eff->GetXaxis()->SetTitle("cos#theta_{EP}");
    h_eff->GetYaxis()->SetTitle("p_{T}[GeV/c]");
    // h_eff->GetYaxis()->SetRangeUser(0., 1.3);
    h_eff->SetMinimum(0);
    h_eff->SetMaximum(1);
    h_eff->Draw("colz");
}

bool IsAcceptable(double pt, double eta)
{
    if ((fabs(eta) < 1.2 && pt > 3.5) || ((fabs(eta) >= 1.2 && fabs(eta) < 2.1) && pt > (5.47 - 1.89 * fabs(eta))) || ((fabs(eta) >= 2.1 && fabs(eta) < 2.4) && pt > 1.5))
        return true;
    //          if ((fabs(eta) < 1.0 && pt>3.4) || ((fabs(eta) >=1.0 && fabs(eta) <1.6) && pt>(5.8-2.4*fabs(eta))) || ((fabs(eta) >=1.6 && fabs(eta) <2.4) && pt>3.3667-7.0/9.0*fabs(eta))) return true;//2010 version for 2.76 TeV

    else
        return false;
}