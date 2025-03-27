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

void test2_2D_correction_numbers(string data_label_ = "test2", int nevt = 100000, int MCtype = 1, int cLow = 0, int cHigh = 180, bool isTnP = true, bool isPtWeight = true)
{
    TStopwatch *t = new TStopwatch;
    t->Start();
    cout << "\n=================================\n";
    cout << "\n Start Acc(x)Eff Calculation\n";
    cout << "\n=================================\n\n";

    using namespace std;
    using namespace hi; // HiEvtPlaneList.h

    // Example of using event plane namespace
    cout << " Index of " << EPNames[HFm2] << " = " << HFm2 << endl;
    cout << " Index of " << EPNames[HFp2] << " = " << HFp2 << endl;
    cout << " Index of " << EPNames[trackmid2] << " = " << trackmid2 << endl;

    // ===== basic setting ===== //
    // gStyle->SetOptStat(0);

    // kinematics
    float ptLow = 0.0, ptHigh = 50.0;
    float yLow = 0.0, yHigh = 2.4;

    int kTrigSel = kTrigJpsi; // jpsi=12,Upsilon=13
    float muPtCut = 0;        // 3.5, 1.8
    int kL2filter = 16;       // jpsi=16,Upsilon=38
    int kL3filter = 17;       // jpsi=17,Upsilon=39

    // jpsi mass
    float massLow = 2.6;
    float massHigh = 3.5;

    // labeling
    TString date_label = data_label_;

    // ===== import oniatree ===== //
    TChain muon_chain("myTree");
    TChain ep_chain("tree");

    if (MCtype == 1)
        for (int i = 1; i <= 2; ++i)
        {
            // PbPb MC prompt
            TString filename = Form("/disk1/Oniatree/Jpsi/OniatreeMC_JPsi_pThat2_TunedCP5_HydjetDrumMB_5p02TeV_Pythia8_part%d.root", i);
            std::cout << "Adding PbPb MC prompt sample: " << filename << std::endl;

            muon_chain.Add(filename);
            ep_chain.Add(filename);
        }
    else if (MCtype == 2)
    {
        // PbPb MC nonprompt
        TString filename = "/disk1/Oniatree/Jpsi/OniatreeMC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8.root";
        std::cout << "Adding PbPb MC nonprompt sample: " << filename << std::endl;

        muon_chain.Add(filename);
        ep_chain.Add(filename);
    }

    // ===== build histograms ===== //
    // pT 빈은 forward, mid 다르게
    vector<double> fwd_pt_bin = {0, 3, 6.5, 9, 12, 50};
    vector<double> mid_pt_bin = {0, 6.5, 9, 12, 15, 20, 25, 50};
    vector<double> cos_bin = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

    double n_fwd_pt = fwd_pt_bin.size() - 1;
    double n_mid_pt = mid_pt_bin.size() - 1;
    double n_cos = cos_bin.size() - 1;

    // 1d histograms
    auto c60_180_fwd_pt_num = std::make_shared<TH1D>("c60_180_fwd_pt_num", "", n_fwd_pt, &fwd_pt_bin[0]);
    auto c60_180_fwd_pt_den = std::make_shared<TH1D>("c60_180_fwd_pt_den", "", n_fwd_pt, &fwd_pt_bin[0]);
    auto eff_c60_180_fwd1_pt = std::make_shared<TH1D>("eff_c60_180_fwd1_pt", "", n_fwd_pt, &fwd_pt_bin[0]);

    auto c60_180_fwd_cos_ep_num1 = std::make_shared<TH1D>("c60_180_fwd_cos_ep_num1", "", n_cos, &cos_bin[0]);
    auto c60_180_fwd_cos_ep_den1 = std::make_shared<TH1D>("c60_180_fwd_cos_ep_den1", "", n_cos, &cos_bin[0]);
    auto c60_180_fwd_cos_ep_num2 = std::make_shared<TH1D>("c60_180_fwd_cos_ep_num2", "no pT weight", n_cos, &cos_bin[0]);
    auto c60_180_fwd_cos_ep_den2 = std::make_shared<TH1D>("c60_180_fwd_cos_ep_den2", "no pT weight", n_cos, &cos_bin[0]);
    auto eff_c60_180_fwd1_cent = std::make_shared<TH1D>("eff_c60_180_fwd1_cent", "", n_cos, &cos_bin[0]);

    // 2d histograms
    auto c60_180_fwd_num1 = std::make_shared<TH2D>("c60_180_fwd_num1", "", n_cos, &cos_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto c60_180_fwd_den1 = std::make_shared<TH2D>("c60_180_fwd_den1", "", n_cos, &cos_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto eff_c60_180_fwd1 = std::make_shared<TH2D>("eff_c60_180_fwd1", "", n_cos, &cos_bin[0], n_fwd_pt, &fwd_pt_bin[0]);

    auto c60_180_fwd_num2 = std::make_shared<TH2D>("c60_180_fwd_num2", "no pT weight", n_cos, &cos_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto c60_180_fwd_den2 = std::make_shared<TH2D>("c60_180_fwd_den2", "no pT weight", n_cos, &cos_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto eff_c60_180_fwd2 = std::make_shared<TH2D>("eff_c60_180_fwd2", "no pT weight", n_cos, &cos_bin[0], n_fwd_pt, &fwd_pt_bin[0]);

    // =====  make hists to save error ===== //
    // 1d hist
    c60_180_fwd_pt_num->Sumw2();
    c60_180_fwd_pt_den->Sumw2();

    c60_180_fwd_cos_ep_num1->Sumw2();
    c60_180_fwd_cos_ep_den1->Sumw2();

    c60_180_fwd_cos_ep_num2->Sumw2();
    c60_180_fwd_cos_ep_den2->Sumw2();

    // 2d hist
    c60_180_fwd_num1->Sumw2();
    c60_180_fwd_den1->Sumw2();
    eff_c60_180_fwd1->Sumw2();

    c60_180_fwd_num2->Sumw2();
    c60_180_fwd_den2->Sumw2();
    eff_c60_180_fwd2->Sumw2();

    // ===== build pT weighting function ===== //
    // import pT reweighting function
    TFile *fPtW = new TFile("WeightedFcN_fit/ratioDataMC_AA_PromptJpsi_DATA_All_y.root", "read");
    if (MCtype == 2)
        TFile *fPtW = new TFile("WeightedFcN_fit/ratioDataMC_AA_btojpsi_DATA_1s.root", "read");
    TF1 *fptw = (TF1 *)fPtW->Get("dataMC_Ratio1");

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
    Int_t Centrality;
    ULong64_t HLTriggers;
    Int_t Gen_QQ_size;
    Int_t Gen_mu_size;
    TClonesArray *Gen_QQ_4mom;
    TClonesArray *Gen_mu_4mom;
    ULong64_t Gen_QQ_trig[maxBranchSize];  //[Gen_QQ_size]
    Float_t Gen_QQ_VtxProb[maxBranchSize]; //[Gen_QQ_size]
    TBranch *b_Centrality;                 //!
    TBranch *b_HLTriggers;                 //!
    TBranch *b_Gen_QQ_size;                //!
    TBranch *b_Gen_mu_size;                //!
    TBranch *b_Gen_QQ_4mom;                //!
    TBranch *b_Gen_mu_4mom;                //!
    TBranch *b_Gen_QQ_trig;                //!
    TBranch *b_Gen_QQ_VtxProb;             //!

    Gen_QQ_4mom = 0;
    Gen_mu_4mom = 0;
    muon_chain.SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
    muon_chain.SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
    muon_chain.SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
    muon_chain.SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

    //  muon id
    Int_t Gen_QQ_mupl_idx[maxBranchSize];
    Int_t Gen_QQ_mumi_idx[maxBranchSize];
    TBranch *b_Gen_QQ_mupl_idx;
    TBranch *b_Gen_QQ_mumi_idx;
    muon_chain.SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
    muon_chain.SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);

    Int_t Gen_mu_charge[maxBranchSize];
    TBranch *b_Gen_mu_charge; //!
    muon_chain.SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);

    Int_t Reco_QQ_size;
    Int_t Reco_mu_size;
    TClonesArray *Reco_QQ_4mom;
    TClonesArray *Reco_mu_4mom;
    ULong64_t Reco_QQ_trig[maxBranchSize];  //[Reco_QQ_size]
    ULong64_t Reco_mu_trig[maxBranchSize];  //[Reco_QQ_size]
    Float_t Reco_QQ_VtxProb[maxBranchSize]; //[Reco_QQ_size]
    TBranch *b_Reco_QQ_size;                //!
    TBranch *b_Reco_mu_size;                //!
    TBranch *b_Reco_QQ_4mom;                //!
    TBranch *b_Reco_mu_4mom;                //!
    TBranch *b_Reco_QQ_trig;                //!
    TBranch *b_Reco_mu_trig;                //!
    TBranch *b_Reco_QQ_VtxProb;             //!

    Reco_QQ_4mom = 0;
    Reco_mu_4mom = 0;
    muon_chain.SetBranchAddress("Centrality", &Centrality, &b_Centrality);
    muon_chain.SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
    muon_chain.SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
    muon_chain.SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
    muon_chain.SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
    muon_chain.SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
    muon_chain.SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
    muon_chain.SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
    muon_chain.SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

    //  muon id
    Int_t Reco_QQ_mupl_idx[maxBranchSize];
    Int_t Reco_QQ_mumi_idx[maxBranchSize];
    TBranch *b_Reco_QQ_mupl_idx;
    TBranch *b_Reco_QQ_mumi_idx;
    muon_chain.SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
    muon_chain.SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);

    Float_t Reco_mu_dxy[maxBranchSize]; //[Reco_mu_size]
    TBranch *b_Reco_mu_dxy;             //!
    muon_chain.SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
    Float_t Reco_mu_dz[maxBranchSize]; //[Reco_mu_size]
    TBranch *b_Reco_mu_dz;             //!
    muon_chain.SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
    Int_t Reco_mu_nTrkWMea[maxBranchSize]; //[Reco_mu_size]
    TBranch *b_Reco_mu_nTrkWMea;           //!
    muon_chain.SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
    Int_t Reco_mu_nPixWMea[maxBranchSize]; //[Reco_mu_size]
    TBranch *b_Reco_mu_nPixWMea;           //!
    muon_chain.SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
    Int_t Reco_QQ_sign[maxBranchSize]; //[Reco_QQ_size]
    TBranch *b_Reco_QQ_sign;           //!
    muon_chain.SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);

    Int_t Reco_mu_SelectionType[maxBranchSize];
    TBranch *b_Reco_mu_SelectionType;
    muon_chain.SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

    Int_t Reco_mu_whichGen[maxBranchSize];
    TBranch *b_Reco_mu_whichGen;
    Float_t Gen_weight;
    TBranch *b_Gen_weight;
    muon_chain.SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
    muon_chain.SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);

    // ===== set event plane branches =====
    const int nEP = 29; // number of event planes in the tree
    double epang[nEP];
    TBranch *b_epang;
    ep_chain.SetBranchAddress("epang", epang, &b_epang);

    // ===== declare Lorentz vectors ===== //
    TLorentzVector *JP_Gen = nullptr;
    TLorentzVector *mupl_Gen = nullptr;
    TLorentzVector *mumi_Gen = nullptr;

    TLorentzVector *JP_Reco = nullptr;
    TLorentzVector *mupl_Reco = nullptr;
    TLorentzVector *mumi_Reco = nullptr;

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
    double tnp_weight = 1;
    double tnp_trig_weight = 1;
    double tnp_trig_weight_mupl = -1;
    double tnp_trig_weight_mumi = -1;
    double tnp_trig_weight_muplL2_num = -1;
    double tnp_trig_weight_muplL3_num = -1;
    double tnp_trig_weight_mumiL2_num = -1;
    double tnp_trig_weight_mumiL3_num = -1;
    double tnp_trig_weight_muplL2_den = -1;
    double tnp_trig_weight_muplL3_den = -1;
    double tnp_trig_weight_mumiL2_den = -1;
    double tnp_trig_weight_mumiL3_den = -1;
    double tnp_trig_weight_num = 1;
    double tnp_trig_weight_den = 1;
    double pt_weight = 1;

    double tnp_trig_dimu = -1;
    // legacy
    TH2D *hpt_tnp_trig_fwd = new TH2D("hpt_tnp_trig_fwd", "hpt_tnp_trig_fwd", n_fwd_pt, &fwd_pt_bin[0], 40, 0, 2);

    // ===== dimuon counters ===== //
    int count = 0;
    int counttnp = 0;
    int countPtW = 0;

    if (nevt == -1)
        nevt = muon_chain.GetEntries();

    cout << "\nTotal Events : " << nevt << endl;

    // ===== main loop start ===== //
    for (int iev = 0; iev < nevt; ++iev)
    {
        if (iev % 100000 == 0)
            cout << ">>>>> EVENT " << iev << " / " << muon_chain.GetEntries() << " (" << (int)(100. * iev / muon_chain.GetEntries()) << "%)" << endl;

        muon_chain.GetEntry(iev);
        ep_chain.GetEntry(iev);

        if (Centrality > cHigh || Centrality < cLow)
            continue;
        weight = findNcoll(Centrality) * Gen_weight;

        // get rpAngles
        float epHFm2 = epang[HFm2];
        float epHFp2 = epang[HFp2];

        // ===== gen loop start - denominator ===== //
        for (int igen = 0; igen < Gen_QQ_size; igen++)
        {
            JP_Gen = (TLorentzVector *)Gen_QQ_4mom->At(igen);
            mupl_Gen = (TLorentzVector *)Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
            mumi_Gen = (TLorentzVector *)Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

            Double_t Rapidity_gen = fabs(JP_Gen->Rapidity());

            pt_weight = 1;
            if (isPtWeight)
                pt_weight = fptw->Eval(JP_Gen->Pt());

            // very simple cut for Polarization
            // other study usually inclues muon kinematic cuts
            if (!(fabs(JP_Gen->Rapidity()) < 2.4))
                continue;

            if (Gen_mu_charge[Gen_QQ_mupl_idx[igen]] * Gen_mu_charge[Gen_QQ_mumi_idx[igen]] > 0)
                continue;

            // ===== calculate polarization angles - Gen ===== //
            // polarization unit vectors in lab
            TVector3 uy_lab(0, 1, 0);
            TVector3 uz_lab(0, 0, 1);
            TVector3 uz_ep = uy_lab;

            // consider auto correleation
            float eta_ = JP_Gen->Eta();
            float epAng = 0;
            if (eta_ > 0)
                epAng = epHFm2;
            else
                epAng = epHFp2;

            // get z-axis in EP frame
            uz_ep.Rotate(epAng, uz_lab); // uz_ep = y rotated as epAng in lab frame

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

            // EP
            TVector3 uy_ep = p1_3d_lab;
            uy_ep = uy_ep.Unit();
            TVector3 ux_ep = uy_ep.Cross(uz_ep);
            ux_ep = ux_ep.Unit();
            float cos_ep_ = u_mupl.Dot(uz_ep);
            float phi_ep_ = TMath::ATan2(u_mupl.Dot(uy_pol), u_mupl.Dot(ux_ep));

            // ===== fill the denominator ===== //
            if (Centrality > 0 && Centrality < 180 && JP_Gen->Pt() > 3)
            {
                c60_180_fwd_pt_den->Fill(JP_Gen->Pt(), weight * pt_weight);
                c60_180_fwd_cos_ep_den1->Fill(cos_ep_, weight * pt_weight);
                c60_180_fwd_cos_ep_den2->Fill(cos_ep_, weight);
                c60_180_fwd_den1->Fill(cos_ep_, JP_Gen->Pt(), weight * pt_weight);
                c60_180_fwd_den2->Fill(cos_ep_, JP_Gen->Pt(), weight);
            }
        } // end of gen loop

        // ===== Reco loop strat - numerator ===== //
        bool HLTPass = false;
        if ((HLTriggers & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))
            HLTPass = true;

        for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
        {
            JP_Reco = (TLorentzVector *)Reco_QQ_4mom->At(irqq);
            mupl_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
            mumi_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

            bool HLTFilterPass = false;
            if ((Reco_QQ_trig[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))
                HLTFilterPass = true;

            if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1)
                continue;
            if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1)
                continue;

            bool passMuonTypePl = true;
            passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 1)));
            passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 3)));

            bool passMuonTypeMi = true;
            passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 1)));
            passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 3)));

            bool muplSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
                (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
                (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
                (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]]) < 0.3) &&
                (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]]) < 20.) &&
                passMuonTypePl //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
            );

            bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
                (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
                (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
                (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]]) < 0.3) &&
                (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]]) < 20.) &&
                passMuonTypeMi //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
            );

            if (!(muplSoft && mumiSoft))
                continue;
            if (Reco_QQ_VtxProb[irqq] < 0.01)
                continue;
            if (Reco_QQ_sign[irqq] != 0)
                continue;

            Double_t Rapidity_reco = fabs(JP_Reco->Rapidity());

            if (Rapidity_reco < 1.2)
                muPtCut = 3.5;
            else if (Rapidity_reco > 1.2 && Rapidity_reco < 2.1)
                muPtCut = (5.47 - 1.89 * Rapidity_reco);
            else if (Rapidity_reco > 2.1 && Rapidity_reco < 2.4)
                muPtCut = 1.5;

            if (!((mupl_Reco->Pt() > muPtCut && fabs(mupl_Reco->Eta()) < 2.4) && (mumi_Reco->Pt() > muPtCut && fabs(mumi_Reco->Eta()) < 2.4) && (fabs(JP_Reco->Rapidity()) < 2.4 && JP_Reco->Pt() < 50 && JP_Reco->Pt() > 3 && JP_Reco->M() > massLow && JP_Reco->M() < massHigh)))
                continue;

            if (HLTPass == true && HLTFilterPass == true)
                count++;
            if (isTnP)
            {
                tnp_weight = 1;
                tnp_trig_weight = 1;
                tnp_trig_weight_mupl = -1;
                tnp_trig_weight_mumi = -1;
                tnp_trig_weight_muplL2_num = -1;
                tnp_trig_weight_muplL3_num = -1;
                tnp_trig_weight_mumiL2_num = -1;
                tnp_trig_weight_mumiL3_num = -1;
                tnp_trig_weight_muplL2_den = -1;
                tnp_trig_weight_muplL3_den = -1;
                tnp_trig_weight_mumiL2_den = -1;
                tnp_trig_weight_mumiL3_den = -1;
                tnp_trig_weight_num = 1;
                tnp_trig_weight_den = 1;

                tnp_weight = tnp_weight * tnp_weight_muid_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0) * tnp_weight_muid_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0); // mu id
                tnp_weight = tnp_weight * tnp_weight_trk_pbpb(mupl_Reco->Eta(), 0) * tnp_weight_trk_pbpb(mumi_Reco->Eta(), 0);                                     // inner tracker

                // Trigger part
                if (!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))))
                {
                    //         cout << "irqq : " << irqq << " - iev : " << iev << endl;
                    //         cout << "TnP ERROR !!!! ::: No matched L2 filter1 " << endl;
                    continue;
                }
                bool mupl_L2Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
                bool mupl_L3Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
                bool mumi_L2Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
                bool mumi_L3Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
                if (mupl_L2Filter == false || mumi_L2Filter == false)
                {
                    cout << "TnP ERROR !!!! ::: No matched L2 filter2 " << endl;
                    cout << endl;
                }

                bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
                bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
                bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
                bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
                bool SelDone = false;

                /*if( mupl_isL2 && mumi_isL3){
              tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
              tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);
              SelDone = true;
              }
              else if( mupl_isL3 && mumi_isL2){
              tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
              tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
              SelDone = true;
              }
              else if( mupl_isL3 && mumi_isL3){
              int t[2] = {-1,1}; // mupl, mumi
              int l = rand() % (2);
                //pick up what will be L2
                if(t[l]==-1){
                tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
                tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);
                }
                else if(t[l]==1){
                tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
                tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
                }
                else {cout << "ERROR :: No random selection done !!!!" << endl; continue;}
                SelDone = true;
                }*/
                if (mupl_isL2 && mumi_isL3)
                {
                    tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
                    tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);
                    SelDone = true;
                    tnp_trig_weight = tnp_trig_weight_mupl * tnp_trig_weight_mumi;
                }
                else if (mupl_isL3 && mumi_isL2)
                {
                    tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
                    tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
                    SelDone = true;
                    tnp_trig_weight = tnp_trig_weight_mupl * tnp_trig_weight_mumi;
                }
                else if (mupl_isL3 && mumi_isL3)
                {
                    // tnp_trig_weight_muplL2_num = -1;
                    // tnp_trig_weight_muplL3_num = -1;
                    // tnp_trig_weight_mumiL2_num = -1;
                    // tnp_trig_weight_mumiL3_num = -1;
                    // tnp_trig_weight_muplL2_den = -1;
                    // tnp_trig_weight_muplL3_den = -1;
                    // tnp_trig_weight_mumiL2_den = -1;
                    // tnp_trig_weight_mumiL3_den = -1;
                    // tnp_trig_weight_num = 1;
                    // tnp_trig_weight_den = 1;

                    tnp_trig_weight_muplL2_num = tnp_weight_trg_pbpb_num(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
                    tnp_trig_weight_muplL3_num = tnp_weight_trg_pbpb_num(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
                    tnp_trig_weight_mumiL2_num = tnp_weight_trg_pbpb_num(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
                    tnp_trig_weight_mumiL3_num = tnp_weight_trg_pbpb_num(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);

                    tnp_trig_weight_muplL2_den = tnp_weight_trg_pbpb_den(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
                    tnp_trig_weight_muplL3_den = tnp_weight_trg_pbpb_den(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
                    tnp_trig_weight_mumiL2_den = tnp_weight_trg_pbpb_den(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
                    tnp_trig_weight_mumiL3_den = tnp_weight_trg_pbpb_den(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);

                    tnp_trig_weight_num = tnp_trig_weight_muplL2_num * tnp_trig_weight_mumiL3_num + tnp_trig_weight_mumiL2_num * tnp_trig_weight_muplL3_num - tnp_trig_weight_muplL3_num * tnp_trig_weight_mumiL3_num;
                    tnp_trig_weight_den = tnp_trig_weight_muplL2_den * tnp_trig_weight_mumiL3_den + tnp_trig_weight_mumiL2_den * tnp_trig_weight_muplL3_den - tnp_trig_weight_muplL3_den * tnp_trig_weight_mumiL3_den;
                    tnp_trig_weight = tnp_trig_weight_num / tnp_trig_weight_den;
                }

                tnp_weight = tnp_weight * tnp_trig_weight;

                // if(SelDone == false){cout << "ERROR :: No muon filter combination selected !!!!" << endl; continue;}
                // if((tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1)){cout << "ERROR :: No trigger muon tnp scale factors selected !!!!" << endl; continue;}
                tnp_weight = tnp_weight * tnp_trig_weight_mupl * tnp_trig_weight_mumi;
                if (HLTPass == true && HLTFilterPass == true)
                {
                    counttnp++;
                    // tnp_trig_dimu = tnp_trig_weight_mupl * tnp_trig_weight_mumi;
                    tnp_trig_dimu = tnp_trig_weight;
                    hpt_tnp_trig_fwd->Fill(JP_Reco->Pt(), tnp_trig_dimu);
                }
            }

            pt_weight = 1;
            if (isPtWeight)
                pt_weight = fptw->Eval(JP_Reco->Pt());

            if (HLTPass == true && HLTFilterPass == true)
            {
                // if( Centrality > 0 && Centrality < 20 && JP_Reco->Pt() > 3) {
                //   if(Rapidity_reco < 0.6)				hpt_reco_1 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                //   else if(Rapidity_reco > 0.6 && Rapidity_reco < 1.2) 	hpt_reco_2 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                //   else if(Rapidity_reco > 1.2 && Rapidity_reco < 1.6) 	hpt_reco_3 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                //   else if(Rapidity_reco > 1.6 && Rapidity_reco < 2.4) 	hpt_reco_4 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                // }

                // if( Centrality > 20 && Centrality < 60 && JP_Reco->Pt() > 3) {
                //   if(Rapidity_reco < 0.6)				hpt_reco_5 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                //   else if(Rapidity_reco > 0.6 && Rapidity_reco < 1.2) 	hpt_reco_6 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                //   else if(Rapidity_reco > 1.2 && Rapidity_reco < 1.6) 	hpt_reco_7 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                //   else if(Rapidity_reco > 1.6 && Rapidity_reco < 2.4) 	hpt_reco_8 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                // }

                // if( Centrality > 60 && Centrality < 180 && JP_Reco->Pt() > 3) {
                //   if(Rapidity_reco < 0.6)				hpt_reco_9 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                //   else if(Rapidity_reco > 0.6 && Rapidity_reco < 1.2) 	hpt_reco_10 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                //   else if(Rapidity_reco > 1.2 && Rapidity_reco < 1.6) 	hpt_reco_11 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                //   else if(Rapidity_reco > 1.6 && Rapidity_reco < 2.4) 	hpt_reco_12 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
                // }
            }

            // ===== calculate polarization angles - Reco ===== //
            // polarization unit vectors in lab
            TVector3 uy_lab(0, 1, 0);
            TVector3 uz_lab(0, 0, 1);
            TVector3 uz_ep = uy_lab;

            // consider auto correleation
            float eta_ = JP_Reco->Eta();
            float epAng = 0;
            if (eta_ >= 0)
                epAng = epHFm2;
            else
                epAng = epHFp2;

            // get z-axis in EP frame
            uz_ep.Rotate(epAng, uz_lab); // uz_ep = y rotated as epAng in lab frame

            // Boost to Quarkonia rest frame
            std::unique_ptr<TLorentzVector> mupl_qrf = std::make_unique<TLorentzVector>(*mupl_Reco);
            mupl_qrf->Boost(-JP_Reco->BoostVector());
            std::unique_ptr<TLorentzVector> mumi_qrf = std::make_unique<TLorentzVector>(*mumi_Reco);
            mumi_qrf->Boost(-JP_Reco->BoostVector());
            std::unique_ptr<TLorentzVector> p1_qrf = std::make_unique<TLorentzVector>(*p1_lab);
            p1_qrf->Boost(-JP_Reco->BoostVector());
            std::unique_ptr<TLorentzVector> p2_qrf = std::make_unique<TLorentzVector>(*p2_lab);
            p2_qrf->Boost(-JP_Reco->BoostVector());

            // u_p1, u_p2
            TVector3 u_p1 = p1_qrf->Vect();
            u_p1 = u_p1.Unit();
            TVector3 u_p2 = p2_qrf->Vect();
            u_p2 = u_p2.Unit();

            // unit vector of mupl
            TVector3 u_mupl = mupl_qrf->Vect(); // get px, py, pz
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
            TVector3 uz_hx = JP_Reco->Vect();
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

            // EP
            TVector3 uy_ep = p1_3d_lab;
            uy_ep = uy_ep.Unit();
            TVector3 ux_ep = uy_ep.Cross(uz_ep);
            ux_ep = ux_ep.Unit();
            float cos_ep_ = u_mupl.Dot(uz_ep);
            float phi_ep_ = TMath::ATan2(u_mupl.Dot(uy_pol), u_mupl.Dot(ux_ep));

            // ===== fill the numerator ===== //
            if (Centrality > 00 && Centrality < 180 && JP_Reco->Pt() > 3)
            {
                c60_180_fwd_pt_num->Fill(JP_Gen->Pt(), weight * pt_weight);
                c60_180_fwd_cos_ep_num1->Fill(cos_ep_, weight * pt_weight);
                c60_180_fwd_cos_ep_num2->Fill(cos_ep_, weight);
                c60_180_fwd_num1->Fill(cos_ep_, JP_Gen->Pt(), weight * pt_weight);
                c60_180_fwd_num2->Fill(cos_ep_, JP_Gen->Pt(), weight);
            }
        } // end of reco loop

        // ===== printing info ===== //
        // do we need it?
        // if(iev%100000==0) cout << ">>>>> Pt_Weight " << pt_weight << " / " << ">>>>> TnP_Weight " <<  tnp_weight << endl;
        // if(pt_weight > 0.8 && pt_weight <1.2) countPtW++;
        // if (iev % 100000 == 0)
        //     cout << ">>>>> tnp_num " << tnp_trig_weight_num << " / " << ">>>>> tnp_den " << tnp_trig_weight_den << ">>>>> SF " << tnp_trig_weight << endl;
        // if (iev % 100000 == 0)
        //     cout << ">>>>> tnp_num_mL2 " << tnp_trig_weight_mumiL2_num << " / " << ">>>>> tnp_num_mL3 " << tnp_trig_weight_mumiL3_num << endl;
        // if (iev % 100000 == 0)
        //     cout << ">>>>> tnp_num_pL2 " << tnp_trig_weight_muplL2_num << " / " << ">>>>> tnp_num_pL3 " << tnp_trig_weight_muplL3_num << endl;
        // if (iev % 100000 == 0)
        //     cout << ">>>>> tnp_den_mL2 " << tnp_trig_weight_mumiL2_den << " / " << ">>>>> tnp_den_mL3 " << tnp_trig_weight_mumiL3_den << endl;
        // if (iev % 100000 == 0)
        //     cout << ">>>>> tnp_den_pL2 " << tnp_trig_weight_muplL2_den << " / " << ">>>>> tnp_den_pL3 " << tnp_trig_weight_muplL3_den << endl;
        // if (iev % 100000 == 0)
        //     cout << ">>>>> mupl_pt " << mupl_Reco->Pt() << " / " << ">>>>> mupl_eta " << mupl_Reco->Eta() << endl;
        // if (iev % 100000 == 0)
        //     cout << ">>>>> mumi_pt " << mumi_Reco->Pt() << " / " << ">>>>> mumi_eta " << mumi_Reco->Eta() << endl;
    } // end of main loop

    cout << "count " << count << endl;
    cout << "counttnp " << counttnp << endl;
    cout << "countptw " << countPtW << endl;

    // ===== compute efficiency ===== //
    // simplly divide num_hist / den_hist
    divide_2d_hist(c60_180_fwd_num1, c60_180_fwd_den1, eff_c60_180_fwd1);

    // ===== pT weight function fitting ===== //
    // 왜 있지?
    // 2D hist라 피팅 불가능.
    // f1->SetLineColor(kBlack);
    // f1->SetLineWidth(2);
    // eff_c0_20_fwd->Fit(f1);
    // eff_c0_20_mid->Fit(f1);
    // eff_c20_60_fwd->Fit(f1);
    // eff_c20_60_mid->Fit(f1);
    // eff_c60_180_fwd1->Fit(f1);
    // eff_c60_180_mid->Fit(f1);

    // ===== draw plots ===== //
    // gROOT->Macro("~/rootlogon.C");
    // gStyle->SetOptFit(0);

    //"|y|: 0.0-0.6, 0-10%"
    draw_2d_hist(eff_c60_180_fwd1, "|y|: 1.6-2.4, 60-90%", MCtype);
    draw_2d_hist(eff_c60_180_fwd2, "|y|: 1.6-2.4, 60-90%", MCtype);

    eff_c60_180_fwd1->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent60to180_absy1p6_2p4_1", isTnP, isPtWeight));
    eff_c60_180_fwd2->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent60to180_absy1p6_2p4_2", isTnP, isPtWeight));

    // ===== save output ===== //
    TString outFileName = Form("roots/mc_eff_vs_pt_cent_rap_prompt_pbpb_Jpsi_PtW%d_tnp%d_%s.root", isPtWeight, isTnP, date_label.Data());
    if (MCtype == 2)
        outFileName = Form("roots/mc_eff_vs_pt_cent_rap_nprompt_pbpb_Jpsi_PtW%d_tnp%d_%s.root", isPtWeight, isTnP, date_label.Data());
    TFile *outFile = new TFile(outFileName, "RECREATE");

    // 2d hist
    eff_c60_180_fwd1->Write();
    eff_c60_180_fwd2->Write();

    c60_180_fwd_num1->Write();
    c60_180_fwd_den1->Write();
    c60_180_fwd_num2->Write();
    c60_180_fwd_den2->Write();

    // 1d
    c60_180_fwd_pt_num->Write();
    c60_180_fwd_cos_ep_num1->Write();
    c60_180_fwd_cos_ep_num2->Write();
    c60_180_fwd_pt_den->Write();
    c60_180_fwd_cos_ep_den1->Write();
    c60_180_fwd_cos_ep_den2->Write();

    // f1->Write();
    if (isTnP)
        hpt_tnp_trig_fwd->Write();
    outFile->Close();

    cout << "\n=================================\n";
    cout << "\n Finish Acc(x)Eff Calculation\n";
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

    // TLegend legend_eff_1(0.6, 0.84);
    // legend_eff_1.AddEntry("h_eff", Form("%s",legend.c_str()), "lep");
    // legend_eff_1.SetBorderSize(0);
    // legend_eff_1.Draw("E");

    // TLatex lt1;
    // lt1.SetNDC();
    // lt1.SetTextSize(0.03);
    // lt1.SetTextSize(0.03);
    // if (MCtype == 1)
    //     lt1.DrawLatex(0.13, 0.84, "Prompt J/#psi (PbPb)");
    // else if (MCtype == 2)
    //     lt1.DrawLatex(0.13, 0.84, "Non-Prompt J/#psi (PbPb)");
}