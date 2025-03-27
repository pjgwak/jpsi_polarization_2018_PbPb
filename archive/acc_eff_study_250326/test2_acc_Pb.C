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

void test2_acc_Pb(string data_label_ = "test2", int nevt = 100000, int MCtype = 1, int cLow = 0, int cHigh = 180, bool isTnP = true, bool isPtWeight = true)
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
    vector<double> phi_bin = {-5.0, -4, -3, -2, -1, 0.0, 1, 2, 3, 4, 5};

    double n_fwd_pt = fwd_pt_bin.size() - 1;
    double n_mid_pt = mid_pt_bin.size() - 1;
    double n_cos = cos_bin.size() - 1;
    double n_phi = phi_bin.size() - 1;

    // 2d histograms
    auto c60_180_fwd_num1 = std::make_shared<TH2D>("c60_180_fwd_num1", "", n_cos, &cos_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto c60_180_fwd_den1 = std::make_shared<TH2D>("c60_180_fwd_den1", "", n_cos, &cos_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto eff_c60_180_fwd_ep = std::make_shared<TH2D>("eff_c60_180_fwd_ep", "", n_cos, &cos_bin[0], n_fwd_pt, &fwd_pt_bin[0]);

    // 3d histograms
    auto c60_180_fwd_ep_num = std::make_shared<TH3D>("c60_180_fwd_ep_num", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto c60_180_fwd_ep_den = std::make_shared<TH3D>("c60_180_fwd_ep_den", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);

    auto c60_180_fwd_hx_num = std::make_shared<TH3D>("c60_180_fwd_hx_num", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto c60_180_fwd_hx_den = std::make_shared<TH3D>("c60_180_fwd_hx_den", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);

    auto c60_180_fwd_cs_num = std::make_shared<TH3D>("c60_180_fwd_cs_num", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);
    auto c60_180_fwd_cs_den = std::make_shared<TH3D>("c60_180_fwd_cs_den", "", n_cos, &cos_bin[0], n_phi, &phi_bin[0], n_fwd_pt, &fwd_pt_bin[0]);

    // =====  make hists to save error ===== //
    // 2d
    c60_180_fwd_num1->Sumw2();
    c60_180_fwd_den1->Sumw2();
    eff_c60_180_fwd_ep->Sumw2();

    // 3d
    c60_180_fwd_ep_num->Sumw2();
    c60_180_fwd_ep_den->Sumw2();
    c60_180_fwd_hx_num->Sumw2();
    c60_180_fwd_hx_den->Sumw2();
    c60_180_fwd_cs_num->Sumw2();
    c60_180_fwd_cs_den->Sumw2();

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


            // ===== Acc denominator cut ===== //
            if (!(fabs(JP_Gen->Rapidity()) < 2.4))
                continue;

            // ===== calculate polarization angles - Gen ===== //
            // polarization unit vectors in lab
            TVector3 uy_lab(0, 1, 0);
            TVector3 uz_lab(0, 0, 1);
            TVector3 uz_ep = uy_lab;

            // consider auto correleation - 어떻게 정하지? 일단 Reco랑 같게 -> 필요하면 바꿔서 측정
            // 그런데 epHFm2랑 epHFp2가 크게 다르지는 않을 것 같기는 한데...
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
                c60_180_fwd_den1->Fill(cos_ep_, JP_Gen->Pt(), weight * pt_weight);
                c60_180_fwd_ep_den->Fill(cos_ep_, phi_ep_, JP_Gen->Pt(), weight * pt_weight);
                c60_180_fwd_hx_den->Fill(cos_hx_, phi_hx_, JP_Gen->Pt(), weight * pt_weight);
                c60_180_fwd_cs_den->Fill(cos_cs_, phi_cs_, JP_Gen->Pt(), weight * pt_weight);
            }


            // ===== numerator cuts ===== //
            if (!(JP_Gen->Pt() > 0 && JP_Gen->Pt() < 50 && fabs(JP_Gen->Rapidity()) < 2.4 && IsAcceptable(mupl_Gen->Pt(), fabs(mupl_Gen->Eta())) && IsAcceptable(mumi_Gen->Pt(), fabs(mumi_Gen->Eta()))))
                continue;

            if (!(fabs(JP_Gen->Rapidity()) < 2.4 && fabs(mupl_Gen->Eta()) < 2.4 && fabs(mumi_Gen->Eta()) < 2.4))
                continue;

            if (Gen_mu_charge[Gen_QQ_mupl_idx[igen]] * Gen_mu_charge[Gen_QQ_mumi_idx[igen]] > 0)
                continue;

            // ===== fill the numerator ===== //
            if (Centrality > 0 && Centrality < 180 && JP_Gen->Pt() > 3)
            {
                c60_180_fwd_num1->Fill(cos_ep_, JP_Gen->Pt(), weight * pt_weight);
                c60_180_fwd_ep_num->Fill(cos_ep_, phi_ep_, JP_Gen->Pt(), weight * pt_weight);
                c60_180_fwd_hx_num->Fill(cos_hx_, phi_hx_, JP_Gen->Pt(), weight * pt_weight);
                c60_180_fwd_cs_num->Fill(cos_cs_, phi_cs_, JP_Gen->Pt(), weight * pt_weight);
            }
        } // end of gen loop
    } // end of main loop

    cout << "count " << count << endl;
    cout << "counttnp " << counttnp << endl;
    cout << "countptw " << countPtW << endl;

    // ===== compute efficiency ===== //
    // simplly divide num_hist / den_hist
    divide_2d_hist(c60_180_fwd_num1, c60_180_fwd_den1, eff_c60_180_fwd_ep);

    // ===== pT weight function fitting ===== //
    // 왜 있지?
    // 2D hist라 피팅 불가능.
    // f1->SetLineColor(kBlack);
    // f1->SetLineWidth(2);
    // eff_c0_20_fwd->Fit(f1);
    // eff_c0_20_mid->Fit(f1);
    // eff_c20_60_fwd->Fit(f1);
    // eff_c20_60_mid->Fit(f1);
    // eff_c60_180_fwd_ep->Fit(f1);
    // eff_c60_180_mid->Fit(f1);

    // ===== draw plots ===== //
    // gROOT->Macro("~/rootlogon.C");
    // gStyle->SetOptFit(0);

    //"|y|: 0.0-0.6, 0-10%"
    draw_2d_hist(eff_c60_180_fwd_ep, "|y|: 1.6-2.4, 60-90%", MCtype);

    eff_c60_180_fwd_ep->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent60to180_absy1p6_2p4_1", isTnP, isPtWeight));

    // ===== save output ===== //
    TString outFileName = Form("roots/mc_eff_vs_pt_cent_rap_prompt_pbpb_Jpsi_PtW%d_tnp%d_%s.root", isPtWeight, isTnP, date_label.Data());
    if (MCtype == 2)
        outFileName = Form("roots/mc_eff_vs_pt_cent_rap_nprompt_pbpb_Jpsi_PtW%d_tnp%d_%s.root", isPtWeight, isTnP, date_label.Data());
    TFile *outFile = new TFile(outFileName, "RECREATE");

    // 2d hist
    eff_c60_180_fwd_ep->Write();
    c60_180_fwd_num1->Write();
    c60_180_fwd_den1->Write();

    // 3d
    c60_180_fwd_ep_num->Write();
    c60_180_fwd_ep_den->Write();
    c60_180_fwd_hx_num->Write();
    c60_180_fwd_hx_den->Write();
    c60_180_fwd_cs_num->Write();
    c60_180_fwd_cs_den->Write();
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

bool IsAcceptable(double pt, double eta)
{
    if ((fabs(eta) < 1.2 && pt > 3.5) || ((fabs(eta) >= 1.2 && fabs(eta) < 2.1) && pt > (5.47 - 1.89 * fabs(eta))) || ((fabs(eta) >= 2.1 && fabs(eta) < 2.4) && pt > 1.5))
        return true;
    //          if ((fabs(eta) < 1.0 && pt>3.4) || ((fabs(eta) >=1.0 && fabs(eta) <1.6) && pt>(5.8-2.4*fabs(eta))) || ((fabs(eta) >=1.6 && fabs(eta) <2.4) && pt>3.3667-7.0/9.0*fabs(eta))) return true;//2010 version for 2.76 TeV

    else
        return false;
}