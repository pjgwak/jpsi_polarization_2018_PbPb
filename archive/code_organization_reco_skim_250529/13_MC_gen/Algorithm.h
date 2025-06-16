#pragma once

#include "TString.h"
#include "TH1D.h"
#include "TTree.h"
#include "Data.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class Algorithm {
  public:
    /**
     * @brief construct a Alogrithm instance
     */
    Algorithm(bool isGen, bool isGenOnly);

    /**
     * @brief initialize the algorithm
     * @param data - pointer to the instance of Data class
     */
    void initialize(Data* data, bool isRun2);

    // excute main logics
    void execute(Long64_t irqq);
    void executeGen(Long64_t igqq);

    TH1D* h_runNb_m = nullptr; // for test
    TTree* m_myTree = nullptr;

    // turn on(off) selection cuts
    bool cut_centrality0_180 = false;
    bool cut_jpsiMass = false;
    bool cut_jpsiRapidity = false;
    bool cut_softMuons = false;
    bool cut_vtxProbability = false;
    bool cut_oppositeSign = false;
    bool cut_singleMuonAcc = false;
    bool cut_HLTriggerPbPbJpsi2018 = false;
    bool cut_recoQQTrigger = false;
    bool cut_L2L3FilterPbPbJpsi2018 = false;
    bool cut_runNb327123 = false;
    bool cut_whichGen = false;
    bool cut_tnp = false; // check L2 and L3 filters
    bool cut_muChargeGen = false;

    /**
     * @brief set centrality bin according to HF's ET
     */
    Int_t getHiBinFromhiHF(const Double_t hiHF);

    /**
     * @brief find number of collisions from the table
     */
    Double_t findNcoll(int hiBin);

    /**
     * @brief methods to make the outputs
     */
    void fillTree();

  protected:
    /**
     * Physics trigger numbers
     */
    int kTrigJpsi2018 = 12;
    int kL2filterJpsi2018 = 16;
    int kL3filterJpsi2018 = 17;

    int kTrigJpsipp2015 = 3;
    

    /**
     * @brief apply the selection
     */
    bool passSelection(Long64_t irqq);
    bool passSoftMuonCut(Long64_t irqq);
    bool passMuonAcc2018();
    bool checkMuonAcc2018(double muPt, double muEta);
    bool passRunNb327123();
    bool passHLTriggerJpsiPbPb2018();
    bool passRecoQQTrigger(Long64_t irqq);
    bool passHLFilterJpsiPbPb2018(Long64_t irqq);
    bool passewhichGen(Long64_t irqq);
    bool passTnpLogic(Long64_t irqq);

    bool passSelectionGen(Long64_t igqq);
    bool passMuChargeGen(Long64_t igqq);
    bool passMuonAcc2018Gen();

    /**
     * @brief fill the histograms and TTree
     */
    void standByFilling(Long64_t irqq);
    void standByFillingGen(Long64_t igqq);
    void fillHists();

    /**
     * @brief variables for polarization
     */
    TVector3 muplHX;
    TVector3 muplCS;
    TVector3 muplEP;
    Double_t m_cosEP;
    Double_t m_phiEP;

    /**
     * @brief polarization transform functions
     */
    void MuPlusVector_Helicity(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab);
    void MuPlusVector_CollinsSoper(const TLorentzVector &QQLV_Lab);
    void MuPlusVector_EventPlane(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab);

    /**
     * @brief instance of the Data class
     * and run2 sample flag
     */
    Data *m_data = nullptr;
    bool m_isRun2 = false;
    bool m_isGen = false;
    bool m_isGenOnly = false;

    /**
     * @brief set output branches
     */
    void setBranches();
    void setBranchesGen();

    /**
     * @brief variables of m_myTree
     */
    int evt;
    int runN;
    int lumi;
    int cBin;
    int nDimu;
    float vz;
    float mass[1500];
    float pt[1500];
    float pt1[1500];
    float pt2[1500];
    float y[1500];
    float phi[1500];
    float phi1[1500];
    float phi2[1500];
    float eta[1500];
    float eta1[1500];
    float eta2[1500];
    float weight0[1500];
    float weight1[1500];
    int recoQQsign[1500];
    float ctau3D[1500];
    float ctau3DErr[1500];
    float ctau3DTrue[1500];
    float ctau3DRes[1500];
    float cosHX[1500];
    float phiHX[1500];
    float cosCS[1500];
    float phiCS[1500];
    float cosEP[1500];
    float phiEP[1500];
    double weight = 1; // HI NColl x Gen weight
    double TnPweight[1500] = {1.}; // HI NColl x Gen weight x TnP 

    /**
     * @brief memober variables to store TnP weights
     */
    Double_t m_tnpWeight = 1;

    /**
     * @brief number of elements of centrality bins
     */
    const Int_t nBins = 200; // table of bin edges

    /**
     * @brief heavy ion centrality table - 2018
     * is 2023 data use same table?
     * I wrote 200 + 1 instaed of 201 to remember
     * that the range of centrality is 0 to 200 (0 ~ 100 %)
     */
    const Double_t binTable[200 + 1] = {0, 10.5072, 11.2099, 11.8364, 12.478, 13.1194, 13.7623, 14.4081, 15.0709, 15.7532, 16.4673, 17.1881, 17.923, 18.673, 19.4865, 20.3033, 21.1536, 22.0086, 22.9046, 23.8196, 24.7924, 25.8082, 26.8714, 27.9481, 29.0828, 30.2757, 31.5043, 32.8044, 34.1572, 35.6142, 37.1211, 38.6798, 40.3116, 42.0398, 43.8572, 45.6977, 47.6312, 49.6899, 51.815, 54.028, 56.3037, 58.7091, 61.2024, 63.8353, 66.5926, 69.3617, 72.2068, 75.2459, 78.3873, 81.5916, 84.9419, 88.498, 92.1789, 95.9582, 99.8431, 103.739, 107.78, 111.97, 116.312, 120.806, 125.46, 130.269, 135.247, 140.389, 145.713, 151.212, 156.871, 162.729, 168.762, 174.998, 181.424, 188.063, 194.907, 201.942, 209.19, 216.683, 224.37, 232.291, 240.43, 248.807, 257.416, 266.256, 275.348, 284.668, 294.216, 304.053, 314.142, 324.488, 335.101, 345.974, 357.116, 368.547, 380.283, 392.29, 404.564, 417.122, 429.968, 443.116, 456.577, 470.357, 484.422, 498.78, 513.473, 528.479, 543.813, 559.445, 575.411, 591.724, 608.352, 625.344, 642.686, 660.361, 678.371, 696.749, 715.485, 734.608, 754.068, 773.846, 794.046, 814.649, 835.608, 856.972, 878.719, 900.887, 923.409, 946.374, 969.674, 993.435, 1017.62, 1042.21, 1067.28, 1092.72, 1118.64, 1144.96, 1171.71, 1198.98, 1226.67, 1254.82, 1283.46, 1312.65, 1342.21, 1372.27, 1402.85, 1433.93, 1465.49, 1497.62, 1530.29, 1563.49, 1597.22, 1631.49, 1666.37, 1701.8, 1737.75, 1774.35, 1811.51, 1849.29, 1887.75, 1926.79, 1966.6, 2006.97, 2047.99, 2089.71, 2132.1, 2175.23, 2219.17, 2263.72, 2309.2, 2355.43, 2402.47, 2450.33, 2499.05, 2548.66, 2599.16, 2650.59, 2703.03, 2756.32, 2810.75, 2866.27, 2922.91, 2980.54, 3039.47, 3099.53, 3160.98, 3223.66, 3287.71, 3353.18, 3420.34, 3489.13, 3559.72, 3632.06, 3706.18, 3782.42, 3860.78, 3941.42, 4024.52, 4110.27, 4199.4, 4292.8, 4394.49, 4519.52, 5199.95};  


    /**
     * TnP tables - copied from tnp_weight_lowptPbPb.h
     */
    double tnp_weight_muid_pbpb(double pt, double eta, int idx = 0);
    double tnp_weight_trk_pbpb(double eta, int idx = 0);
    double tnp_weight_trg_pbpb(double pt, double eta, int filterId = 0, int idx = 0, bool getRawDen = false);
};