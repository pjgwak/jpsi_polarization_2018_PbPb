//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 12 15:32:50 2025 by ROOT version 6.32.02
// from TTree myTree/My TTree of dimuons
// found on file: /disk1/Oniatree/polarization/oniatree_5p36/2023_MC/Oniatree_MC_miniAOD_250612.root
//////////////////////////////////////////////////////////

#ifndef Run3_MC_branches_h
#define Run3_MC_branches_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "vector"

class Run3_MC_branches {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          eventNb;
   Float_t         zVtx;
   Short_t         nPV;
   Int_t           Centrality;
   Short_t         Npix;
   Short_t         NpixelTracks;
   Short_t         Ntracks;
   Int_t           trigPrescale[27];
   ULong64_t       HLTriggers;
   Float_t         SumET_HF;
   Float_t         SumET_HFplus;
   Float_t         SumET_HFminus;
   Float_t         SumET_HFplusEta4;
   Float_t         SumET_HFminusEta4;
   Float_t         SumET_ET;
   Float_t         SumET_EE;
   Float_t         SumET_EB;
   Float_t         SumET_EEplus;
   Float_t         SumET_EEminus;
   Float_t         SumET_ZDC;
   Float_t         SumET_ZDCplus;
   Float_t         SumET_ZDCminus;
   Int_t           nEP;
   Float_t         rpAng_origin[12];   //[nEP]
   Float_t         rpSin_origin[12];   //[nEP]
   Float_t         rpCos_origin[12];   //[nEP]
   Float_t         rpAng[12];   //[nEP]
   Float_t         rpSin[12];   //[nEP]
   Float_t         rpCos[12];   //[nEP]
   Short_t         Reco_QQ_size;
   Short_t         Reco_QQ_type[616];   //[Reco_QQ_size]
   Short_t         Reco_QQ_sign[616];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_4mom;
   Short_t         Reco_QQ_mupl_idx[616];   //[Reco_QQ_size]
   Short_t         Reco_QQ_mumi_idx[616];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_trig[616];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau[616];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr[616];   //[Reco_QQ_size]
   Float_t         Reco_QQ_cosAlpha[616];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau3D[616];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr3D[616];   //[Reco_QQ_size]
   Float_t         Reco_QQ_cosAlpha3D[616];   //[Reco_QQ_size]
   Short_t         Reco_QQ_whichGen[616];   //[Reco_QQ_size]
   Float_t         Reco_QQ_VtxProb[616];   //[Reco_QQ_size]
   Float_t         Reco_QQ_dca[616];   //[Reco_QQ_size]
   Float_t         Reco_QQ_MassErr[616];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_vtx;
   Short_t         Reco_mu_size;
   Short_t         Reco_mu_type[50];   //[Reco_mu_size]
   Short_t         Reco_mu_whichGen[50];   //[Reco_mu_size]
   Int_t           Reco_mu_SelectionType[50];   //[Reco_mu_size]
   Short_t         Reco_mu_charge[50];   //[Reco_mu_size]
   TClonesArray    *Reco_mu_4mom;
   TClonesArray    *Reco_mu_L1_4mom;
   ULong64_t       Reco_mu_trig[50];   //[Reco_mu_size]
   Bool_t          Reco_mu_InTightAcc[50];   //[Reco_mu_size]
   Bool_t          Reco_mu_InLooseAcc[50];   //[Reco_mu_size]
   Bool_t          Reco_mu_highPurity[50];   //[Reco_mu_size]
   Bool_t          Reco_mu_TMOneStaTight[50];   //[Reco_mu_size]
   Bool_t          Reco_mu_isPF[50];   //[Reco_mu_size]
   Bool_t          Reco_mu_isTracker[50];   //[Reco_mu_size]
   Bool_t          Reco_mu_isGlobal[50];   //[Reco_mu_size]
   Bool_t          Reco_mu_isSoftCutBased[50];   //[Reco_mu_size]
   Bool_t          Reco_mu_isHybridSoft[50];   //[Reco_mu_size]
   Bool_t          Reco_mu_isMediumCutBased[50];   //[Reco_mu_size]
   Bool_t          Reco_mu_isTightCutBased[50];   //[Reco_mu_size]
   Short_t         Reco_mu_candType[50];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixValHits[50];   //[Reco_mu_size]
   Int_t           Reco_mu_nMuValHits[50];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkHits[50];   //[Reco_mu_size]
   Float_t         Reco_mu_segmentComp[50];   //[Reco_mu_size]
   Float_t         Reco_mu_kink[50];   //[Reco_mu_size]
   Float_t         Reco_mu_localChi2[50];   //[Reco_mu_size]
   Float_t         Reco_mu_validFraction[50];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_bestTracker[50];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_inner[50];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_global[50];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixWMea[50];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkWMea[50];   //[Reco_mu_size]
   Int_t           Reco_mu_StationsMatched[50];   //[Reco_mu_size]
   Float_t         Reco_mu_dxy[50];   //[Reco_mu_size]
   Float_t         Reco_mu_dxyErr[50];   //[Reco_mu_size]
   Float_t         Reco_mu_dz[50];   //[Reco_mu_size]
   Float_t         Reco_mu_dzErr[50];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_inner[50];   //[Reco_mu_size]
   Float_t         Reco_mu_iso[50];   //[Reco_mu_size]
   Short_t         Reco_trk_size;
   vector<float>   *Reco_trk_4mom_pt;
   vector<float>   *Reco_trk_4mom_eta;
   vector<float>   *Reco_trk_4mom_phi;
   Float_t         Gen_weight;
   Float_t         Gen_pthat;
   Short_t         Gen_QQ_size;
   TClonesArray    *Gen_QQ_4mom;
   Float_t         Gen_QQ_ctau[3];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau3D[3];   //[Gen_QQ_size]
   Short_t         Gen_QQ_mupl_idx[3];   //[Gen_QQ_size]
   Short_t         Gen_QQ_mumi_idx[3];   //[Gen_QQ_size]
   Short_t         Gen_QQ_whichRec[3];   //[Gen_QQ_size]
   Short_t         Gen_mu_size;
   Short_t         Gen_mu_charge[20];   //[Gen_mu_size]
   TClonesArray    *Gen_mu_4mom;
   Short_t         Gen_mu_whichRec[20];   //[Gen_mu_size]

   // List of branches
   TBranch        *b_eventNb;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_Centrality;   //!
   TBranch        *b_Npix;   //!
   TBranch        *b_NpixelTracks;   //!
   TBranch        *b_Ntracks;   //!
   TBranch        *b_trigPrescale;   //!
   TBranch        *b_HLTriggers;   //!
   TBranch        *b_SumET_HF;   //!
   TBranch        *b_SumET_HFplus;   //!
   TBranch        *b_SumET_HFminus;   //!
   TBranch        *b_SumET_HFplusEta4;   //!
   TBranch        *b_SumET_HFminusEta4;   //!
   TBranch        *b_SumET_ET;   //!
   TBranch        *b_SumET_EE;   //!
   TBranch        *b_SumET_EB;   //!
   TBranch        *b_SumET_EEplus;   //!
   TBranch        *b_SumET_EEminus;   //!
   TBranch        *b_SumET_ZDC;   //!
   TBranch        *b_SumET_ZDCplus;   //!
   TBranch        *b_SumET_ZDCminus;   //!
   TBranch        *b_nEP;   //!
   TBranch        *b_rpAng_origin;   //!
   TBranch        *b_rpSin_origin;   //!
   TBranch        *b_rpCos_origin;   //!
   TBranch        *b_rpAng;   //!
   TBranch        *b_rpSin;   //!
   TBranch        *b_rpCos;   //!
   TBranch        *b_Reco_QQ_size;   //!
   TBranch        *b_Reco_QQ_type;   //!
   TBranch        *b_Reco_QQ_sign;   //!
   TBranch        *b_Reco_QQ_4mom;   //!
   TBranch        *b_Reco_QQ_mupl_idx;   //!
   TBranch        *b_Reco_QQ_mumi_idx;   //!
   TBranch        *b_Reco_QQ_trig;   //!
   TBranch        *b_Reco_QQ_ctau;   //!
   TBranch        *b_Reco_QQ_ctauErr;   //!
   TBranch        *b_Reco_QQ_cosAlpha;   //!
   TBranch        *b_Reco_QQ_ctau3D;   //!
   TBranch        *b_Reco_QQ_ctauErr3D;   //!
   TBranch        *b_Reco_QQ_cosAlpha3D;   //!
   TBranch        *b_Reco_QQ_whichGen;   //!
   TBranch        *b_Reco_QQ_VtxProb;   //!
   TBranch        *b_Reco_QQ_dca;   //!
   TBranch        *b_Reco_QQ_MassErr;   //!
   TBranch        *b_Reco_QQ_vtx;   //!
   TBranch        *b_Reco_mu_size;   //!
   TBranch        *b_Reco_mu_type;   //!
   TBranch        *b_Reco_mu_whichGen;   //!
   TBranch        *b_Reco_mu_SelectionType;   //!
   TBranch        *b_Reco_mu_charge;   //!
   TBranch        *b_Reco_mu_4mom;   //!
   TBranch        *b_Reco_mu_L1_4mom;   //!
   TBranch        *b_Reco_mu_trig;   //!
   TBranch        *b_Reco_mu_InTightAcc;   //!
   TBranch        *b_Reco_mu_InLooseAcc;   //!
   TBranch        *b_Reco_mu_highPurity;   //!
   TBranch        *b_Reco_mu_TMOneStaTight;   //!
   TBranch        *b_Reco_mu_isPF;   //!
   TBranch        *b_Reco_mu_isTracker;   //!
   TBranch        *b_Reco_mu_isGlobal;   //!
   TBranch        *b_Reco_mu_isSoftCutBased;   //!
   TBranch        *b_Reco_mu_isHybridSoft;   //!
   TBranch        *b_Reco_mu_isMediumCutBased;   //!
   TBranch        *b_Reco_mu_isTightCutBased;   //!
   TBranch        *b_Reco_mu_candType;   //!
   TBranch        *b_Reco_mu_nPixValHits;   //!
   TBranch        *b_Reco_mu_nMuValHits;   //!
   TBranch        *b_Reco_mu_nTrkHits;   //!
   TBranch        *b_Reco_mu_segmentComp;   //!
   TBranch        *b_Reco_mu_kink;   //!
   TBranch        *b_Reco_mu_localChi2;   //!
   TBranch        *b_Reco_mu_validFraction;   //!
   TBranch        *b_Reco_mu_normChi2_bestTracker;   //!
   TBranch        *b_Reco_mu_normChi2_inner;   //!
   TBranch        *b_Reco_mu_normChi2_global;   //!
   TBranch        *b_Reco_mu_nPixWMea;   //!
   TBranch        *b_Reco_mu_nTrkWMea;   //!
   TBranch        *b_Reco_mu_StationsMatched;   //!
   TBranch        *b_Reco_mu_dxy;   //!
   TBranch        *b_Reco_mu_dxyErr;   //!
   TBranch        *b_Reco_mu_dz;   //!
   TBranch        *b_Reco_mu_dzErr;   //!
   TBranch        *b_Reco_mu_ptErr_inner;   //!
   TBranch        *b_Reco_mu_iso;   //!
   TBranch        *b_Reco_trk_size;   //!
   TBranch        *b_Reco_trk_4mom_pt;   //!
   TBranch        *b_Reco_trk_4mom_eta;   //!
   TBranch        *b_Reco_trk_4mom_phi;   //!
   TBranch        *b_Gen_weight;   //!
   TBranch        *b_Gen_pthat;   //!
   TBranch        *b_Gen_QQ_size;   //!
   TBranch        *b_Gen_QQ_4mom;   //!
   TBranch        *b_Gen_QQ_ctau;   //!
   TBranch        *b_Gen_QQ_ctau3D;   //!
   TBranch        *b_Gen_QQ_mupl_idx;   //!
   TBranch        *b_Gen_QQ_mumi_idx;   //!
   TBranch        *b_Gen_QQ_whichRec;   //!
   TBranch        *b_Gen_mu_size;   //!
   TBranch        *b_Gen_mu_charge;   //!
   TBranch        *b_Gen_mu_4mom;   //!
   TBranch        *b_Gen_mu_whichRec;   //!

   Run3_MC_branches(TTree *tree=0);
   virtual ~Run3_MC_branches();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Run3_MC_branches_cxx
Run3_MC_branches::Run3_MC_branches(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/disk1/Oniatree/polarization/oniatree_5p36/2023_MC/Oniatree_MC_miniAOD_250612.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/disk1/Oniatree/polarization/oniatree_5p36/2023_MC/Oniatree_MC_miniAOD_250612.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/disk1/Oniatree/polarization/oniatree_5p36/2023_MC/Oniatree_MC_miniAOD_250612.root:/hionia");
      dir->GetObject("myTree",tree);

   }
   Init(tree);
}

Run3_MC_branches::~Run3_MC_branches()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Run3_MC_branches::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Run3_MC_branches::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Run3_MC_branches::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Reco_QQ_4mom = 0;
   Reco_QQ_vtx = 0;
   Reco_mu_4mom = 0;
   Reco_mu_L1_4mom = 0;
   Reco_trk_4mom_pt = 0;
   Reco_trk_4mom_eta = 0;
   Reco_trk_4mom_phi = 0;
   Gen_QQ_4mom = 0;
   Gen_mu_4mom = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
   fChain->SetBranchAddress("Npix", &Npix, &b_Npix);
   fChain->SetBranchAddress("NpixelTracks", &NpixelTracks, &b_NpixelTracks);
   fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
   fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
   fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
   fChain->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
   fChain->SetBranchAddress("SumET_HFplus", &SumET_HFplus, &b_SumET_HFplus);
   fChain->SetBranchAddress("SumET_HFminus", &SumET_HFminus, &b_SumET_HFminus);
   fChain->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4, &b_SumET_HFplusEta4);
   fChain->SetBranchAddress("SumET_HFminusEta4", &SumET_HFminusEta4, &b_SumET_HFminusEta4);
   fChain->SetBranchAddress("SumET_ET", &SumET_ET, &b_SumET_ET);
   fChain->SetBranchAddress("SumET_EE", &SumET_EE, &b_SumET_EE);
   fChain->SetBranchAddress("SumET_EB", &SumET_EB, &b_SumET_EB);
   fChain->SetBranchAddress("SumET_EEplus", &SumET_EEplus, &b_SumET_EEplus);
   fChain->SetBranchAddress("SumET_EEminus", &SumET_EEminus, &b_SumET_EEminus);
   fChain->SetBranchAddress("SumET_ZDC", &SumET_ZDC, &b_SumET_ZDC);
   fChain->SetBranchAddress("SumET_ZDCplus", &SumET_ZDCplus, &b_SumET_ZDCplus);
   fChain->SetBranchAddress("SumET_ZDCminus", &SumET_ZDCminus, &b_SumET_ZDCminus);
   fChain->SetBranchAddress("nEP", &nEP, &b_nEP);
   fChain->SetBranchAddress("rpAng_origin", rpAng_origin, &b_rpAng_origin);
   fChain->SetBranchAddress("rpSin_origin", rpSin_origin, &b_rpSin_origin);
   fChain->SetBranchAddress("rpCos_origin", rpCos_origin, &b_rpCos_origin);
   fChain->SetBranchAddress("rpAng", rpAng, &b_rpAng);
   fChain->SetBranchAddress("rpSin", rpSin, &b_rpSin);
   fChain->SetBranchAddress("rpCos", rpCos, &b_rpCos);
   fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
   fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
   fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
   fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
   fChain->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
   fChain->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
   fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
   fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
   fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
   fChain->SetBranchAddress("Reco_QQ_cosAlpha", Reco_QQ_cosAlpha, &b_Reco_QQ_cosAlpha);
   fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
   fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
   fChain->SetBranchAddress("Reco_QQ_cosAlpha3D", Reco_QQ_cosAlpha3D, &b_Reco_QQ_cosAlpha3D);
   fChain->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen, &b_Reco_QQ_whichGen);
   fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
   fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
   fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
   fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
   fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
   fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
   fChain->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
   fChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
   fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
   fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
   fChain->SetBranchAddress("Reco_mu_L1_4mom", &Reco_mu_L1_4mom, &b_Reco_mu_L1_4mom);
   fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
   fChain->SetBranchAddress("Reco_mu_InTightAcc", Reco_mu_InTightAcc, &b_Reco_mu_InTightAcc);
   fChain->SetBranchAddress("Reco_mu_InLooseAcc", Reco_mu_InLooseAcc, &b_Reco_mu_InLooseAcc);
   fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
   fChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
   fChain->SetBranchAddress("Reco_mu_isPF", Reco_mu_isPF, &b_Reco_mu_isPF);
   fChain->SetBranchAddress("Reco_mu_isTracker", Reco_mu_isTracker, &b_Reco_mu_isTracker);
   fChain->SetBranchAddress("Reco_mu_isGlobal", Reco_mu_isGlobal, &b_Reco_mu_isGlobal);
   fChain->SetBranchAddress("Reco_mu_isSoftCutBased", Reco_mu_isSoftCutBased, &b_Reco_mu_isSoftCutBased);
   fChain->SetBranchAddress("Reco_mu_isHybridSoft", Reco_mu_isHybridSoft, &b_Reco_mu_isHybridSoft);
   fChain->SetBranchAddress("Reco_mu_isMediumCutBased", Reco_mu_isMediumCutBased, &b_Reco_mu_isMediumCutBased);
   fChain->SetBranchAddress("Reco_mu_isTightCutBased", Reco_mu_isTightCutBased, &b_Reco_mu_isTightCutBased);
   fChain->SetBranchAddress("Reco_mu_candType", Reco_mu_candType, &b_Reco_mu_candType);
   fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
   fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
   fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
   fChain->SetBranchAddress("Reco_mu_segmentComp", Reco_mu_segmentComp, &b_Reco_mu_segmentComp);
   fChain->SetBranchAddress("Reco_mu_kink", Reco_mu_kink, &b_Reco_mu_kink);
   fChain->SetBranchAddress("Reco_mu_localChi2", Reco_mu_localChi2, &b_Reco_mu_localChi2);
   fChain->SetBranchAddress("Reco_mu_validFraction", Reco_mu_validFraction, &b_Reco_mu_validFraction);
   fChain->SetBranchAddress("Reco_mu_normChi2_bestTracker", Reco_mu_normChi2_bestTracker, &b_Reco_mu_normChi2_bestTracker);
   fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
   fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
   fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
   fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
   fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
   fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
   fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
   fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
   fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
   fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
   fChain->SetBranchAddress("Reco_mu_iso", Reco_mu_iso, &b_Reco_mu_iso);
   fChain->SetBranchAddress("Reco_trk_size", &Reco_trk_size, &b_Reco_trk_size);
   fChain->SetBranchAddress("Reco_trk_4mom_pt", &Reco_trk_4mom_pt, &b_Reco_trk_4mom_pt);
   fChain->SetBranchAddress("Reco_trk_4mom_eta", &Reco_trk_4mom_eta, &b_Reco_trk_4mom_eta);
   fChain->SetBranchAddress("Reco_trk_4mom_phi", &Reco_trk_4mom_phi, &b_Reco_trk_4mom_phi);
   fChain->SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);
   fChain->SetBranchAddress("Gen_pthat", &Gen_pthat, &b_Gen_pthat);
   fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
   fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
   fChain->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
   fChain->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
   fChain->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
   fChain->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);
   fChain->SetBranchAddress("Gen_QQ_whichRec", Gen_QQ_whichRec, &b_Gen_QQ_whichRec);
   fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
   fChain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
   fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
   fChain->SetBranchAddress("Gen_mu_whichRec", Gen_mu_whichRec, &b_Gen_mu_whichRec);
   Notify();
}

bool Run3_MC_branches::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void Run3_MC_branches::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Run3_MC_branches::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Run3_MC_branches_cxx
