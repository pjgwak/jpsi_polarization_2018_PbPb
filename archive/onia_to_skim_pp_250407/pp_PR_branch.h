//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  8 20:04:26 2025 by ROOT version 6.32.02
// from TTree myTree/My TTree of dimuons
// found on file: /disk1/Oniatree/Jpsi/OniaTreeMC_Jpsi_Pythia8_nonUL_5p02TeV_merged.root
//////////////////////////////////////////////////////////

#ifndef pp_PR_branch_h
#define pp_PR_branch_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"

class pp_PR_branch {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          eventNb;
   Float_t         zVtx;
   Short_t         nPV;
   Short_t         Ntracks;
   Int_t           trigPrescale[26];
   ULong64_t       HLTriggers;
   Short_t         Reco_QQ_size;
   Short_t         Reco_QQ_type[10];   //[Reco_QQ_size]
   Short_t         Reco_QQ_sign[10];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_4mom;
   Short_t         Reco_QQ_mupl_idx[10];   //[Reco_QQ_size]
   Short_t         Reco_QQ_mumi_idx[10];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_trig[10];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_isCowboy[10];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau[10];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr[10];   //[Reco_QQ_size]
   Float_t         Reco_QQ_cosAlpha[10];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau3D[10];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr3D[10];   //[Reco_QQ_size]
   Float_t         Reco_QQ_cosAlpha3D[10];   //[Reco_QQ_size]
   Short_t         Reco_QQ_whichGen[10];   //[Reco_QQ_size]
   Float_t         Reco_QQ_VtxProb[10];   //[Reco_QQ_size]
   Float_t         Reco_QQ_dca[10];   //[Reco_QQ_size]
   Float_t         Reco_QQ_MassErr[10];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_vtx;
   Short_t         Reco_mu_size;
   Short_t         Reco_mu_type[5];   //[Reco_mu_size]
   Short_t         Reco_mu_whichGen[5];   //[Reco_mu_size]
   Int_t           Reco_mu_SelectionType[5];   //[Reco_mu_size]
   Short_t         Reco_mu_charge[5];   //[Reco_mu_size]
   TClonesArray    *Reco_mu_4mom;
   ULong64_t       Reco_mu_trig[5];   //[Reco_mu_size]
   Bool_t          Reco_mu_isPF[5];   //[Reco_mu_size]
   Bool_t          Reco_mu_InTightAcc[5];   //[Reco_mu_size]
   Bool_t          Reco_mu_InLooseAcc[5];   //[Reco_mu_size]
   Bool_t          Reco_mu_highPurity[5];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixValHits[5];   //[Reco_mu_size]
   Int_t           Reco_mu_nMuValHits[5];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkHits[5];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_inner[5];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_global[5];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixWMea[5];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkWMea[5];   //[Reco_mu_size]
   Int_t           Reco_mu_StationsMatched[5];   //[Reco_mu_size]
   Float_t         Reco_mu_dxy[5];   //[Reco_mu_size]
   Float_t         Reco_mu_dxyErr[5];   //[Reco_mu_size]
   Float_t         Reco_mu_dz[5];   //[Reco_mu_size]
   Float_t         Reco_mu_dzErr[5];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_inner[5];   //[Reco_mu_size]
   Int_t           Reco_mu_simExtType[5];   //[Reco_mu_size]
   Float_t         Gen_weight;
   Float_t         Gen_pthat;
   Short_t         Gen_QQ_size;
   TClonesArray    *Gen_QQ_4mom;
   Float_t         Gen_QQ_ctau[3];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau3D[3];   //[Gen_QQ_size]
   Short_t         Gen_QQ_mupl_idx[3];   //[Gen_QQ_size]
   Short_t         Gen_QQ_mumi_idx[3];   //[Gen_QQ_size]
   Short_t         Gen_QQ_whichRec[3];   //[Gen_QQ_size]
   Int_t           Gen_QQ_momId[3];   //[Gen_QQ_size]
   Short_t         Gen_mu_size;
   Short_t         Gen_mu_charge[9];   //[Gen_mu_size]
   TClonesArray    *Gen_mu_4mom;
   Short_t         Gen_mu_whichRec[9];   //[Gen_mu_size]

   // List of branches
   TBranch        *b_eventNb;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_Ntracks;   //!
   TBranch        *b_trigPrescale;   //!
   TBranch        *b_HLTriggers;   //!
   TBranch        *b_Reco_QQ_size;   //!
   TBranch        *b_Reco_QQ_type;   //!
   TBranch        *b_Reco_QQ_sign;   //!
   TBranch        *b_Reco_QQ_4mom;   //!
   TBranch        *b_Reco_QQ_mupl_idx;   //!
   TBranch        *b_Reco_QQ_mumi_idx;   //!
   TBranch        *b_Reco_QQ_trig;   //!
   TBranch        *b_Reco_QQ_isCowboy;   //!
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
   TBranch        *b_Reco_mu_trig;   //!
   TBranch        *b_Reco_mu_isPF;   //!
   TBranch        *b_Reco_mu_InTightAcc;   //!
   TBranch        *b_Reco_mu_InLooseAcc;   //!
   TBranch        *b_Reco_mu_highPurity;   //!
   TBranch        *b_Reco_mu_nPixValHits;   //!
   TBranch        *b_Reco_mu_nMuValHits;   //!
   TBranch        *b_Reco_mu_nTrkHits;   //!
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
   TBranch        *b_Reco_mu_simExtType;   //!
   TBranch        *b_Gen_weight;   //!
   TBranch        *b_Gen_pthat;   //!
   TBranch        *b_Gen_QQ_size;   //!
   TBranch        *b_Gen_QQ_4mom;   //!
   TBranch        *b_Gen_QQ_ctau;   //!
   TBranch        *b_Gen_QQ_ctau3D;   //!
   TBranch        *b_Gen_QQ_mupl_idx;   //!
   TBranch        *b_Gen_QQ_mumi_idx;   //!
   TBranch        *b_Gen_QQ_whichRec;   //!
   TBranch        *b_Gen_QQ_momId;   //!
   TBranch        *b_Gen_mu_size;   //!
   TBranch        *b_Gen_mu_charge;   //!
   TBranch        *b_Gen_mu_4mom;   //!
   TBranch        *b_Gen_mu_whichRec;   //!

   pp_PR_branch(TTree *tree=0);
   virtual ~pp_PR_branch();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef pp_PR_branch_cxx
pp_PR_branch::pp_PR_branch(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/disk1/Oniatree/Jpsi/OniaTreeMC_Jpsi_Pythia8_nonUL_5p02TeV_merged.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/disk1/Oniatree/Jpsi/OniaTreeMC_Jpsi_Pythia8_nonUL_5p02TeV_merged.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/disk1/Oniatree/Jpsi/OniaTreeMC_Jpsi_Pythia8_nonUL_5p02TeV_merged.root:/hionia");
      dir->GetObject("myTree",tree);

   }
   Init(tree);
}

pp_PR_branch::~pp_PR_branch()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pp_PR_branch::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pp_PR_branch::LoadTree(Long64_t entry)
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

void pp_PR_branch::Init(TTree *tree)
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
   fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
   fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
   fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
   fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
   fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
   fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
   fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
   fChain->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
   fChain->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
   fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
   fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
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
   fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
   fChain->SetBranchAddress("Reco_mu_isPF", Reco_mu_isPF, &b_Reco_mu_isPF);
   fChain->SetBranchAddress("Reco_mu_InTightAcc", Reco_mu_InTightAcc, &b_Reco_mu_InTightAcc);
   fChain->SetBranchAddress("Reco_mu_InLooseAcc", Reco_mu_InLooseAcc, &b_Reco_mu_InLooseAcc);
   fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
   fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
   fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
   fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
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
   fChain->SetBranchAddress("Reco_mu_simExtType", Reco_mu_simExtType, &b_Reco_mu_simExtType);
   fChain->SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);
   fChain->SetBranchAddress("Gen_pthat", &Gen_pthat, &b_Gen_pthat);
   fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
   fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
   fChain->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
   fChain->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
   fChain->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
   fChain->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);
   fChain->SetBranchAddress("Gen_QQ_whichRec", Gen_QQ_whichRec, &b_Gen_QQ_whichRec);
   fChain->SetBranchAddress("Gen_QQ_momId", Gen_QQ_momId, &b_Gen_QQ_momId);
   fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
   fChain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
   fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
   fChain->SetBranchAddress("Gen_mu_whichRec", Gen_mu_whichRec, &b_Gen_mu_whichRec);
   Notify();
}

bool pp_PR_branch::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void pp_PR_branch::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pp_PR_branch::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef pp_PR_branch_cxx
