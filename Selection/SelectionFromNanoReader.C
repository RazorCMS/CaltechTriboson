#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>
#include <TMath.h>                  // ROOT math library
#include <TRandom3.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <sstream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <TChain.h>
#include <TH1D.h>
#include <TLorentzVector.h>

#include "NanoReader.hh"

#endif

struct LepCand {
  TLorentzVector v;
  int pid;
};

void SelectionFromNanoReader(TString infile, TString outfile) {

  const float EL_MASS = 0.000511;
  const float MU_MASS = 0.105658;
  const float Z_MASS = 91.1876;

  const float LEP_ETA_CUT = 2.4;
  const float LEP_PT_CUT = 10;
  const float LEADLEP_PT_CUT = 30;

  TChain *ch = new TChain("Events");
  ch->Add(infile);

  NanoReader nr(ch);
  int nEntries = nr.GetEntries(); 
  
  TFile *of = new TFile(outfile,"recreate");
  TTree *ot = new TTree("wwz","wwz");

  TH1D *hwt = new TH1D("hwt","hwt",1,0.5,1.5);

  uint run, ls, evt;
  float genWeight;
  int npu, npv;

  uint nE=0, nM=0;
  uint nJ15=0;

  LepCand l1, l2, l3, l4;

  float metPt=0, metPhi=0;
  float puppiPt=0, puppiPhi=0;

  ot->Branch("run", &run, "run/i");
  ot->Branch("ls", &ls, "ls/i");
  ot->Branch("evt", &evt, "evt/i");

  ot->Branch("genWeight", &genWeight, "genWeight/F");
  ot->Branch("npu", &npu, "npu/I");
  ot->Branch("npv", &npv, "npv/I");

  ot->Branch("nE", &nE, "nE/i");
  ot->Branch("nM", &nM, "nM/i");
  ot->Branch("nJ15", &nJ15, "nJ15/i");

  ot->Branch("lep1", &l1);
  ot->Branch("lep2", &l2);
  ot->Branch("lep3", &l3);
  ot->Branch("lep4", &l4);

  ot->Branch("metPt", &metPt);
  ot->Branch("metPhi", &metPhi);

  ot->Branch("puppiPt", &puppiPt);
  ot->Branch("puppiPhi", &puppiPhi);

  for (int i=0; i<nEntries; i++) {
    //for (int i=0; i<100; i++) {
    nr.LoadTree(i); nr.GetEntry(i);

    if (i%10000==0) std::cout << "Processing event " << i << " of " << nEntries << std::endl;

    hwt->SetBinContent(1, nr.genWeight+hwt->GetBinContent(1));

    l1.pid=0; l2.pid=0; l3.pid=0; l4.pid=0;
    l1.v.SetPtEtaPhiM(0,0,0,0); l2.v.SetPtEtaPhiM(0,0,0,0);
    l3.v.SetPtEtaPhiM(0,0,0,0); l4.v.SetPtEtaPhiM(0,0,0,0);

    //trigger requirements
    if (!nr.HLT_Ele27_WPTight_Gsf && !nr.HLT_IsoMu22) continue;

    nM=0;
    for (uint im=0; im<nr.nMuon; im++) {
      if (nr.Muon_pt[im]<LEP_PT_CUT || abs(nr.Muon_eta[im])>LEP_ETA_CUT) continue; //acceptance
      if (!nr.Muon_mediumId[im] && !nr.Muon_softId[im]) continue; //muon ID
      
      if (nr.Muon_miniPFRelIso_all[im]>0.2) continue; //muon isolation

      if (nr.Muon_pt[im]>l1.v.Pt()) {
	l4.v = l3.v; l4.pid=l3.pid;
	l3.v = l2.v; l3.pid=l2.pid;
	l2.v = l1.v; l2.pid=l1.pid;
	l1.v.SetPtEtaPhiM(nr.Muon_pt[im], nr.Muon_eta[im], nr.Muon_phi[im], MU_MASS); l1.pid=nr.Muon_pdgId[im];
      }
      else if (nr.Muon_pt[im]>l2.v.Pt()) {
	l4.v = l3.v; l4.pid=l3.pid;
	l3.v = l2.v; l3.pid=l2.pid;
	l2.v.SetPtEtaPhiM(nr.Muon_pt[im], nr.Muon_eta[im], nr.Muon_phi[im], MU_MASS); l2.pid=nr.Muon_pdgId[im];
      }
      else if (nr.Muon_pt[im]>l3.v.Pt()) {
	l4.v = l3.v; l4.pid=l3.pid;
	l3.v.SetPtEtaPhiM(nr.Muon_pt[im], nr.Muon_eta[im], nr.Muon_phi[im], MU_MASS); l3.pid=nr.Muon_pdgId[im];
      }
      else if (nr.Muon_pt[im]>l4.v.Pt()) {
	l4.v.SetPtEtaPhiM(nr.Muon_pt[im], nr.Muon_eta[im], nr.Muon_phi[im], MU_MASS); l4.pid=nr.Muon_pdgId[im];
      }

      nM++;
    }
    nE=0;
    for (uint ie=0; ie<nr.nElectron; ie++) {
      if (nr.Electron_pt[ie]<LEP_PT_CUT || abs(nr.Electron_eta[ie])>LEP_ETA_CUT) continue; //acceptance
      if (!nr.Electron_mvaSpring16HZZ_WPL[ie]) continue; //electron ID

      //electron isolation
      if (abs(nr.Electron_eta[ie]+nr.Electron_deltaEtaSC[ie])<1.479 && nr.Electron_miniPFRelIso_all[ie]>0.175) continue;
      else if ( nr.Electron_miniPFRelIso_all[ie]>0.159) continue;

      // d0 cut
      if (! ((abs(nr.Electron_eta[ie])<1.5  && abs(nr.Electron_dxy[ie])<0.0564) ||
	     (abs(nr.Electron_eta[ie])>=1.5 && abs(nr.Electron_dxy[ie])<0.222))) continue;

      if (nr.Electron_pt[ie]>l1.v.Pt()) {
	l4.v = l3.v; l4.pid=l3.pid;
	l3.v = l2.v; l3.pid=l2.pid;
	l2.v = l1.v; l2.pid=l1.pid;
	l1.v.SetPtEtaPhiM(nr.Electron_pt[ie], nr.Electron_eta[ie], nr.Electron_phi[ie], EL_MASS); l1.pid=nr.Electron_pdgId[ie];
      }
      else if (nr.Electron_pt[ie]>l2.v.Pt()) {
	l4.v = l3.v; l4.pid=l3.pid;
	l3.v = l2.v; l3.pid=l2.pid;
	l2.v.SetPtEtaPhiM(nr.Electron_pt[ie], nr.Electron_eta[ie], nr.Electron_phi[ie], EL_MASS); l2.pid=nr.Electron_pdgId[ie];
      }
      else if (nr.Electron_pt[ie]>l3.v.Pt()) {
	l4.v = l3.v; l4.pid=l3.pid;
	l3.v.SetPtEtaPhiM(nr.Electron_pt[ie], nr.Electron_eta[ie], nr.Electron_phi[ie], EL_MASS); l3.pid=nr.Electron_pdgId[ie];
      }
      else if (nr.Electron_pt[ie]>l4.v.Pt()) {
	l4.v.SetPtEtaPhiM(nr.Electron_pt[ie], nr.Electron_eta[ie], nr.Electron_phi[ie], EL_MASS); l4.pid=nr.Electron_pdgId[ie];
      }

      nE++;

    }

    if (l1.v.Pt()<30) continue;

    run=nr.run;
    ls=nr.luminosityBlock;
    evt=nr.event;

    genWeight=nr.genWeight;

    npu=nr.Pileup_nPU;
    npv=nr.PV_npvs;

    nJ15=nr.nJet;

    metPt=nr.MET_pt;
    metPhi=nr.MET_phi;

    puppiPt=nr.PuppiMET_pt;
    puppiPhi=nr.PuppiMET_phi;

    ot->Fill();

  }

  of->Write();
  of->Close();

}
