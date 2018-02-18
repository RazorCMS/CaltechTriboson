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

#endif

void Selection(TString infile="/eos/cms/store/group/phys_susy/razor/run2/Run2RazorNtupleV3.15/MC_Summer16/RunIISpring16/v5/sixie/WWZJetsTo4L2Nu_4f_TuneCUETP8M1_13TeV_aMCatNLOFxFx_pythia8/Run2RazorNtuplerV3p15_ToCERN_MC_Summer16_25ns_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2_v5_v1/170815_180614/0000/razorNtuple_10.root",
		  TString outfile="test.root"
		  ) {

  const float EL_MASS = 0.000511;
  const float MU_MASS = 0.105658;
  const float Z_MASS = 91.1876;

  const float LEP_ETA_CUT = 2.4;
  const float LEP_PT_CUT = 10;

  TFile *f = new TFile(infile,"read");

  TTree *t = (TTree*) f->Get("ntuples/RazorEvents");

  Bool_t          isData;
  UInt_t          runNum;
  UInt_t          lumiNum;
  UInt_t          eventNum;
  Int_t           nMuons;
  Float_t         muonE[700];   //[nMuons]
  Float_t         muonPt[700];   //[nMuons]
  Float_t         muonEta[700];   //[nMuons]
  Float_t         muonPhi[700];   //[nMuons]
  Int_t           muonCharge[700];   //[nMuons]
  Int_t           nElectrons;
  Float_t         eleE[700];   //[nElectrons]
  Float_t         elePt[700];   //[nElectrons]
  Float_t         eleEta[700];   //[nElectrons]
  Float_t         elePhi[700];   //[nElectrons]
  Float_t         eleCharge[700];   //[nElectrons]
  Int_t           nJets;
  Float_t         jetE[900];   //[nJets]
  Float_t         jetPt[900];   //[nJets]
  Float_t         jetEta[900];   //[nJets]
  Float_t         jetPhi[900];   //[nJets]
  Float_t         jetCSV[900];   //[nJets]
  Float_t         jetCISV[900];   //[nJets]
  Float_t         metPt;
  Float_t         metPhi;

  Int_t           nGenParticle;
  Int_t           gParticleMotherId[4000];   //[nGenParticle]
  Int_t           gParticleMotherIndex[4000];   //[nGenParticle]
  Int_t           gParticleId[4000];   //[nGenParticle]
  Int_t           gParticleStatus[4000];   //[nGenParticle]
  Float_t         gParticleE[4000];   //[nGenParticle]
  Float_t         gParticlePt[4000];   //[nGenParticle]
  Float_t         gParticleEta[4000];   //[nGenParticle]
  Float_t         gParticlePhi[4000]; //[nGenParticle]

  t->SetBranchAddress("isData",	         &isData);
  t->SetBranchAddress("runNum",	       	 &runNum);
  t->SetBranchAddress("lumiNum",	 &lumiNum);
  t->SetBranchAddress("eventNum",	 &eventNum);
  t->SetBranchAddress("nMuons",	       	 &nMuons);
  t->SetBranchAddress("muonE",      	 &muonE);
  t->SetBranchAddress("muonPt",     	 &muonPt);
  t->SetBranchAddress("muonEta",    	 &muonEta);
  t->SetBranchAddress("muonPhi",    	 &muonPhi);
  t->SetBranchAddress("muonCharge",    	 &muonCharge);
  t->SetBranchAddress("nElectrons",   	 &nElectrons);
  t->SetBranchAddress("eleE",       	 &eleE);
  t->SetBranchAddress("elePt",      	 &elePt);
  t->SetBranchAddress("eleEta",     	 &eleEta);
  t->SetBranchAddress("elePhi",     	 &elePhi);
  t->SetBranchAddress("eleCharge",  	 &eleCharge);
  t->SetBranchAddress("nJets",	       	 &nJets);
  t->SetBranchAddress("jetE",       	 &jetE);
  t->SetBranchAddress("jetPt",      	 &jetPt);
  t->SetBranchAddress("jetEta",     	 &jetEta);
  t->SetBranchAddress("jetPhi",     	 &jetPhi);
  t->SetBranchAddress("jetCSV",     	 &jetCSV);
  t->SetBranchAddress("jetCISV",   	 &jetCISV);
  t->SetBranchAddress("metPt",	       	 &metPt);
  t->SetBranchAddress("metPhi",          &metPhi);

  //t->SetBranchAddress("nGenParticle",         &nGenParticle);
  //t->SetBranchAddress("gParticleMotherId",    &gParticleMotherId);
  //t->SetBranchAddress("gParticleMotherIndex", &gParticleMotherIndex);
  //t->SetBranchAddress("gParticleId",          &gParticleId);
  //t->SetBranchAddress("gParticleStatus",      &gParticleStatus);
  //t->SetBranchAddress("gParticleE",           &gParticleE);
  //t->SetBranchAddress("gParticlePt",          &gParticlePt);
  //t->SetBranchAddress("gParticleEta",         &gParticleEta);
  //t->SetBranchAddress("gParticlePhi",         &gParticlePhi);

  TFile *of = new TFile(outfile,"recreate");
  TTree *ot = new TTree("wwz","wwz");

  //
  //int wp1_f=0, wp2_f=0, wm1_f=0, wm2_f=0, z1_f=0, z2_f=0;
  //TLorentzVector z(0,0,0,0),   wp(0,0,0,0), wm(0,0,0,0);
  //TLorentzVector z1(0,0,0,0),  z2(0,0,0,0);
  //TLorentzVector wp1(0,0,0,0), wp2(0,0,0,0);
  //TLorentzVector wm1(0,0,0,0), wm2(0,0,0,0);
  //
  //float genZ_m=0, genWp_m=0, genWm_m=0;
  //float genWW_m=0, genWWZ_m=0;
  //float genLL_m=0, genL1Z_m=0, genL2Z_m=0, genLLZ_m=0;

  //ot->Branch("genZ",  "TLorentzVector", &z);
  //ot->Branch("genWp", "TLorentzVector", &wp);
  //ot->Branch("genWm", "TLorentzVector", &wm);
  //ot->Branch("genZ1",  "TLorentzVector", &z1);
  //ot->Branch("genZ2",  "TLorentzVector", &z2);
  //ot->Branch("genWp1", "TLorentzVector", &wp1);
  //ot->Branch("genWp2", "TLorentzVector", &wp2);
  //ot->Branch("genWm1", "TLorentzVector", &wm1);
  //ot->Branch("genWm2", "TLorentzVector", &wm2);

  //ot->Branch("z1_f",  &z1_f,  "z1_f/I");
  //ot->Branch("z2_f",  &z2_f,  "z2_f/I");
  //ot->Branch("wp1_f", &wp1_f, "wp1_f/I");
  //ot->Branch("wp2_f", &wp2_f, "wp2_f/I");
  //ot->Branch("wm1_f", &wm1_f, "wm1_f/I");
  //ot->Branch("wm2_f", &wm2_f, "wm2_f/I");
  //
  //ot->Branch("genZ_m",   &genZ_m,   "genZ_m/F");
  //ot->Branch("genWp_m",  &genWp_m,  "genWp_m/F");
  //ot->Branch("genWm_m",  &genWm_m,  "genWm_m/F");
  //ot->Branch("genWW_m",  &genWW_m,  "genWW_m/F");
  //ot->Branch("genWWZ_m", &genWWZ_m, "genWWZ_m/F");
  //ot->Branch("genLL_m",  &genLL_m,  "genLL_m/F");
  //ot->Branch("genL1Z_m", &genL1Z_m, "genL1Z_m/F");
  //ot->Branch("genL2Z_m", &genL2Z_m, "genL2Z_m/F");
  //ot->Branch("genLLZ_m", &genLLZ_m, "genLLZ_m/F");

  int bestZ1=-1;
  int bestZ2=-1;

  float L1_pt=0, L1_eta=0, L1_phi=0, L1_m=0; int L1_pid=0;
  float L2_pt=0, L2_eta=0, L2_phi=0, L2_m=0; int L2_pid=0;
  float L3_pt=0, L3_eta=0, L3_phi=0, L3_m=0; int L3_pid=0;
  float L4_pt=0, L4_eta=0, L4_phi=0, L4_m=0; int L4_pid=0;

  float MET_pt=0, MET_phi=0;

  float Z_m=0, WW_mvis=0, WWZ_mvis=0;

  ot->Branch("L1_pt",  &L1_pt,  "L1_pt/F"); //highest pT lepton from Z
  ot->Branch("L1_eta", &L1_eta, "L1_eta/F");
  ot->Branch("L1_phi", &L1_phi, "L1_phi/F");
  ot->Branch("L1_m",   &L1_m,   "L1_m/F");
  ot->Branch("L1_pid", &L1_pid, "L1_pid/I");

  ot->Branch("L2_pt",  &L2_pt,  "L2_pt/F"); // 2nd lepton from Z
  ot->Branch("L2_eta", &L2_eta, "L2_eta/F");
  ot->Branch("L2_phi", &L2_phi, "L2_phi/F");
  ot->Branch("L2_m",   &L2_m,   "L2_m/F");
  ot->Branch("L2_pid", &L2_pid, "L2_pid/I");

  ot->Branch("L3_pt",  &L3_pt,  "L3_pt/F"); //highest pT lepton from W
  ot->Branch("L3_eta", &L3_eta, "L3_eta/F");
  ot->Branch("L3_phi", &L3_phi, "L3_phi/F");
  ot->Branch("L3_m",   &L1_m,   "L3_m/F");
  ot->Branch("L3_pid", &L3_pid, "L3_pid/I");

  ot->Branch("L4_pt",  &L4_pt,  "L4_pt/F"); // 2nd lepton from W
  ot->Branch("L4_eta", &L4_eta, "L4_eta/F");
  ot->Branch("L4_phi", &L4_phi, "L4_phi/F");
  ot->Branch("L4_m",   &L4_m,   "L4_m/F");
  ot->Branch("L4_pid", &L4_pid, "L4_pid/I");

  ot->Branch("MET_pt",  &MET_pt,  "MET_pt/F");
  ot->Branch("MET_phi", &MET_phi, "MET_phi/F");

  ot->Branch("Z_m", &Z_m, "Z_m/F");
  ot->Branch("WW_mvis",  &WW_mvis,  "WW_mvis/F");
  ot->Branch("WWZ_mvis", &WWZ_mvis, "WWZ_mvis/F");

  for (uint i=0; i<t->GetEntries(); i++) {
  //for (uint i=0; i<50; i++) {
    t->GetEntry(i);

    std::vector<int> eleIndex;
    std::vector<int> muIndex;
    for (int j=0; j<nElectrons; j++) {
      if (abs(eleEta[j])<LEP_ETA_CUT && elePt[j]>LEP_PT_CUT) {
	eleIndex.push_back(j);
      }
    }

    for (int j=0; j<nMuons; j++) {
      if (abs(muonEta[j])<LEP_ETA_CUT && muonPt[j]>LEP_PT_CUT) {
	muIndex.push_back(j);
      }
    }

    if(muIndex.size()+eleIndex.size()<4) continue;

    float bestZMassDelta=99999;
    float bestZMass=0;
    bestZ1=-1;
    bestZ2=-1;

    for (uint j=0; j<eleIndex.size(); j++) {
      for (uint k=j+1; k<eleIndex.size(); k++) {
	if (eleCharge[eleIndex[j]]==eleCharge[eleIndex[k]]) continue;

	TLorentzVector e1(0,0,0,0); e1.SetPtEtaPhiM(elePt[eleIndex[j]], eleEta[eleIndex[j]], elePhi[eleIndex[j]], EL_MASS);
	TLorentzVector e2(0,0,0,0); e2.SetPtEtaPhiM(elePt[eleIndex[k]], eleEta[eleIndex[k]], elePhi[eleIndex[k]], EL_MASS);

	if ( abs((e1+e2).M()-Z_MASS)< bestZMassDelta ) {
	  bestZMassDelta=abs((e1+e2).M()-Z_MASS);
	  bestZMass=(e1+e2).M();
	  bestZ1=eleIndex[j];
	  bestZ2=eleIndex[k];
	}	
      }
    }

    for (uint j=0; j<muIndex.size(); j++) {
      for (uint k=j+1; k<muIndex.size(); k++) {
	if (muonCharge[muIndex[j]]==muonCharge[muIndex[k]]) continue;

	TLorentzVector m1(0,0,0,0); m1.SetPtEtaPhiM(muonPt[muIndex[j]], muonEta[muIndex[j]], muonPhi[muIndex[j]], MU_MASS);
	TLorentzVector m2(0,0,0,0); m2.SetPtEtaPhiM(muonPt[muIndex[k]], muonEta[muIndex[k]], muonPhi[muIndex[k]], MU_MASS);

	if ( abs((m1+m2).M()-Z_MASS)< bestZMassDelta ) {
	  bestZMassDelta=abs((m1+m2).M()-Z_MASS);
	  bestZMass=(m1+m2).M();
	  bestZ1=muIndex[j]+100;
	  bestZ2=muIndex[k]+100;
	}
      }
    }

    if (bestZMassDelta==99999) continue;

    bool zee=false;    
    if (bestZ1>-1 && bestZ1<100) {
      zee=true;
      if (elePt[bestZ1]>elePt[bestZ2]) {
	L1_pt =elePt[bestZ1];
	L1_eta=eleEta[bestZ1];
	L1_phi=elePhi[bestZ1];
	L1_m  =EL_MASS;
	L1_pid=-eleCharge[bestZ1]*11;

	L2_pt =elePt[bestZ2];
	L2_eta=eleEta[bestZ2];
	L2_phi=elePhi[bestZ2];
	L2_m  =EL_MASS;
	L2_pid=-eleCharge[bestZ2]*11;
      } else {
	L1_pt =elePt[bestZ2];
	L1_eta=eleEta[bestZ2];
	L1_phi=elePhi[bestZ2];
	L1_m  =EL_MASS;
	L1_pid=-eleCharge[bestZ2]*11;

	L2_pt =elePt[bestZ1];
	L2_eta=eleEta[bestZ1];
	L2_phi=elePhi[bestZ1];
	L2_m  =EL_MASS;
	L2_pid=-eleCharge[bestZ1]*11;
      } 
    } else {
      if (muonPt[bestZ1-100]>muonPt[bestZ2-100]) {
	L1_pt =muonPt[bestZ1-100];
	L1_eta=muonEta[bestZ1-100];
	L1_phi=muonPhi[bestZ1-100];
	L1_m  =MU_MASS;
	L1_pid=-muonCharge[bestZ1-100]*13;

	L2_pt =muonPt[bestZ2-100];
	L2_eta=muonEta[bestZ2-100];
	L2_phi=muonPhi[bestZ2-100];
	L2_m  =MU_MASS;
	L2_pid=-muonCharge[bestZ2-100]*13;
      } else {
	L1_pt =muonPt[bestZ2-100];
	L1_eta=muonEta[bestZ2-100];
	L1_phi=muonPhi[bestZ2-100];
	L1_m  =MU_MASS;
	L1_pid=-muonCharge[bestZ2-100]*13;

	L2_pt =muonPt[bestZ1-100];
	L2_eta=muonEta[bestZ1-100];
	L2_phi=muonPhi[bestZ1-100];
	L2_m  =MU_MASS;
	L2_pid=-muonCharge[bestZ1-100]*13;
      }
    }

    vector<int> eleIndex2;
    vector<int> muIndex2;

    if (bestZ1>-1 && bestZ1<100) {
      for (uint j=0; j<eleIndex.size(); j++) {
	if (eleIndex[j]==bestZ1||eleIndex[j]==bestZ2) continue;
	eleIndex2.push_back(eleIndex[j]);
      }
      muIndex2=muIndex;
    } else { 
      for (uint j=0; j<muIndex.size(); j++) {
	if (muIndex[j]+100==bestZ1||muIndex[j]+100==bestZ2) continue;
	muIndex2.push_back(muIndex[j]);
      }
      eleIndex2=eleIndex;
    }

    if (eleIndex.size()!=eleIndex2.size() && muIndex.size()!=muIndex2.size()) {
      std::cout << "both lepton collections changed when removing Z's!" << std::endl;
      continue;
    }

    int bestPlep=-1, bestMlep=-1;

    for (uint j=0; j<eleIndex2.size(); j++) {
      if (eleCharge[eleIndex2[j]]==1 && (bestPlep==-1 || elePt[eleIndex2[j]]>L3_pt)) {
	bestPlep=eleIndex2[j];
	L3_pt =elePt [bestPlep];
	L3_eta=eleEta[bestPlep];
	L3_phi=elePhi[bestPlep];
	L3_m  =EL_MASS;
	L3_pid=-eleCharge[bestPlep]*11;
      } else if (eleCharge[eleIndex2[j]]==-1 && (bestMlep==-1 || elePt[eleIndex2[j]]>L4_pt)) {
	bestMlep=eleIndex2[j];
	L4_pt =elePt [bestMlep];
	L4_eta=eleEta[bestMlep];
	L4_phi=elePhi[bestMlep];
	L4_m  =EL_MASS;
	L4_pid=-eleCharge[bestMlep]*11;
      }
    }

    for (uint j=0; j<muIndex2.size(); j++) {
      if (muonCharge[muIndex2[j]]==1 && (bestPlep==-1 ||muonPt[muIndex2[j]]>L3_pt)) {
	bestPlep=muIndex2[j];
	L3_pt =muonPt [bestPlep];
	L3_eta=muonEta[bestPlep];
	L3_phi=muonPhi[bestPlep];
	L3_m  =MU_MASS;
	L3_pid=-muonCharge[bestPlep]*13;
      } else if (muonCharge[muIndex2[j]]==-1 && (bestMlep==-1 || muonPt[muIndex2[j]]>L4_pt)) {
	bestMlep=muIndex2[j];
	L4_pt =muonPt [bestMlep];
	L4_eta=muonEta[bestMlep];
	L4_phi=muonPhi[bestMlep];
	L4_m  =MU_MASS;
	L4_pid=-muonCharge[bestMlep]*13;
      }
    }

    if (bestPlep==-1 || bestMlep==-1) continue;

    if (L4_pt>L3_pt) {
      float tmpPt, tmpEta, tmpPhi, tmpM; int tmpPid;
      tmpPt=L3_pt; tmpEta=L3_eta; tmpPhi=L3_phi; tmpM=L3_m; tmpPid=L3_pid;
      L3_pt=L4_pt; L3_eta=L4_eta; L3_phi=L4_phi; L3_m=L4_m; L3_pid=L4_pid;
      L4_pt=tmpPt; L4_eta=tmpEta; L4_phi=tmpPhi; L4_m=tmpM; L4_pid=tmpPid;
    }

    TLorentzVector lep1(0,0,0,0); lep1.SetPtEtaPhiM(L1_pt, L1_eta, L1_phi, L1_m);
    TLorentzVector lep2(0,0,0,0); lep2.SetPtEtaPhiM(L2_pt, L2_eta, L2_phi, L2_m);
    TLorentzVector lep3(0,0,0,0); lep3.SetPtEtaPhiM(L3_pt, L3_eta, L3_phi, L3_m);
    TLorentzVector lep4(0,0,0,0); lep4.SetPtEtaPhiM(L4_pt, L4_eta, L4_phi, L4_m);

    Z_m = (lep1+lep2).M();
    WW_mvis = (lep3+lep4).M();
    WWZ_mvis = (lep1+lep2+lep3+lep4).M();

    MET_pt=metPt;
    MET_phi=metPhi;
    
    ot->Fill();

    //std::cout << "best Z mass: " << bestZMass << ", " << bestZ1 << ", " << bestZ2 << std::endl;

    //gen-level info
    /*    wp1_f=0; wp2_f=0; wm1_f=0; wm2_f=0; z1_f=0; z2_f=0;

    z .SetPtEtaPhiE(0,0,0,0);
    wp.SetPtEtaPhiE(0,0,0,0);
    wm.SetPtEtaPhiE(0,0,0,0);
    z1 .SetPtEtaPhiE(0,0,0,0);
    wp1.SetPtEtaPhiE(0,0,0,0);
    wm1.SetPtEtaPhiE(0,0,0,0);
    z2 .SetPtEtaPhiE(0,0,0,0);
    wp2.SetPtEtaPhiE(0,0,0,0);
    wm2.SetPtEtaPhiE(0,0,0,0);

    for (int j=0; j<nGenParticle; j++) {
      if (gParticleStatus[j]==22 && gParticleId[j]==23) {
	z.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
      }
      else if (gParticleStatus[j]==22 && gParticleId[j]==24) {
	wp.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
      }
      else if (gParticleStatus[j]==22 && gParticleId[j]==-24) {
	wm.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
      }
      else if ((abs(gParticleId[j])>10 && abs(gParticleId[j])<17) && gParticleStatus[j]==1) { 
	if (gParticleMotherId[j]==23 && gParticleId[j]>0) {
	  z1.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
	  z1_f = gParticleId[j];
	}
	else if (gParticleMotherId[j]==23) {
	  z2.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
	  z2_f = gParticleId[j];
	}
	else if (gParticleMotherId[j]==24 && gParticleId[j]>0) {
	  wp1.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
	  wp1_f = gParticleId[j];
	}
	else if (gParticleMotherId[j]==24) {
	  wp2.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
	  wp2_f = gParticleId[j];
	}
	else if (gParticleMotherId[j]==-24 && gParticleId[j]>0) {
	  wm1.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
	  wm1_f = gParticleId[j];
	}
	else if (gParticleMotherId[j]==-24) {
	  wm2.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
	  wm2_f = gParticleId[j];
	}
      } 
    }

    if (z.M()==0 && z1_f==0 && z2_f==0) {
      int c=0;
      for (int j=0; j<nGenParticle; j++) {
	if (abs(gParticleId[j])<11 || abs(gParticleId[j])>16) continue;
	if (abs(gParticleMotherId[j])<6 && gParticleMotherIndex[j]>0 && gParticleStatus[j]==1) {
	  if (c>1) {
	    std::cout << "too many compatible leptons for Z/gamma! " << i <<  std::endl;
	    //std::cout << gParticleId[j] << " from " << gParticleMotherId[j] << std::endl;
	  }

	  c++;
	  if (gParticleId[j]>0) {
	    z1.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
	    z1_f = gParticleId[j];
	  }
	  else {
	    z2.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
	    z2_f = gParticleId[j];
	  }
	}
      }
      z=z1+z2;
    }

    if (abs(z1.Eta())>LEP_ETA_CUT ||abs(z2.Eta())>LEP_ETA_CUT ||abs(wp2.Eta())>LEP_ETA_CUT ||abs(wm1.Eta())>LEP_ETA_CUT) continue;
    if (z1.Pt()<LEP_PT_CUT||z2.Pt()<LEP_PT_CUT||wp2.Pt()<LEP_PT_CUT||wm1.Pt()<LEP_PT_CUT) continue;


    genZ_m  = z.M();
    genWp_m = wp.M();
    genWm_m = wm.M();

    genWW_m  = (wp+wm).M();
    genWWZ_m = (wp+wm+z).M();
    genLL_m  = (wp2+wm1).M();
    genLLZ_m = (wp2+wm1+z).M();
    genL1Z_m = (wm1+z).M();
    genL2Z_m = (wp2+z).M();
    */    
    //ot->Fill();
  }
  
  of->Write();
  of->Close();

}
