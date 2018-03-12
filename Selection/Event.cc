//#include "CaltechTriboson/Selection/WWZEventClass.hh"
#include "/afs/cern.ch/user/j/jlawhorn/CMSSW_9_3_2/src/CaltechTriboson/Selection/Event.hh"

WWZ::Event::Event() {
  Reset();
};

void WWZ::Event::Reset() {

  weight = 1.0;
  
  pileupWeight      = 1.0;
  pileupWeightUp    = 1.0;
  pileupWeightDown  = 1.0;
  triggerEffWeight  = 1.0;
  triggerEffSFWeight  = 1.0;

  lep1.pid = 0; lep1.dZ = -999; lep1.index=-1;
  lep1.pt = -999; lep1.eta = -999; lep1.phi = -999; lep1.m = -999; 
  lep1.isPrompt = false; lep1.passLooseMVA = false;
  lep1.genPt = -999; lep1.genEta = -999; lep1.genPhi = -999;

  lep2.pid = 0; lep2.dZ = -999; lep2.index=-1;
  lep2.pt = -999; lep2.eta = -999; lep2.phi = -999; lep2.m = -999; 
  lep2.isPrompt = false; lep2.passLooseMVA = false;
  lep2.genPt = -999; lep2.genEta = -999; lep2.genPhi = -999;

  lep3.pid = 0; lep3.dZ = -999; lep3.index=-1;
  lep3.pt = -999; lep3.eta = -999; lep3.phi = -999; lep3.m = -999; 
  lep3.isPrompt = false; lep3.passLooseMVA = false;
  lep3.genPt = -999; lep3.genEta = -999; lep3.genPhi = -999;

  lep4.pid = 0; lep4.dZ = -999; lep4.index=-1;
  lep4.pt = -999; lep4.eta = -999; lep4.phi = -999; lep4.m = -999; 
  lep4.isPrompt = false; lep4.passLooseMVA = false;
  lep4.genPt = -999; lep4.genEta = -999;  lep4.genPhi = -999;

  ZMass = -999;
  ZPt = -999;
  lep3MT = -999;
  lep4MT = -999;   
  phi0 = -999;
  theta0 = -999;
  phi = -999;
  theta1 = -999;
  theta2 = -999;   
  phiH = -999;
  MET = -999;
  MET_JESUp = -999;
  MET_JESDown = -999;
  NJet20 = 0;

  //struct Jet {
  //  float pt, eta, phi;
  //  float genPt, genEta, genPhi;
  //  float csv;
  //}

  jet1.pt=-999; jet1.eta=-999; jet1.phi=-999; jet1.e=-999;
  jet1.genPt=-999; jet1.genEta=-999; jet1.genPhi=-999;
  jet1.csv=-999;

  jet2.pt=-999; jet2.eta=-999; jet2.phi=-999; jet2.e=-999;
  jet2.genPt=-999; jet2.genEta=-999; jet2.genPhi=-999;
  jet2.csv=-999;

  jet3.pt=-999; jet3.eta=-999; jet3.phi=-999; jet3.e=-999;
  jet3.genPt=-999; jet3.genEta=-999; jet3.genPhi=-999;
  jet3.csv=-999;

  jet4.pt=-999; jet4.eta=-999; jet4.phi=-999; jet4.e=-999;
  jet4.genPt=-999; jet4.genEta=-999; jet4.genPhi=-999;
  jet4.csv=-999;

  NJet30 = 0;
  NBJet20 = 0;
  NBJet30 = 0;     
  minDRJetToLep3 = 9999;
  minDRJetToLep4 = 9999;
  GenZLepton1Pt = -999;
  GenZLepton1Eta = -999;
  GenZLepton2Pt = -999;
  GenZLepton2Eta = -999;
  GenZPt = -999;
  GenZEta = -999;
  GenWPlusLeptonPt = -999;
  GenWPlusLeptonEta = -999;
  GenWPlusPt = -999;
  GenWPlusEta = -999;
  GenWMinusLeptonPt = -999;
  GenWMinusLeptonEta = -999;
  GenWMinusPt = -999;
  GenWMinusEta = -999;
  GenMET = -999;
  pt_zeta = -999;
  pt_zeta_vis = -999.;
  
};

ClassImp(WWZ::Event);
