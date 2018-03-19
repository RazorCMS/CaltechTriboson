#include "CaltechTriboson/Utils/Event.hh"

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

  lep1.ResetLep();
  lep2.ResetLep();
  lep3.ResetLep();
  lep4.ResetLep();

  jet1.ResetJet();
  jet2.ResetJet();
  jet3.ResetJet();
  jet4.ResetJet();

  NJet30 = 0;
  NBJet20 = 0;
  NBJet30 = 0;     
  minDRJetToLep3 = 9999;

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
