#ifndef CALTECHTRIBOSON_SELECTION_EVENTS_HH
#define CALTECHTRIBOSON_SELECTION_EVENTS_HH

#include "TLorentzVector.h"
#include "CaltechTriboson/Utils/Lep.hh"
#include "CaltechTriboson/Utils/Jet.hh"

namespace WWZ {

  class Event : public TObject {
  public:
    Event();
    ~Event() {};
    
    void Reset();

    float weight;
    float pileupWeight, pileupWeightUp, pileupWeightDown;
    float triggerEffWeight;
    float triggerEffSFWeight;

    Lep lep1, lep2, lep3, lep4;
    Jet jet1, jet2, jet3, jet4;

    float phi0, theta0, phi, theta1, theta2, phiH;
    
    int NPU;
    float ZMass, ZPt;
    float lep3MT, lep4MT;
    float lep34MT;
    
    int NJet20;
    int NJet30;
    int NBJet20;
    int NBJet30;
    float minDRJetToLep3;
    float minDRJetToLep4;
    
    float MET, MET_JESUp, MET_JESDown; 
    float METPhi;
    float METPuppiPt;
    float METPuppiPhi;
    float pt_zeta;
    float pt_zeta_vis;
    unsigned int run, lumi, event;
  
    float GenZLepton1Pt,GenZLepton1Eta;
    float GenZLepton2Pt,GenZLepton2Eta;
    float GenZPt,GenZEta;
    float GenWPlusLeptonPt,GenWPlusLeptonEta;
    float GenWPlusPt,GenWPlusEta;
    float GenWMinusLeptonPt,GenWMinusLeptonEta;
    float GenWMinusPt,GenWMinusEta;
    float GenMET;

    ClassDef(Event,1);

  };
}
#endif
