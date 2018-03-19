#ifndef CALTECHTRIBOSON_SELECTION_OBJECT_HH
#define CALTECHTRIBOSON_SELECTION_OBJECT_HH

#include "TLorentzVector.h"
#include "TObject.h"

namespace WWZ {

  class Obj : public TObject {
  public:
    Obj() {
      fisReco=true; fisGen=false;
      Reset();
    };
    
    Obj(bool isReco, bool isGen) {
      fisReco=isReco; fisGen=isGen;
      Reset();
    };
    
    ~Obj() {}
    
    void Reset() {
      recV.SetPtEtaPhiM(0,0,0,0); genV.SetPtEtaPhiM(0,0,0,0); 
      fisMatched=false;
    };
    
    TLorentzVector recV;
    TLorentzVector genV;
    bool fisGen, fisReco, fisMatched;
    
    ClassDef(Obj,1);
  };
  
}
#endif
