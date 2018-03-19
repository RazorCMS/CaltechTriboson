#ifndef CALTECHTRIBOSON_SELECTION_LEPTON_HH
#define CALTECHTRIBOSON_SELECTION_LEPTON_HH

#include "CaltechTriboson/Utils/Obj.hh"

namespace WWZ {

  class Lep: public Obj {
  public:

    Lep() : Obj() {};

    Lep(float pt, float eta, float phi, float m, float dZ, int pid, bool isReco=true, bool isGen=false) 
      : Obj(isReco, isGen)
    {
      recV.SetPtEtaPhiM(pt, eta, phi, m);
      fdZ=dZ; fpid=pid;
      isPrompt=false; passLoose=false; passTight=false;
    };
    ~Lep() {}

    void ResetLep() {
      Reset();
      fdZ=-999; fpid=0;
      isPrompt=false; passLoose=false; passTight=false;
    }

    float fdZ;
    int fpid, findex;
    bool isPrompt, passLoose, passTight;

    ClassDef(Lep,1);
  };

}
#endif
