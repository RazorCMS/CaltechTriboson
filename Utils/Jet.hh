#ifndef CALTECHTRIBOSON_SELECTION_JET_HH
#define CALTECHTRIBOSON_SELECTION_JET_HH

#include "CaltechTriboson/Utils/Obj.hh"

namespace WWZ {

  class Jet : public Obj {
  public:
    Jet() : Obj() {};
    Jet(float pt, float eta, float phi, float e, bool isReco=true, bool isGen=false)
      : Obj(isReco, isGen)
    {
      recV.SetPtEtaPhiE(pt, eta, phi, e);
      csv=-999; jec=-1.0; isBtag=false;
    }
    ~Jet() {}

    void ResetJet() {
      Reset();
      csv=-999; jec=-1.0; isBtag=false;
    }

    float csv, jec;
    bool isBtag;

    ClassDef(Jet,1);
  };

}
#endif
