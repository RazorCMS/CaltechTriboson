float eleEffAreaMean(float scEta);
int dataTakingPeriod(int run);

// MUON ID and isolation
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
bool muLooseID(bool isLoose, float ip3dsig) {
  if (isLoose && fabs(ip3dsig) > 4) return true;
  return false;
};

bool muLooseIso(float muPt, float chIso, float neIso, float nhIso, float puIso) {
  if ((chIso + fmax(0.0, neIso + nhIso - 0.5 * puIso))/muPt < 0.2) return true;
  return false;
};

// ELECTRON ID and isolation
//
// Giovanni Zevi Della Porta
// tried to match working points of 80X electron MVA with 74X electron MVA
// pT < 10 bin should use "HZZ" MVA
// pT > 10 bin should use "GeneralPurpose" MVA
// https://indico.cern.ch/event/590228/contributions/2380031/attachments/1375541/2088587/EGMSUS_newIDs_17Nov16.pdf

bool eleVetoMvaID(float pt, float scEta, float mvaVal, float hzzMvaVal) {

  bool ans=false;
  if (pt<10) { //use hZZ MVA
    if      (fabs(scEta)<0.8   && hzzMvaVal>0.46)  ans=true;
    else if (fabs(scEta)<1.479 && hzzMvaVal>-0.03) ans=true;
    else if (                     hzzMvaVal>0.06)  ans=true;
  }
  else if (pt<15) { //use general MVA
    if      (fabs(scEta)<0.8   && mvaVal>-0.48) ans=true;
    else if (fabs(scEta)<1.479 && mvaVal>-0.67) ans=true;
    else if (                     mvaVal>-0.49) ans=true;
  }
  else if (pt<25) {
    if      (fabs(scEta)<0.8   && mvaVal>(-0.48-0.037*(pt-15))) ans=true;
    else if (fabs(scEta)<1.479 && mvaVal>(-0.67-0.024*(pt-15))) ans=true;
    else if (                     mvaVal>(-0.49-0.034*(pt-15))) ans=true;
  }
  else {
    if      (fabs(scEta)<0.8   && mvaVal>-0.85) ans=true;
    else if (fabs(scEta)<1.479 && mvaVal>-0.85) ans=true;
    else if (                     mvaVal>-0.85) ans=true;
  }

  return ans;
};

bool eleVetoIso(float pt, float scEta, float chIso, float neIso, float nhIso, float rho) {

  if (fabs(scEta)<1.479 && 
      (chIso + fmax(0.0, neIso + nhIso - eleEffAreaMean(scEta)*rho))/pt < 0.175) return true;

  else if ((chIso + fmax(0.0, neIso + nhIso - eleEffAreaMean(scEta)*rho))/pt < 0.159) return true;

  return false;
};

float eleEffAreaMean(float scEta) {
  float effArea=0.0;

  if (fabs(scEta) < 1.0) {
    effArea = 0.1752;
  } else if (fabs(scEta) < 1.479) {
    effArea = 0.1862;
  } else if (fabs(scEta) < 2.0) {
    effArea = 0.1411;
  } else if (fabs(scEta) < 2.2) {
    effArea = 0.1534;
  } else if (fabs(scEta) < 2.3) {
    effArea = 0.1903;
  } else if (fabs(scEta) < 2.4) {
    effArea = 0.2243;
  } else if (fabs(scEta) < 2.5) {
    effArea = 0.2687;
  }
  return effArea;
};

//Jet B-tagging
//From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X

bool jetLooseCISV(float jetCISV, int run) {
  if      (dataTakingPeriod(run)==2015||run==2015) {
    return jetCISV > 0.605;
  }
  else if (dataTakingPeriod(run)==2016||run==2016) {
    return jetCISV > 0.5426;
  }
  else if (dataTakingPeriod(run)==2017||run==2017) {
    return jetCISV > 0.5803;
  }

  return false;
}


float deltaR(float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1-eta2;
  float dphi = phi1-phi2;
  while (dphi>TMath::Pi()) dphi-=2*TMath::Pi();
  while (dphi<=-TMath::Pi()) dphi+=2*TMath::Pi();

  return sqrt(deta*deta+dphi*dphi);
}

int dataTakingPeriod(int run) {
  //from cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification
  if      (run>=254231 && run<=260627) return 2015;
  else if (run>=271036 && run<=284044) return 2016;
  else if (run>=294927 && run<=306462) return 2017;
  else return 0;
}
