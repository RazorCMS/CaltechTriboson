{

  //gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  //gROOT->ProcessLine(".include $CMSSW_BASE/src");

  //TString path = gSystem->GetIncludePath();
  //path += "-I. -I$ROOTSYS/src -I$CMSSW_BASE/src";
  //gSystem->SetIncludePath(path.Data());

  gROOT->ProcessLine(".L $CMSSW_BASE/src/CaltechTriboson/Utils/Obj.cc+");
  gROOT->ProcessLine(".L $CMSSW_BASE/src/CaltechTriboson/Utils/Lep.cc+");
  gROOT->ProcessLine(".L $CMSSW_BASE/src/CaltechTriboson/Utils/Jet.cc+");
  gROOT->ProcessLine(".L $CMSSW_BASE/src/CaltechTriboson/Utils/Event.cc+");

  //gSystem->Load("$CMSSW_BASE/lib/slc6_amd64_gcc491/libBaconAnaDataFormats.so");
    
  //gROOT->Macro("$CMSSW_BASE/src/BaconAna/macros/setRootEnv.C+");  

}
