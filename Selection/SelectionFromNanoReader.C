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

#include "NanoReader.hh"

#endif

void SelectionFromNanoReader() {

  TChain *ch = new TChain("Events");
  ch->Add("DYJetsToLL_nano.root");

  NanoReader nr(ch);
  int nEntries = nr.GetEntries(); 

  for (int i=0; i<100; i++) {
    nr.LoadTree(i); nr.GetEntry(i);

    if ((nr.nElectron + nr.nMuon)<2) continue;
    
    std::cout << nr.nElectron << ", " << nr.nMuon << std::endl;

  }

}
