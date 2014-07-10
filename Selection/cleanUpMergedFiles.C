//-------------------------------------------------------------------
// Clean up merged files from lxbtch
//
// execute with:
// root -l -q cleanUpMergedFiles(_infile_, _outfile_)
//
// Code based on Jay Lawhorn 11/4/13
//-------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TMath.h>
#include <TChain.h>
#include <TH1.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Math/LorentzVector.h"

#include "modules/Delphes.h"                   // delphes
#include "ExRootAnalysis/ExRootTreeReader.h"   // delphes
#include "classes/DelphesClasses.h"            // delphes

#endif

void cleanUpMergedFiles(TString infilename="/afs/cern.ch/work/k/klawhorn/SnowmassSamples/PhaseII/Configuration4v2/Working/LL-4p-0-100-v1510_14TEV_temp.root",
			  TString outfilename="/afs/cern.ch/work/a/ariostas/SnowmassSamples/PhaseII/Configuration4v2/Working/LL-4p-0-100-v1510_14TEV_clean.root") {

  // set up input/output variables and file
  UInt_t nEvents, nLeptons, nJets;
  Double_t weight;
  
  Photon *photon1=0, *photon2=0;
  Jet *bJet1=0, *bJet2=0;

  TFile* infile = new TFile(infilename); assert(infile);
  TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("weight",		&weight);
  intree->SetBranchAddress("nJets",			&nJets);
  intree->SetBranchAddress("nLeptons",		&nLeptons);
  intree->SetBranchAddress("photon1",		&photon1);
  intree->SetBranchAddress("photon2",		&photon2);
  intree->SetBranchAddress("bJet1",          &bJet1);
  intree->SetBranchAddress("bJet2",          &bJet2);

  TTree* infotree = (TTree*) infile->Get("Info"); assert(infotree);
  infotree->SetBranchAddress("nEvents",      &nEvents);

  Long64_t totalEvents=0;

  for (UInt_t iEntry=0; iEntry<infotree->GetEntries(); iEntry++) {
    infotree->GetEntry(iEntry);
    totalEvents+=nEvents;
  }

  TFile *outFile = new TFile(outfilename, "RECREATE");

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");
	outTree->Branch("weight",	&weight,    "weight/D");
	outTree->Branch("nJets",	&nJets,    "nJets/i");
	outTree->Branch("nLeptons",    &nLeptons,    "nLeptons/i");  
	outTree->Branch("photon1",   	"Photon", &photon1);     
	outTree->Branch("photon2",    	"Photon", &photon2); 
	outTree->Branch("bJet1",   		"Jet", &bJet1); 
	outTree->Branch("bJet2",   		"Jet", &bJet2); 

  for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
    intree->GetEntry(iEntry);
	weight/=Double_t(totalEvents);

    outTree->Fill();

  }
  outFile->Write();
  outFile->Close();
  
	cout << totalEvents << endl;

}
