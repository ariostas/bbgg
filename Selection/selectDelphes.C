//-------------------------------------------------------------------
// Select bbtautau events from delphes
//
// execute with:
// root -l -q selectDelphes.C+\(\"_inputfile_name_\",_file_cross_section_,\"_outputfile_name\"\)
//
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
#include "JEC/JECHelper/interface/jcorr.h"

#include "modules/Delphes.h"                   // delphes
#include "ExRootAnalysis/ExRootTreeReader.h"   // delphes
#include "classes/DelphesClasses.h"            // delphes

#endif

int puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t betastar);
double doJetcorr(mithep::jcorr *corrector,Jet* ijet,double rho_2,double rho_1,double rho_0);

void selectDelphes(const TString inputfile="root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/HHToGGBB_14TeV/HHToGGBB_14TeV_0.root",
const Double_t xsec=1341.36923) {
	
	mithep::jcorr *corrector = new mithep::jcorr;
	char *PATH = getenv("CMSSW_BASE"); assert(PATH);
	TString path(TString::Format("%s/src/JEC/JECHelper/data/JetCorrections_phase2/", PATH));
	corrector->AddJetCorr(path+"Delphes_V1_MC_L1FastJet_AK4PF.txt");
	corrector->AddJetCorr(path+"Delphes_V1_MC_L2Relative_AK4PF.txt");
	corrector->AddJetCorr(path+"Delphes_V1_MC_L3Absolute_AK4PFl1.txt");
	corrector->setup();

	const TString inputFile = inputfile + ".txt";
	ifstream ifs(inputFile);

	assert(ifs.is_open());
	
	TString filename;
	
	cout << "Reading " << inputfile << endl;

	while(ifs >> filename){

	// read input input file
	TChain chain("Delphes");
	chain.Add("root://eoscms.cern.ch/" + (inputfile == "HHToGGBB_14TeV" ? "/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/" + inputfile + "/" + filename : "/store/group/upgrade/delphes/ProdJun14/" + inputfile + "/" + filename));
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	// set up branches to read in from file
	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchRho = treeReader->UseBranch("Rho");
	TClonesArray *branchJet = treeReader->UseBranch("RawJet");

	if (!(branchJet)) {
		cout << "File broken" << endl;
		return;
	}

	// set up loop variables
	Photon *photon;
	Jet *jet;

	// set up storage variables
	Photon *photon1=0, *photon2=0;
	Jet *bJet1=0, *bJet2=0;
	
	Int_t nPhoton1, nPhoton2, nbJet1, nbJet2;
	Double_t weight=3000000*xsec;

	// set up output variables and file
	UInt_t nEvents, nJets, nLeptons;
	Double_t corr1=0, corr2=0;

	const TString outfile = "/afs/cern.ch/work/a/ariostas/private/bbgg_temp/" + inputfile + "/" + filename;

	TFile *outFile = new TFile(outfile, "RECREATE");

	// tree to hold the number of events in the file before selection
	TTree *sampTree = new TTree("Info", "Info");
	sampTree->Branch("nEvents",       &nEvents,        "nEvents/i");
	nEvents=numberOfEntries;
	sampTree->Fill();

	// tree to hold information about selected events
	TTree *outTree = new TTree("Events", "Events");
	outTree->Branch("weight",		&weight,    "weight/D");  // number of jets
	outTree->Branch("nJets",		&nJets,    "nJets/i");  // number of jets
	outTree->Branch("nLeptons",		&nLeptons,    "nLeptons/i");  // number of jets
	outTree->Branch("photon1",   	"Photon", &photon1);      // 4-vector for leading lepton
	outTree->Branch("photon2",    	"Photon", &photon2);      // 4-vector for second lepton
	outTree->Branch("bJet1",   		"Jet", &bJet1); // 4-vector for leading jet
	outTree->Branch("bJet2",   		"Jet", &bJet2); // 4-vector for second jet
	outTree->Branch("corr1",		&corr1,    "corr1/D");  // number of jets
	outTree->Branch("corr2",		&corr2,    "corr2/D");  // number of jets


	for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
		treeReader->ReadEntry(iEntry);

		// Reset index variables
		nPhoton1=nPhoton2=nbJet1=nbJet2=-1;
		
		// Store the number of leptons
		nLeptons = branchElectron->GetEntries() + branchMuon->GetEntries();
		nJets = corr1 = corr2 = 0;
		
		bool isoVar = true;
			
		// Select photons
		for (Int_t iPhoton=0; iPhoton<branchPhoton->GetEntries(); iPhoton++) { // Photon loop
			photon = (Photon*) branchPhoton->At(iPhoton);
		
			if(photon->IsolationVar > 0.5) isoVar = false;
			
			if((photon->PT > 25) && (fabs(photon->Eta) < 2.5)){
				
				if(nPhoton1 == -1){
					
					nPhoton1 = iPhoton;
					photon1 = (Photon*) branchPhoton->At(nPhoton1);
					
				}
				else if(photon->PT > photon1->PT){
										
					nPhoton2 = nPhoton1;
					photon2 = (Photon*) branchPhoton->At(nPhoton2);
				
					nPhoton1 = iPhoton;
					photon1 = (Photon*) branchPhoton->At(nPhoton1);
				
				}
				else if(nPhoton2 == -1){
				
					nPhoton2 = iPhoton;
					photon2 = (Photon*) branchPhoton->At(nPhoton2);
				
				}
				else if(photon->PT > photon2->PT){
				
					nPhoton2 = iPhoton;
					photon2 = (Photon*) branchPhoton->At(nPhoton2);
				
				}
				
			}
			
		}// End photon loop

		// Select bjets and count jets
		for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // bjet loop
			jet = (Jet*) branchJet->At(iJet);
			
			Double_t corr = doJetcorr(corrector,jet,((Rho*) branchRho->At(2))->Rho,((Rho*) branchRho->At(1))->Rho,((Rho*) branchRho->At(0))->Rho);
			
			if ((fabs(jet->Eta)<2.5) && (corr*jet->PT > 30)) nJets++;
		
			if((jet->BTag) && (corr*jet->PT > 30) && (fabs(jet->Eta) < 2.4) && (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar) == 0)){
				if(nbJet1 == -1){
					corr1 = corr;
					nbJet1 = iJet;
					bJet1 = (Jet*) branchJet->At(nbJet1);
				}
				else if(corr*jet->PT > corr1*bJet1->PT){
					
					corr2 = corr1;
					nbJet2 = nbJet1;
					bJet2 = (Jet*) branchJet->At(nbJet2);
					
					corr1 = corr;
					nbJet1 = iJet;
					bJet1 = (Jet*) branchJet->At(nbJet1);
				
				}
				else if(nbJet2 == -1){
					
					corr2 = corr;
					nbJet2 = iJet;
					bJet2 = (Jet*) branchJet->At(nbJet2);
				
				}
				else if(corr*jet->PT > corr2*bJet2->PT){
					
					corr2 = corr;
					nbJet2 = iJet;
					bJet2 = (Jet*) branchJet->At(nbJet2);
				
				}
				
			}
			
		}// End bjet loop
		
		if((nPhoton1!=-1) && (nPhoton2!=-1) && (nbJet1!=-1) && (nbJet2!=-1) && (isoVar)) outTree->Fill();

	} // end event loop

	// save file
	outFile->Write();
	// close file
	outFile->Close();

	cout << "----SUMMARY----" << endl;
	cout << " input file " << inputfile << " selection done " << endl;

	}

}

int puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t betastar) {
  
	Float_t MeanSqDeltaRMaxBarrel=0.07;
	Float_t BetaMinBarrel=0.87;
	Float_t MeanSqDeltaRMaxEndcap=0.07;
	Float_t BetaMinEndcap=0.85;

	//cout << eta << ", " << meanSqDeltaR << ", " << betastar << ": ";

	if (fabs(eta)<1.5) {
		if ((meanSqDeltaR<MeanSqDeltaRMaxBarrel)&&(betastar<BetaMinBarrel)) {
			//cout << "barrel 0" << endl;
			return 0;
		}
		else {
			//cout << "barrel 1" << endl;
			return 1;
		}
	}
	else if (fabs(eta)<4.0) {
		if ((meanSqDeltaR<MeanSqDeltaRMaxEndcap)&&(betastar<BetaMinEndcap)) {
			//cout << "endcap 0" << endl;
			return 0;
		}
		else {
			//cout << "endcap 1" << endl;
			return 1;
		}
	}
	//cout << "forward 1" << endl;
	return 1;

}

double doJetcorr(mithep::jcorr *corrector,Jet* ijet,double rho_2,double rho_1,double rho_0)
{
  double area = TMath::Sqrt(ijet->AreaX*ijet->AreaX+ijet->AreaY*ijet->AreaY);
  //std::cout << "Area pT " << area << std::endl;
  double rrho=0;
  if(fabs(ijet->Eta)>4.0) rrho=rho_2;
  else if(fabs(ijet->Eta)>2.5) rrho=rho_1;
  else rrho=rho_0;
  double corr = corrector->getCorrection(ijet->PT,ijet->Eta,rrho,area);
  //std::cout << ijet->PT << " " << corr << std::endl;
  //return 1.0;
  return corr;
}
