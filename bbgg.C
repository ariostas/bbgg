/*
 * bbgg.C
 * 
 * This macro analyses bbgg signal and background samples.
 * More specifically, it analyses the ttH and ZH backgrounds.
 * To run use "root bbgg.C"
 * 
 * The event selection consists in requiring a leading photon with pt > 40,
 * and a subleading photon with pt > 25, both with |eta| < 2.5.
 * It also requires two btagged jets with pt > 30 and |eta| < 2.4.
 * Finally it requires no reconstructed leptons and less than four jets
 * with pt > 30 and |eta| < 2.5
 * 
 * The fit mass window consists in requiring a reconstructed mass
 * 120 < M_gg < 130 and 70 < M_bb < 200
 * 
 * The kinematic cut requires deltaR_gg < 2, deltaR_bb < 2 and that
 * min(deltaR_bg) > 1.5
 * 
 * The mass cut requires 120 < M_gg < 130 and 105 < M_bb < 145
 * 
 * Coded by: Andres Rios
 */


// include statements for all needed dependencies
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>

// include statements for Delphes
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

// Declare variables
const bool Signal=true, Background=false;

// Declare functions
void histogram(TH1D*, TH1D*, TCanvas*, const char*, const char*, const char*);
void histogram(TH1D*, TCanvas*, const char*, const char*, const char*);
Float_t deltaR( const Float_t, const Float_t, const Float_t, const Float_t);
void analyze(TString, Double_t, bool);
void saveResults();


// Initialize histograms



// Set up variables for event yields
Double_t selectionSignal=0, kinCutSignal=0, massCutSignal=0;
Double_t selectionBackground=0, kinCutBackground=0, massCutBackground=0;

/*
 * MAIN FUNCTION
 */

void bbgg(){
	
	// Analyse signal
	//analyze("HHToGGBB_14TeV", 0.089, Signal);
	
	// Analyze ttH background
	analyze("ttB-4p-0-900-v1510_14TEV", 2.6673*1000, Background);
	analyze("ttB-4p-900-1600-v1510_14TEV", 0.250469*1000, Background);
	analyze("ttB-4p-1600-2500-v1510_14TEV", 0.0237441*1000, Background);
	analyze("ttB-4p-2500-100000-v1510_14TEV", 0.00208816*1000, Background);
	
	// Analyze ZH background
	analyze("BB-4p-0-300-v1510_14TEV", 249.97710*1000, Background);
	analyze("BB-4p-300-700-v1510_14TEV", 35.23062*1000, Background);
	analyze("BB-4p-700-1300-v1510_14TEV ", 4.13743*1000, Background);
	analyze("BB-4p-1300-2100-v1510_14TEV", 0.41702*1000, Background);
	analyze("BB-4p-2100-100000_14TEV", 0.04770*1000, Background);
	
	// Save results
	saveResults();
	
}


/*
 * ANALYSIS
 */

void analyze(TString inputfile, Double_t crossSection, bool SorB)
{	
	
	// Set up temporal vatiables
	Double_t tempSelection=0, tempKinCut=0, tempMassCut=0;
	
	const TString inputFilet = inputfile + ".txt";
	ifstream ifs(inputFilet);

	assert(ifs.is_open());

	TString filename;
	
	cout << "Reading from folder " << inputfile << endl;

	TChain chain("Delphes");

	while(ifs >> filename){
		
		TString filenamef;
		if(SorB) filenamef = "root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/" + inputfile + "/" + filename;
		else filenamef = "root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/" + inputfile + "/" + filename;
		
		cout << "Reading " << filenamef << endl;
		chain.Add(filenamef);
	}
		
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();
	
	// Set up loop variables
	GenParticle *particle;
	Photon *photon;
	Jet *jet;
	
	// Set up storage variables
	Int_t nPhoton1, nPhoton2, nbJet1, nbJet2;
	Photon *photon1, *photon2;
	Jet *bJet1, *bJet2;
	LHEFEvent *event;
	UInt_t nJets, nLeptons; 
	Double_t weight, dRPhotons, dRBJets, dRPhotonBJet;
	TLorentzVector vPhoton1, vPhoton2, vBJet1, vBJet2, photonSystem, bJetSystem;
	
	
	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchEvent;
	(!SorB ? branchEvent = treeReader->UseBranch("Event"): branchEvent = 0);
			
	for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // Event loop
		treeReader->ReadEntry(iEntry);		
			
		// If it's a background check whether the event contains a Higgs, top, or Z
		bool H=false, t=false, Z=false;
		if(!SorB){
		
			for(Int_t iParticle=0;iParticle<branchParticle->GetEntries();iParticle++){
				particle = (GenParticle*) branchParticle->At(iParticle);
			
				if(std::abs(particle->PID)==25)H=true;
				if(std::abs(particle->PID)==6)t=true;
				if(std::abs(particle->PID)==23)Z=true;
			
			}
		}
		
		// Check if the event is Signal, ttH or ZH
		if(SorB || (H && t) || (H && Z)){
			
			// Reset index variables
			nPhoton1=nPhoton2=nbJet1=nbJet2=-1;
			
			// Store the number of leptons
			nLeptons = branchElectron->GetEntries() + branchMuon->GetEntries();
			nJets = 0;
				
			// Select photons
			for (Int_t iPhoton=0; iPhoton<branchPhoton->GetEntries(); iPhoton++) { // Photon loop
				photon = (Photon*) branchPhoton->At(iPhoton);
			
				if((photon->PT > 25) && (fabs(photon->Eta)) < 2.5){
					
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
					else if(jet->PT > bJet2->PT){
					
						nPhoton2 = iPhoton;
						photon2 = (Photon*) branchPhoton->At(nPhoton2);
					
					}
					
				}
				
			}// End photon loop
	
			// Select bjets and count jets
			for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // bjet loop
				jet = (Jet*) branchJet->At(iJet);
				
				if ((fabs(jet->Eta)<2.5) && (jet->PT > 30)) nJets++;
			
				if((jet->BTag) && (jet->PT > 30) && (fabs(jet->Eta) < 2.4)){
					if(nbJet1 == -1){
						nbJet1 = iJet;
						bJet1 = (Jet*) branchJet->At(nbJet1);
					}
					else if(jet->PT > bJet1->PT){
					
						nbJet2 = nbJet1;
						bJet2 = (Jet*) branchJet->At(nbJet2);
					
						nbJet1 = iJet;
						bJet1 = (Jet*) branchJet->At(nbJet1);
					
					}
					else if(nbJet2 == -1){
					
						nbJet2 = iJet;
						bJet2 = (Jet*) branchJet->At(nbJet2);
					
					}
					else if(jet->PT > bJet2->PT){
					
						nbJet2 = iJet;
						bJet2 = (Jet*) branchJet->At(nbJet2);
					
					}
					
				}
				
			}// End bjet loop
	
			// Check if there are two selected photons and bjets
			if((nPhoton1!=-1) && (nPhoton2!=-1) && (nbJet1!=-1) && (nbJet2!=-1)){
				
				if((photon1->PT > 40) && (nLeptons == 0) && (nJets < 4)){
					
					// Set up four-vectors
					vPhoton1.SetPtEtaPhiE(photon1->PT, photon1->Eta, photon1->Phi, photon1->E);
					vPhoton2.SetPtEtaPhiE(photon2->PT, photon2->Eta, photon2->Phi, photon2->E);
					vBJet1.SetPtEtaPhiM(bJet1->PT, bJet1->Eta, bJet1->Phi, bJet1->Mass);
					vBJet2.SetPtEtaPhiM(bJet2->PT, bJet2->Eta, bJet2->Phi, bJet2->Mass);
				
					photonSystem = vPhoton1 + vPhoton2;
					bJetSystem = vBJet1 + vBJet2;
					
					// Check for fit mass window requirements
					if((photonSystem.M() > 120) && (photonSystem.M() < 130) && (bJetSystem.M() > 70) && (bJetSystem.M() < 200)){
						
						weight = 1;
						if(!SorB){
							
							event = (LHEFEvent*) branchEvent->At(0);
							weight = event->Weight;
							
						}
						
						tempSelection += weight;
						
						// Calculate all required deltaRs
						dRPhotons = deltaR(photon1->Eta, photon2->Eta, photon1->Phi, photon2->Phi);
						dRBJets = deltaR(bJet1->Eta, bJet2->Eta, bJet1->Phi, bJet2->Phi);
						Float_t dR1 = deltaR(bJet1->Eta, photon1->Eta, bJet1->Phi, photon1->Phi);
						Float_t dR2 = deltaR(bJet1->Eta, photon2->Eta, bJet1->Phi, photon2->Phi);
						Float_t dR3 = deltaR(bJet2->Eta, photon1->Eta, bJet2->Phi, photon1->Phi);
						Float_t dR4 = deltaR(bJet2->Eta, photon2->Eta, bJet2->Phi, photon2->Phi);

						dRPhotonBJet = dR1;
						if (dRPhotonBJet > dR2) dRPhotonBJet = dR2;
						if (dRPhotonBJet > dR3) dRPhotonBJet = dR3;
						if (dRPhotonBJet > dR4) dRPhotonBJet = dR4;
						
						// Check for kinematic requirements
						if((dRPhotons < 2) && (dRPhotonBJet > 1.5) && (dRBJets < 2)){
						
							tempKinCut += weight;
						
							// Check for mass cut requirement
							if((bJetSystem.M() > 105) && (bJetSystem.M() < 145)){
							
								tempMassCut += weight;
							
							}
							
						}
					
					}
				
				}
			
			}
			
		}
		
	} // end event loop
	
	// Update global event yield variables
	(SorB ? selectionSignal : selectionBackground) += crossSection*3000*tempSelection/numberOfEntries;
	(SorB ? kinCutBackground : kinCutBackground) += crossSection*3000*tempKinCut/numberOfEntries;
	(SorB ? massCutSignal : massCutBackground) += crossSection*3000*tempMassCut/numberOfEntries;

	ifs.close();
}

/*
 * SAVE RESULTS
 */
 
void saveResults(){
	
	cout << "Signal" << endl;
	cout << "Events after selection and mass window: " << selectionSignal << endl;
	cout << "Events after kinematic cut: " << kinCutSignal << endl;
	cout << "Events after mass cut: " << massCutSignal << endl;
	
	cout << "\n\n\nBackground" << endl;
	cout << "Events after selection and mass window: " << selectionBackground << endl;
	cout << "Events after kinematic cut: " << kinCutBackground << endl;
	cout << "Events after mass cut: " << massCutBackground << endl;
	
}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

	const Float_t pi = 3.14159265358979;

	Float_t etaDiff = (eta1-eta2);
	Float_t phiDiff = fabs(phi1-phi2);
	while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

	Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

	return TMath::Sqrt(deltaRSquared);

}



void histogram(TH1D *histoS, TH1D *histoB, TCanvas *can, const char* xTitle, const char* yTitle, const char* name){
	Double_t nS=1, nB=1;
	
	nS/=histoS->Integral();
	nB/=histoB->Integral();
	histoS->Scale(nS);
	histoB->Scale(nB);
	
	Double_t max;
	if((histoS->GetMaximum()) > (histoB->GetMaximum())){max=1.1*(histoS->GetMaximum());}
	else if((histoS->GetMaximum()) <= (histoB->GetMaximum())){max=1.1*(histoB->GetMaximum());}
	
	histoS->SetMaximum(max);
	histoB->SetMaximum(max);
	
	histoS->Draw();
	// add axis labels
	histoS->GetXaxis()->SetTitle(xTitle);
	histoS->GetYaxis()->SetTitle(yTitle);
	histoS->SetTitle(""); // title on top
	
	
	histoB->SetLineColor(kRed);
	histoB->Draw("same");
	

	can->SaveAs(name);
}

void histogram(TH1D *histo, TCanvas *can, const char* xTitle, const char* yTitle, const char* name){
	Double_t norm=1;
	norm/=histo->Integral();
	histo->Scale(norm);
	histo->Draw();
	// add axis labels
	histo->GetXaxis()->SetTitle(xTitle);
	histo->GetYaxis()->SetTitle(yTitle);
	histo->SetTitle(""); // title on top

	can->SaveAs(name);
}
