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
#include <TStyle.h>
#include <TLegend.h>
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

TString BackgroundSample;

// Declare functions
void histogram(TH1D*, TH1D*, TCanvas*, const char*, const char*, const char*);
void histogram(TH1D*, TCanvas*, const char*, const char*, const char*);
Float_t deltaR( const Float_t, const Float_t, const Float_t, const Float_t);
void analyze(TString, Double_t, bool);
void saveResults();
void saveResultsS();


// Initialize histograms
TH1D *hnLeptonsS = new TH1D("hnLeptonsS", "hnLeptonsS", 5, -0.5, 5);				TH1D *hnLeptonsB = new TH1D("hnLeptonsB", "hnLeptonsB", 5, -0.5, 5);
TH1D *hnJetsS = new TH1D("hnJetsS", "hnJetsS", 20, -0.5, 20);						TH1D *hnJetsB = new TH1D("hnJetsB", "hnJetsB", 20, -0.5, 20);
TH1D *hdRPhotonsS = new TH1D("hdRPhotonsS", "hdRPhotonsS", 25, 0, 5);				TH1D *hdRPhotonsB = new TH1D("hdRPhotonsB", "hdRPhotonsB", 25, 0, 5);
TH1D *hdRBJetsS = new TH1D("hdRBJetsS", "hdRBJetsS", 25, 0, 5);						TH1D *hdRBJetsB = new TH1D("hdRBJetsB", "hdRBJetsB", 25, 0, 5);
TH1D *hdRPhotonBJetS = new TH1D("hdRPhotonBJetS", "hdRPhotonBJetS", 25, 0, 5);		TH1D *hdRPhotonBJetB = new TH1D("hdRPhotonBJetB", "hdRPhotonBJetB", 25, 0, 5);
TH1D *hmPhotonsS = new TH1D("hmPhotonsS", "hmPhotonsS", 25, 100, 150);				TH1D *hmPhotonsB = new TH1D("hmPhotonsB", "hmPhotonsB", 25, 100, 150);
TH1D *hmBJetsS = new TH1D("hmBJetsS", "hmBJetsS", 50, 70, 170);						TH1D *hmBJetsB = new TH1D("hmBJetsB", "hmBJetsB", 50, 70, 170);
TH1D *hPhotonIsoS = new TH1D("hPhotonIsoS", "hPhotonIsoS", 50, -5, 25);				TH1D *hPhotonIsoB = new TH1D("hPhotonIsoB", "hPhotonIsoB", 50, -5, 25);


// Set up variables for event yields
Double_t selectionSignal=0, kinCutSignal=0, massCutSignal=0;
Double_t selectionBackground=0, kinCutBackground=0, massCutBackground=0;

Double_t errorSelectionSignal=0, errorKinCutSignal=0, errorMassCutSignal=0;
Double_t errorSelectionBackground=0, errorKinCutBackground=0, errorMassCutBackground=0;

/*
 * MAIN FUNCTION
 */

void bbgg(TString backgroundSample){
	
	BackgroundSample = backgroundSample;
	
	// Analyze signal
	analyze("HHToGGBB_14TeV", 0.089, Signal);
	
	if(BackgroundSample == "B"){
		// Analyze B background
		analyze("B-4p-0-1-v1510_14TEV", 200944.68129*1000, Background);
	}
	
	else if(BackgroundSample == "BB"){
		// Analyse BB background
		analyze("BB-4p-0-300-v1510_14TEV", 249.97710*1000, Background);
		analyze("BB-4p-300-700-v1510_14TEV", 35.23062*1000, Background);
		analyze("BB-4p-700-1300-v1510_14TEV", 4.13743*1000, Background);
		analyze("BB-4p-1300-2100-v1510_14TEV", 0.41702*1000, Background);
		analyze("BB-4p-2100-100000-v1510_14TEV", 0.04770*1000, Background);
	}
	
	else if(BackgroundSample == "BBB"){
		// Analyze BBB background
		analyze("BBB-4p-0-600-v1510_14TEV", 2.57304*1000, Background);
		analyze("BBB-4p-600-1300-v1510_14TEV", 0.14935*1000, Background);
		analyze("BBB-4p-1300-100000-v1510_14TEV", 0.01274*1000, Background);
	}
	
	else if(BackgroundSample == "Bj"){
		// Analyze Bj background
		analyze("Bj-4p-0-300-v1510_14TEV", 34409.92339*1000, Background);
		analyze("Bj-4p-300-600-v1510_14TEV", 2642.85309*1000, Background);
		analyze("Bj-4p-600-1100-v1510_14TEV", 294.12311*1000, Background);
		analyze("Bj-4p-1100-1800-v1510_14TEV", 25.95000*1000, Background);
		analyze("Bj-4p-1800-2700-v1510_14TEV", 2.42111*1000, Background);
		analyze("Bj-4p-2700-3700-v1510_14TEV", 0.22690*1000, Background);
		analyze("Bj-4p-3700-100000-v1510_14TEV", 0.02767*1000, Background);
	}
	
	else if(BackgroundSample == "Bjj"){
		// Analyze Bjj-vbf background
		analyze("Bjj-vbf-4p-0-700-v1510_14TEV", 86.45604*1000, Background);
		analyze("Bjj-vbf-4p-700-1400-v1510_14TEV", 4.34869*1000, Background);
		analyze("Bjj-vbf-4p-1400-2300-v1510_14TEV", 0.32465*1000, Background);
		analyze("Bjj-vbf-4p-2300-3400-v1510_14TEV", 0.03032*1000, Background);
		//analyze("Bjj-vbf-4p-3400-100000-v1510_14TEV", 0.00313*1000, Background);
	}
	
	else if(BackgroundSample == "H"){
		// Analyze H background
		analyze("H-4p-0-300-v1510_14TEV", 21.55990*1000, Background);
		analyze("H-4p-300-800-v1510_14TEV", 1.11282*1000, Background);
		analyze("H-4p-800-1500-v1510_14TEV", 0.09188*1000, Background);
		analyze("H-4p-1500-100000-v1510_14TEV", 0.01009*1000, Background);
	}
	
	else if(BackgroundSample == "LL"){
		// Analyze LL background
		analyze("LL-4p-0-100-v1510_14TEV", 1341.36923*1000, Background);
		analyze("LL-4p-100-200-v1510_14TEV", 156.29534*1000, Background);
		analyze("LL-4p-200-500-v1510_14TEV", 42.40132*1000, Background);
		analyze("LL-4p-500-900-v1510_14TEV", 2.84373*1000, Background);
		analyze("LL-4p-900-1400-v1510_14TEV", 0.20914*1000, Background);
		analyze("LL-4p-1400-100000-v1510_14TEV", 0.02891*1000, Background);
	}
	
	else if(BackgroundSample == "LLB"){
		// Analyze LLB background
		analyze("LLB-4p-0-400-v1510_14TEV", 2.97380*1000, Background);
		//analyze("LLB-4p-400-900-v1510_14TEV", 0.22854*1000, Background);
		analyze("LLB-4p-900-100000-v1510_14TEV", 0.02080*1000, Background);
	}
	
	else if(BackgroundSample == "tB"){
		// Analyze tB background
		analyze("tB-4p-0-500-v1510_14TEV", 63.88923*1000, Background);
		analyze("tB-4p-500-900-v1510_14TEV", 7.12172*1000, Background);
		analyze("tB-4p-900-1500-v1510_14TEV", 0.98030*1000, Background);
		analyze("tB-4p-1500-2200-v1510_14TEV", 0.08391*1000, Background);
		analyze("tB-4p-2200-100000-v1510_14TEV", 0.00953*1000, Background);
	}
	
	else if(BackgroundSample == "tj"){
		// Analyze tj background
		analyze("tj-4p-0-500-v1510_14TEV", 109.73602*1000, Background);
		analyze("tj-4p-500-1000-v1510_14TEV", 5.99325*1000, Background);
		analyze("tj-4p-1000-1600-v1510_14TEV", 0.37680*1000, Background);
		analyze("tj-4p-1600-2400-v1510_14TEV", 0.03462*1000, Background);
		analyze("tj-4p-2400-100000-v1510_14TEV", 0.00312*1000, Background);
	}
	
	else if(BackgroundSample == "tt"){
		// Analyse tt background
		analyze("tt-4p-0-600-v1510_14TEV", 530.89358*1000, Background);
		analyze("tt-4p-600-1100-v1510_14TEV", 42.55351*1000, Background);
		analyze("tt-4p-1100-1700-v1510_14TEV", 4.48209*1000, Background);
		analyze("tt-4p-1700-2500-v1510_14TEV", 0.52795*1000, Background);
		analyze("tt-4p-2500-100000-v1510_14TEV", 0.05449*1000, Background);
	}
	
	else if(BackgroundSample == "ttB"){
		// Analyze ttB background
		analyze("ttB-4p-0-900-v1510_14TEV", 2.6673*1000, Background);
		analyze("ttB-4p-900-1600-v1510_14TEV", 0.250469*1000, Background);
		analyze("ttB-4p-1600-2500-v1510_14TEV", 0.0237441*1000, Background);
		analyze("ttB-4p-2500-100000-v1510_14TEV", 0.00208816*1000, Background);
	}
	
	else if(BackgroundSample == "Signal"){
		cout << "Processing only signal events" << endl;
	}
	
	else {
		cout << "Background sample not found" << endl;
		assert(false);
	}
	
	// Save results
	if(BackgroundSample == "Signal") saveResultsS();
	else saveResults();
	
}


/*
 * ANALYSIS
 */

void analyze(TString inputfile, Double_t crossSection, bool SorB)
{	
	
	// Set up temporal vatiables
	Double_t tempSelection=0, tempKinCut=0, tempMassCut=0;
	
	Double_t tempErrorSelection=0, tempErrorKinCut=0, tempErrorMassCut=0;
	
	const TString inputFilet = inputfile + ".txt";
	ifstream ifs(inputFilet);

	assert(ifs.is_open());

	TString filename;
	
	cout << "Reading from folder " << inputfile << endl;

	TChain chain("Delphes");

	while(ifs >> filename){
		
		TString filenamef;
		if(SorB) filenamef = "root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/" + inputfile + "/" + filename;
		else filenamef = "root://eoscms.cern.ch//store/group/upgrade/delphes/ProdJun14/" + inputfile + "/" + filename;
		
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
	//((!SorB) && (BackgroundSample != "tt") ? branchEvent = treeReader->UseBranch("Event"): branchEvent = 0);
			
	for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // Event loop
		treeReader->ReadEntry(iEntry);		
			
		// Reset index variables
		nPhoton1=nPhoton2=nbJet1=nbJet2=-1;
		
		// Store the number of leptons
		nLeptons = branchElectron->GetEntries() + branchMuon->GetEntries();
		nJets = 0;
			
		// Select photons
		for (Int_t iPhoton=0; iPhoton<branchPhoton->GetEntries(); iPhoton++) { // Photon loop
			photon = (Photon*) branchPhoton->At(iPhoton);
		
			(SorB ? hPhotonIsoS : hPhotonIsoB)->Fill(photon->IsolationVar);
		
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
		
		weight = 1;
		/*if(!SorB){
						
			event = (LHEFEvent*) branchEvent->At(0);
			weight = event->Weight;
						
		}*/

		(SorB ? hnLeptonsS : hnLeptonsB)->Fill(nLeptons, weight);
		(SorB ? hnJetsS : hnJetsB)->Fill(nJets, weight);

		// Check if there are two selected photons and bjets
		if(!((nPhoton1!=-1) && (nPhoton2!=-1) && (nbJet1!=-1) && (nbJet2!=-1))) continue;
		
		(SorB ? hPhotonIsoS : hPhotonIsoB)->Fill(photon1->IsolationVar, weight);
		(SorB ? hPhotonIsoS : hPhotonIsoB)->Fill(photon2->IsolationVar, weight);
		
		
		if(!((photon1->PT > 40) && (nLeptons == 0) && (nJets < 4))) continue;
			
		// Set up four-vectors
		vPhoton1.SetPtEtaPhiE(photon1->PT, photon1->Eta, photon1->Phi, photon1->E);
		vPhoton2.SetPtEtaPhiE(photon2->PT, photon2->Eta, photon2->Phi, photon2->E);
		vBJet1.SetPtEtaPhiM(bJet1->PT, bJet1->Eta, bJet1->Phi, bJet1->Mass);
		vBJet2.SetPtEtaPhiM(bJet2->PT, bJet2->Eta, bJet2->Phi, bJet2->Mass);
	
		photonSystem = vPhoton1 + vPhoton2;
		bJetSystem = vBJet1 + vBJet2;
		
		// Check for fit mass window requirements
		if((photonSystem.M() > 120) && (photonSystem.M() < 130) && (bJetSystem.M() > 70) && (bJetSystem.M() < 200)){
					
			tempSelection += weight;
			tempErrorSelection += weight*weight;
		}
		
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
		
		(SorB ? hdRPhotonsS : hdRPhotonsB)->Fill(dRPhotons, weight);
		(SorB ? hdRPhotonBJetS : hdRPhotonBJetB)->Fill(dRPhotonBJet, weight);
		
		// Check for kinematic requirements
		if(!((dRPhotons < 2) && (dRPhotonBJet > 1.5))) continue;
		
		(SorB ? hdRBJetsS : hdRBJetsB)->Fill(dRBJets, weight);
		
		if(!(dRBJets < 2)) continue;
		
		(SorB ? hmPhotonsS : hmPhotonsB)->Fill(photonSystem.M(), weight);
		(SorB ? hmBJetsS : hmBJetsB)->Fill(bJetSystem.M(), weight);
		
		if(!((photonSystem.M() > 120) && (photonSystem.M() < 130) && (bJetSystem.M() > 70) && (bJetSystem.M() < 200))) continue;
					
		tempKinCut += weight;
		tempErrorKinCut += weight*weight;
	
		// Check for mass cut requirement
		if(!((bJetSystem.M() > 105) && (bJetSystem.M() < 145))) continue;
		
		tempMassCut += weight;
		tempErrorMassCut += weight*weight;
		
	} // end event loop
	
	// Update global event yield variables
	(SorB ? selectionSignal : selectionBackground) += crossSection*3000*tempSelection/numberOfEntries;
	(SorB ? kinCutSignal : kinCutBackground) += crossSection*3000*tempKinCut/numberOfEntries;
	(SorB ? massCutSignal : massCutBackground) += crossSection*3000*tempMassCut/numberOfEntries;
	
	(SorB ? errorSelectionSignal : errorSelectionBackground) += crossSection*3000*sqrtf(tempErrorSelection)/numberOfEntries;
	(SorB ? errorKinCutSignal : errorKinCutBackground) += crossSection*3000*sqrtf(tempErrorKinCut)/numberOfEntries;
	(SorB ? errorMassCutSignal : errorMassCutBackground) += crossSection*3000*sqrtf(tempErrorMassCut)/numberOfEntries;

	ifs.close();
}

/*
 * SAVE RESULTS
 */
 
void saveResults(){
	
	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);
	
	gStyle->SetOptStat(kFALSE);
	
	histogram(hnLeptonsS, hnLeptonsB, c1, "Number of Leptons", "Count", "./Histograms/histogramnL_" + BackgroundSample + ".jpg");
	histogram(hnJetsS, hnJetsB, c1, "Number of Jets with |Eta|<2.5", "Count", "./Histograms/histogramnJ_" + BackgroundSample + ".jpg");
	histogram(hdRPhotonsS, hdRPhotonsB, c1, "delta R of Photons", "Count", "./Histograms/histogramdRP_" + BackgroundSample + ".jpg");
	histogram(hdRBJetsS, hdRBJetsB, c1, "delta R of BJets", "Count", "./Histograms/histogramdRBJ_" + BackgroundSample + ".jpg");
	histogram(hdRPhotonBJetS, hdRPhotonBJetB, c1, "min delta R of Photon and BJet", "Count", "./Histograms/histogramdRPBJ_" + BackgroundSample + ".jpg");
	histogram(hmPhotonsS, hmPhotonsB, c1, "Reconstructed mass from Photons", "Count", "./Histograms/histogramdmP_" + BackgroundSample + ".jpg");
	histogram(hmBJetsS, hmBJetsB, c1, "Reconstructed mass from BJets", "Count", "./Histograms/histogramdmBJ_" + BackgroundSample + ".jpg");
	histogram(hPhotonIsoS, hPhotonIsoB, c1, "Photon isolation variable", "Count", "./Histograms/histogramPhotonIso_" + BackgroundSample + ".jpg");
	
	cout << "\nSignal\n" << endl;
	cout << "Events after selection and mass window: " << selectionSignal << " +- " << errorSelectionSignal << endl;
	cout << "Events after kinematic cut: " << kinCutSignal << " +- " << errorKinCutSignal << endl;
	cout << "Events after mass cut: " << massCutSignal << " +- " << errorMassCutSignal << endl;
	
	cout << "\n\n" << BackgroundSample << " Background\n" << endl;
	cout << "Events after selection and mass window: " << selectionBackground << " +- " << errorSelectionBackground << endl;
	cout << "Events after kinematic cut: " << kinCutBackground << " +- " << errorKinCutBackground << endl;
	cout << "Events after mass cut: " << massCutBackground << " +- " << errorMassCutBackground << endl;
	
}

void saveResultsS(){
	
	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);
	
	gStyle->SetOptStat(kFALSE);
	
	histogram(hnLeptonsS, c1, "Number of Leptons", "Count", "./Histograms/histogramnL_Signal.jpg");
	histogram(hnJetsS, c1, "Number of Jets with |Eta|<2.5", "Count", "./Histograms/histogramnJ_Signal.jpg");
	histogram(hdRPhotonsS, c1, "delta R of Photons", "Count", "./Histograms/histogramdRP_Signal.jpg");
	histogram(hdRBJetsS, hdRBJetsB, c1, "delta R of BJets", "Count", "./Histograms/histogramdRBJ_Signal.jpg");
	histogram(hdRPhotonBJetS, c1, "min delta R of Photon and BJet", "Count", "./Histograms/histogramdRPBJ_Signal.jpg");
	histogram(hmPhotonsS, c1, "Reconstructed mass from Photons", "Count", "./Histograms/histogramdmP_Signal.jpg");
	histogram(hmBJetsS, c1, "Reconstructed mass from BJets", "Count", "./Histograms/histogramdmBJ_Signal.jpg");
	histogram(hPhotonIsoS, c1, "Photon isolation variable", "Count", "./Histograms/histogramPhotonIso_Signal.jpg");
	
	cout << "\nSignal\n" << endl;
	cout << "Events after selection and mass window: " << selectionSignal << " +- " << errorSelectionSignal << endl;
	cout << "Events after kinematic cut: " << kinCutSignal << " +- " << errorKinCutSignal << endl;
	cout << "Events after mass cut: " << massCutSignal << " +- " << errorMassCutSignal << endl;
	
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
	
	/*TLegend *leg = new TLegend();
	leg->AddEntry(histoS,"Signal","f");
	leg->AddEntry(histoB,BackgroundSample + " background","l");
	leg->Draw();*/


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
