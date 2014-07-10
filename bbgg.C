/*
 * bbgg.C
 * 
 * This macro analyses bbgg signal and background samples.
 * More specifically, it analyses the ttH and ZH backgrounds.
 * To run use root bbgg.C+\(\"name of background\"\)
 * 
 * The event selection consists in requiring a leading photon with pt > 90,
 * and a subleading photon with pt > 40, both with |eta| < 2.5.
 * It also requires two btagged jets with |eta| < 2.4, the leading
 * bjet with pt > 60 and the subleading with pt > 40.
 * Finally it requires no reconstructed leptons and less than four jets
 * with pt > 30 and |eta| < 2.5
 * 
 * The fit mass window consists in requiring a reconstructed mass
 * 120 < M_gg < 130 and 70 < M_bb < 200
 * 
 * The kinematic cut requires deltaR_gg < 2, deltaR_bb < 2 and that
 * min(deltaR_bg) > 1.5
 * 
 * The mass cut requires 120 < M_gg < 130,  105 < M_bb < 145
 * and 300 < M_HH < 900
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

int puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t betastar);

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
TH1D *hmBJetsS = new TH1D("hmBJetsS", "hmBJetsS", 250, -5, 5);						TH1D *hmBJetsB = new TH1D("hmBJetsB", "hmBJetsB", 250, -5, 5);


// Set up variables for event yields
Double_t totalSignal=0, selectionSignal=0, kinCutSignal=0, massCutSignal=0;
Double_t totalBackground=0, selectionBackground=0, kinCutBackground=0, massCutBackground=0;

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
	
	else if(BackgroundSample == "all"){
		
		// Analyze B background
		analyze("B-4p-0-1-v1510_14TEV", 200944.68129*1000, Background);

		// Analyse BB background
		analyze("BB-4p-0-300-v1510_14TEV", 249.97710*1000, Background);
		analyze("BB-4p-300-700-v1510_14TEV", 35.23062*1000, Background);
		analyze("BB-4p-700-1300-v1510_14TEV", 4.13743*1000, Background);
		analyze("BB-4p-1300-2100-v1510_14TEV", 0.41702*1000, Background);
		analyze("BB-4p-2100-100000-v1510_14TEV", 0.04770*1000, Background);

		// Analyze BBB background
		analyze("BBB-4p-0-600-v1510_14TEV", 2.57304*1000, Background);
		analyze("BBB-4p-600-1300-v1510_14TEV", 0.14935*1000, Background);
		analyze("BBB-4p-1300-100000-v1510_14TEV", 0.01274*1000, Background);

		// Analyze Bj background
		analyze("Bj-4p-0-300-v1510_14TEV", 34409.92339*1000, Background);
		analyze("Bj-4p-300-600-v1510_14TEV", 2642.85309*1000, Background);
		analyze("Bj-4p-600-1100-v1510_14TEV", 294.12311*1000, Background);
		analyze("Bj-4p-1100-1800-v1510_14TEV", 25.95000*1000, Background);
		analyze("Bj-4p-1800-2700-v1510_14TEV", 2.42111*1000, Background);
		analyze("Bj-4p-2700-3700-v1510_14TEV", 0.22690*1000, Background);
		analyze("Bj-4p-3700-100000-v1510_14TEV", 0.02767*1000, Background);

		// Analyze Bjj-vbf background
		analyze("Bjj-vbf-4p-0-700-v1510_14TEV", 86.45604*1000, Background);
		analyze("Bjj-vbf-4p-700-1400-v1510_14TEV", 4.34869*1000, Background);
		analyze("Bjj-vbf-4p-1400-2300-v1510_14TEV", 0.32465*1000, Background);
		analyze("Bjj-vbf-4p-2300-3400-v1510_14TEV", 0.03032*1000, Background);
		//analyze("Bjj-vbf-4p-3400-100000-v1510_14TEV", 0.00313*1000, Background);

		// Analyze H background
		analyze("H-4p-0-300-v1510_14TEV", 21.55990*1000, Background);
		analyze("H-4p-300-800-v1510_14TEV", 1.11282*1000, Background);
		analyze("H-4p-800-1500-v1510_14TEV", 0.09188*1000, Background);
		analyze("H-4p-1500-100000-v1510_14TEV", 0.01009*1000, Background);

		// Analyze LL background
		analyze("LL-4p-0-100-v1510_14TEV", 1341.36923*1000, Background);
		analyze("LL-4p-100-200-v1510_14TEV", 156.29534*1000, Background);
		analyze("LL-4p-200-500-v1510_14TEV", 42.40132*1000, Background);
		analyze("LL-4p-500-900-v1510_14TEV", 2.84373*1000, Background);
		analyze("LL-4p-900-1400-v1510_14TEV", 0.20914*1000, Background);
		analyze("LL-4p-1400-100000-v1510_14TEV", 0.02891*1000, Background);

		// Analyze LLB background
		analyze("LLB-4p-0-400-v1510_14TEV", 2.97380*1000, Background);
		//analyze("LLB-4p-400-900-v1510_14TEV", 0.22854*1000, Background);
		analyze("LLB-4p-900-100000-v1510_14TEV", 0.02080*1000, Background);

		// Analyze tB background
		analyze("tB-4p-0-500-v1510_14TEV", 63.88923*1000, Background);
		analyze("tB-4p-500-900-v1510_14TEV", 7.12172*1000, Background);
		analyze("tB-4p-900-1500-v1510_14TEV", 0.98030*1000, Background);
		analyze("tB-4p-1500-2200-v1510_14TEV", 0.08391*1000, Background);
		analyze("tB-4p-2200-100000-v1510_14TEV", 0.00953*1000, Background);

		// Analyze tj background
		analyze("tj-4p-0-500-v1510_14TEV", 109.73602*1000, Background);
		analyze("tj-4p-500-1000-v1510_14TEV", 5.99325*1000, Background);
		analyze("tj-4p-1000-1600-v1510_14TEV", 0.37680*1000, Background);
		analyze("tj-4p-1600-2400-v1510_14TEV", 0.03462*1000, Background);
		analyze("tj-4p-2400-100000-v1510_14TEV", 0.00312*1000, Background);

		// Analyse tt background
		analyze("tt-4p-0-600-v1510_14TEV", 530.89358*1000, Background);
		analyze("tt-4p-600-1100-v1510_14TEV", 42.55351*1000, Background);
		analyze("tt-4p-1100-1700-v1510_14TEV", 4.48209*1000, Background);
		analyze("tt-4p-1700-2500-v1510_14TEV", 0.52795*1000, Background);
		analyze("tt-4p-2500-100000-v1510_14TEV", 0.05449*1000, Background);

		// Analyze ttB background
		analyze("ttB-4p-0-900-v1510_14TEV", 2.6673*1000, Background);
		analyze("ttB-4p-900-1600-v1510_14TEV", 0.250469*1000, Background);
		analyze("ttB-4p-1600-2500-v1510_14TEV", 0.0237441*1000, Background);
		analyze("ttB-4p-2500-100000-v1510_14TEV", 0.00208816*1000, Background);
		
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
	
	const TString inputFile = "/afs/cern.ch/work/a/ariostas/public/bbgg/" + inputfile + ".root";
	
	cout << "Reading from file " << inputfile << endl;
	
	// Set up storage variables
	Photon *photon1=0, *photon2=0;
	Jet *bJet1=0, *bJet2=0;
	UInt_t nJets, nLeptons; 
	Double_t weight, dRPhotons, dRBJets, dRPhotonBJet;
	TLorentzVector vPhoton1, vPhoton2, vBJet1, vBJet2, photonSystem, bJetSystem, diHiggsSystem;

	TFile* infile = new TFile(inputFile); assert(infile);
	TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

	intree->SetBranchAddress("weight",		&weight);
	intree->SetBranchAddress("nJets",			&nJets);
	intree->SetBranchAddress("nLeptons",		&nLeptons);
	intree->SetBranchAddress("photon1",		&photon1);
	intree->SetBranchAddress("photon2",		&photon2);
	intree->SetBranchAddress("bJet1",          &bJet1);
	intree->SetBranchAddress("bJet2",          &bJet2);

			
	for (Long64_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // Event loop
		intree->GetEntry(iEntry);

		(SorB ? hnLeptonsS : hnLeptonsB)->Fill(nLeptons, weight);
		(SorB ? hnJetsS : hnJetsB)->Fill(nJets, weight);
		
		if(!((photon1->PT > 90) && (photon2->PT > 40) && (bJet1->PT > 60) && (bJet2->PT > 40) && (nLeptons == 0) && (nJets < 4))) continue;
			
		// Set up four-vectors
		vPhoton1.SetPtEtaPhiE(photon1->PT, photon1->Eta, photon1->Phi, photon1->E);
		vPhoton2.SetPtEtaPhiE(photon2->PT, photon2->Eta, photon2->Phi, photon2->E);
		vBJet1.SetPtEtaPhiM(bJet1->PT, bJet1->Eta, bJet1->Phi, bJet1->Mass);
		vBJet2.SetPtEtaPhiM(bJet2->PT, bJet2->Eta, bJet2->Phi, bJet2->Mass);
	
		photonSystem = vPhoton1 + vPhoton2;
		bJetSystem = vBJet1 + vBJet2;
		
		diHiggsSystem = vPhoton1 + vPhoton2 + vBJet1 + vBJet2;
		
		// Check for fit mass window requirements
		if((photonSystem.M() > 120) && (photonSystem.M() < 130) && (bJetSystem.M() > 70) && (bJetSystem.M() < 200)){
					
			tempSelection += weight;
			tempErrorSelection++;
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
		(SorB ? hmBJetsS : hmBJetsB)->Fill(diHiggsSystem.M(), weight);
		
		if(!((photonSystem.M() > 120) && (photonSystem.M() < 130) && (bJetSystem.M() > 70) && (bJetSystem.M() < 200))) continue;
					
		tempKinCut += weight;
		tempErrorKinCut++;
	
		// Check for mass cut requirement
		if(!((bJetSystem.M() > 105) && (bJetSystem.M() < 145) && (diHiggsSystem.M() > 300) && (diHiggsSystem.M() < 900))) continue;
		
		tempMassCut += weight;
		tempErrorMassCut++;
		
	} // end event loop
	
	// Update global event yield variables
	(SorB ? totalSignal : totalBackground) += 3000*crossSection;
	(SorB ? selectionSignal : selectionBackground) += tempSelection;
	(SorB ? kinCutSignal : kinCutBackground) += tempKinCut;
	(SorB ? massCutSignal : massCutBackground) += tempMassCut;
	
	(SorB ? errorSelectionSignal : errorSelectionBackground) += sqrtf(tempErrorSelection)*weight;
	(SorB ? errorKinCutSignal : errorKinCutBackground) += sqrtf(tempErrorKinCut)*weight;
	(SorB ? errorMassCutSignal : errorMassCutBackground) += sqrtf(tempErrorMassCut)*weight;

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
	
	cout << "\nSignal\n" << endl;
	cout << "Total Events " << totalSignal << endl;
	cout << "Events after selection and mass window: " << selectionSignal << " +- " << errorSelectionSignal << endl;
	cout << "Events after kinematic cut: " << kinCutSignal << " +- " << errorKinCutSignal << endl;
	cout << "Events after mass cut: " << massCutSignal << " +- " << errorMassCutSignal << endl;
	
	cout << "\n\n" << BackgroundSample << " Background\n" << endl;
	cout << "Total Events " << totalBackground << endl;
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
	
	cout << "\nSignal\n" << endl;
	cout << "Total Events " << totalSignal << endl;
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
	if((histoS->GetMaximum()) > (histoB->GetMaximum())) max=1.1*(histoS->GetMaximum());
	else max=1.1*(histoB->GetMaximum());
	
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
