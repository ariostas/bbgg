/*
 * bbgg.C
 * 
 * This macro analyses bbgg signal and background samples.
 * To run use  root -b -q bbgg.C+  for all backgrounds and with the High-Pt
 * category or root -b -q bbgg.C+\(\"background\",\"category\"\) 
 * for a specific background and category.
 * 
 * Code by: Andres Rios
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

TString BackgroundSample, Category;

// Declare functions
void histogram(TH1D*, TH1D*, TCanvas*, const char*, const char*, const char*);
void histogram(TH1D*, TCanvas*, const char*, const char*, const char*);
Float_t deltaR( const Float_t, const Float_t, const Float_t, const Float_t);
void analyze(TString, Double_t, bool, TString);
void analyzeQCD(const TString, const Double_t, TString);
void saveResults();
void saveResultsS();


// Initialize histograms
TH1D *hnLeptonsS = new TH1D("hnLeptonsS", "hnLeptonsS", 5, -0.5, 5);				TH1D *hnLeptonsB = new TH1D("hnLeptonsB", "hnLeptonsB", 5, -0.5, 5);
TH1D *hnJetsS = new TH1D("hnJetsS", "hnJetsS", 20, -0.5, 20);						TH1D *hnJetsB = new TH1D("hnJetsB", "hnJetsB", 20, -0.5, 20);
TH1D *hdRPhotonsS = new TH1D("hdRPhotonsS", "hdRPhotonsS", 250, 0, 5);				TH1D *hdRPhotonsB = new TH1D("hdRPhotonsB", "hdRPhotonsB", 250, 0, 5);
TH1D *hdRBJetsS = new TH1D("hdRBJetsS", "hdRBJetsS", 250, 0, 5);					TH1D *hdRBJetsB = new TH1D("hdRBJetsB", "hdRBJetsB", 250, 0, 5);
TH1D *hdRPhotonBJetS = new TH1D("hdRPhotonBJetS", "hdRPhotonBJetS", 250, 0, 5);		TH1D *hdRPhotonBJetB = new TH1D("hdRPhotonBJetB", "hdRPhotonBJetB", 250, 0, 5);
TH1D *hmPhotonsS = new TH1D("hmPhotonsS", "hmPhotonsS", 250, 100, 150);				TH1D *hmPhotonsB = new TH1D("hmPhotonsB", "hmPhotonsB", 250, 100, 150);
TH1D *hmBJetsS = new TH1D("hmBJetsS", "hmBJetsS", 250, 0, 200);						TH1D *hmBJetsB = new TH1D("hmBJetsB", "hmBJetsB", 250, 0, 200);
TH1D *hmDiHiggsS = new TH1D("hmDiHiggsS", "hmDiHiggsS", 250, 200, 1300);			TH1D *hmDiHiggsB = new TH1D("hmDiHiggsB", "hmDiHiggsB", 250, 200, 1300);

TH1D *hPhoton1PtS = new TH1D("hPhoton1PtS", "hPhoton1PtS", 150, 0, 150);			TH1D *hPhoton1PtB = new TH1D("hPhoton1PtB", "hPhoton1PtB", 150, 0, 150);
TH1D *hPhoton2PtS = new TH1D("hPhoton2PtS", "hPhoton2PtS", 150, 0, 150);			TH1D *hPhoton2PtB = new TH1D("hPhoton2PtB", "hPhoton2PtB", 150, 0, 150);
TH1D *hBJet1PtS = new TH1D("hBJet1PtS", "hBJet1PtS", 150, 0, 150);					TH1D *hBJet1PtB = new TH1D("hBJet1PtB", "hBJet1PtB", 150, 0, 150);
TH1D *hBJet2PtS = new TH1D("hBJet2PtS", "hBJet2PtS", 150, 0, 150);					TH1D *hBJet2PtB = new TH1D("hBJet2PtB", "hBJet2PtB", 150, 0, 150);

// Set up variables for event yields
Double_t totalSignal=0, selectionSignal=0, kinCutSignal=0, massCutSignal=0;
Double_t totalBackground=0, selectionBackground=0, kinCutBackground=0, massCutBackground=0;

Double_t errorSelectionSignal=0, errorKinCutSignal=0, errorMassCutSignal=0;
Double_t errorSelectionBackground=0, errorKinCutBackground=0, errorMassCutBackground=0;

/*
 * MAIN FUNCTION
 */

void bbgg(TString backgroundSample = "All", TString category = "HighPt"){
	
	BackgroundSample = backgroundSample;
	Category = category;
	
	if(Category != "HighPt" && Category != "LowPt") return;
		
	// Analyze signal
	analyze("HHToGGBB_14TeV", 40*0.561*0.00228*2, Signal, Category);
	
	if(BackgroundSample == "B"){
		// Analyze B background
		analyze("B-4p-0-1-v1510_14TEV", 200944.68129*1000, Background, Category);
	}
	
	else if(BackgroundSample == "BB"){
		// Analyse BB background
		analyze("BB-4p-0-300-v1510_14TEV", 249.97710*1000, Background, Category);
		analyze("BB-4p-300-700-v1510_14TEV", 35.23062*1000, Background, Category);
		analyze("BB-4p-700-1300-v1510_14TEV", 4.13743*1000, Background, Category);
		analyze("BB-4p-1300-2100-v1510_14TEV", 0.41702*1000, Background, Category);
		analyze("BB-4p-2100-100000-v1510_14TEV", 0.04770*1000, Background, Category);
	}
	
	else if(BackgroundSample == "BBB"){
		// Analyze BBB background
		analyze("BBB-4p-0-600-v1510_14TEV", 2.57304*1000, Background, Category);
		analyze("BBB-4p-600-1300-v1510_14TEV", 0.14935*1000, Background, Category);
		analyze("BBB-4p-1300-100000-v1510_14TEV", 0.01274*1000, Background, Category);
	}
	
	else if(BackgroundSample == "Bj"){
		// Analyze Bj background
		analyze("Bj-4p-0-300-v1510_14TEV", 34409.92339*1000, Background, Category);
		analyze("Bj-4p-300-600-v1510_14TEV", 2642.85309*1000, Background, Category);
		analyze("Bj-4p-600-1100-v1510_14TEV", 294.12311*1000, Background, Category);
		analyze("Bj-4p-1100-1800-v1510_14TEV", 25.95000*1000, Background, Category);
		analyze("Bj-4p-1800-2700-v1510_14TEV", 2.42111*1000, Background, Category);
		analyze("Bj-4p-2700-3700-v1510_14TEV", 0.22690*1000, Background, Category);
		analyze("Bj-4p-3700-100000-v1510_14TEV", 0.02767*1000, Background, Category);
	}
	
	else if(BackgroundSample == "Bjj"){
		// Analyze Bjj-vbf background
		analyze("Bjj-vbf-4p-0-700-v1510_14TEV", 86.45604*1000, Background, Category);
		analyze("Bjj-vbf-4p-700-1400-v1510_14TEV", 4.34869*1000, Background, Category);
		analyze("Bjj-vbf-4p-1400-2300-v1510_14TEV", 0.32465*1000, Background, Category);
		analyze("Bjj-vbf-4p-2300-3400-v1510_14TEV", 0.03032*1000, Background, Category);
		//analyze("Bjj-vbf-4p-3400-100000-v1510_14TEV", 0.00313*1000, Background, Category);
	}
	
	else if(BackgroundSample == "H"){
		// Analyze H background
		analyze("H-4p-0-300-v1510_14TEV", 21.55990*1000, Background, Category);
		analyze("H-4p-300-800-v1510_14TEV", 1.11282*1000, Background, Category);
		analyze("H-4p-800-1500-v1510_14TEV", 0.09188*1000, Background, Category);
		analyze("H-4p-1500-100000-v1510_14TEV", 0.01009*1000, Background, Category);
	}
	
	else if(BackgroundSample == "LL"){
		// Analyze LL background
		analyze("LL-4p-0-100-v1510_14TEV", 1341.36923*1000, Background, Category);
		analyze("LL-4p-100-200-v1510_14TEV", 156.29534*1000, Background, Category);
		analyze("LL-4p-200-500-v1510_14TEV", 42.40132*1000, Background, Category);
		analyze("LL-4p-500-900-v1510_14TEV", 2.84373*1000, Background, Category);
		analyze("LL-4p-900-1400-v1510_14TEV", 0.20914*1000, Background, Category);
		analyze("LL-4p-1400-100000-v1510_14TEV", 0.02891*1000, Background, Category);
	}
	
	else if(BackgroundSample == "LLB"){
		// Analyze LLB background
		analyze("LLB-4p-0-400-v1510_14TEV", 2.97380*1000, Background, Category);
		//analyze("LLB-4p-400-900-v1510_14TEV", 0.22854*1000, Background, Category);
		analyze("LLB-4p-900-100000-v1510_14TEV", 0.02080*1000, Background, Category);
	}
	
	else if(BackgroundSample == "tB"){
		// Analyze tB background
		analyze("tB-4p-0-500-v1510_14TEV", 63.88923*1000, Background, Category);
		analyze("tB-4p-500-900-v1510_14TEV", 7.12172*1000, Background, Category);
		analyze("tB-4p-900-1500-v1510_14TEV", 0.98030*1000, Background, Category);
		analyze("tB-4p-1500-2200-v1510_14TEV", 0.08391*1000, Background, Category);
		analyze("tB-4p-2200-100000-v1510_14TEV", 0.00953*1000, Background, Category);
	}
	
	else if(BackgroundSample == "tj"){
		// Analyze tj background
		analyze("tj-4p-0-500-v1510_14TEV", 109.73602*1000, Background, Category);
		analyze("tj-4p-500-1000-v1510_14TEV", 5.99325*1000, Background, Category);
		analyze("tj-4p-1000-1600-v1510_14TEV", 0.37680*1000, Background, Category);
		analyze("tj-4p-1600-2400-v1510_14TEV", 0.03462*1000, Background, Category);
		analyze("tj-4p-2400-100000-v1510_14TEV", 0.00312*1000, Background, Category);
	}
	
	else if(BackgroundSample == "tt"){
		// Analyse tt background
		analyze("tt-4p-0-600-v1510_14TEV", 530.89358*1000, Background, Category);
		analyze("tt-4p-600-1100-v1510_14TEV", 42.55351*1000, Background, Category);
		analyze("tt-4p-1100-1700-v1510_14TEV", 4.48209*1000, Background, Category);
		analyze("tt-4p-1700-2500-v1510_14TEV", 0.52795*1000, Background, Category);
		analyze("tt-4p-2500-100000-v1510_14TEV", 0.05449*1000, Background, Category);
	}
	
	else if(BackgroundSample == "ttB"){
		// Analyze ttB background
		analyze("ttB-4p-0-900-v1510_14TEV", 2.6673*1000, Background, Category);
		analyze("ttB-4p-900-1600-v1510_14TEV", 0.250469*1000, Background, Category);
		analyze("ttB-4p-1600-2500-v1510_14TEV", 0.0237441*1000, Background, Category);
		analyze("ttB-4p-2500-100000-v1510_14TEV", 0.00208816*1000, Background, Category);
	}
	
	else if(BackgroundSample == "Signal"){
		cout << "Processing only signal events" << endl;
	}
	
	else if(BackgroundSample == "bbgg"){
		// Analyze QCD bbgg background
		analyzeQCD("bbgg", 136, Category);
	}
	
	else if(BackgroundSample == "bbjj"){
		// Analyze QCD bbjj background
		analyzeQCD("bbjj", 6.1e+8, Category);
	}
	
	else if(BackgroundSample == "ccgg"){
		// Analyze QCD ccgg background
		analyzeQCD("ccgg", 1294.2, Category);
	}
	
	else if(BackgroundSample == "ccjg"){
		// Analyze QCD ccjg background
		analyzeQCD("ccjg", 2.5e+6, Category);
	}
	
	else if(BackgroundSample == "ccjj"){
		// Analyze QCD ccjj background
		analyzeQCD("ccjj", 6.5e+8, Category);
	}
	
	else if(BackgroundSample == "jjgg"){
		// Analyze QCD jjgg background
		analyzeQCD("jjgg", 2.2e+4, Category);
	}
	
	else if(BackgroundSample == "jjjj"){
		// Analyze QCD jjjj background
		analyzeQCD("jjjj", 2e+10, Category);
	}
	
	else if(BackgroundSample == "QCD"){
		// Analyze QCD backgrounds
		analyzeQCD("bbgg", 136, Category);
		analyzeQCD("bbjj", 6.1e+8, Category);
		analyzeQCD("ccgg", 1294.2, Category);
		analyzeQCD("ccjg", 2.5e+6, Category);
		analyzeQCD("ccjj", 6.5e+8, Category);
		analyzeQCD("jjgg", 2.2e+4, Category);
		analyzeQCD("jjjj", 2e+10, Category);
	}
	
	else if(BackgroundSample == "All"){
		
		// Analyze B background
		analyze("B-4p-0-1-v1510_14TEV", 200944.68129*1000, Background, Category);

		// Analyse BB background
		analyze("BB-4p-0-300-v1510_14TEV", 249.97710*1000, Background, Category);
		analyze("BB-4p-300-700-v1510_14TEV", 35.23062*1000, Background, Category);
		analyze("BB-4p-700-1300-v1510_14TEV", 4.13743*1000, Background, Category);
		analyze("BB-4p-1300-2100-v1510_14TEV", 0.41702*1000, Background, Category);
		analyze("BB-4p-2100-100000-v1510_14TEV", 0.04770*1000, Background, Category);

		// Analyze BBB background
		analyze("BBB-4p-0-600-v1510_14TEV", 2.57304*1000, Background, Category);
		analyze("BBB-4p-600-1300-v1510_14TEV", 0.14935*1000, Background, Category);
		analyze("BBB-4p-1300-100000-v1510_14TEV", 0.01274*1000, Background, Category);

		// Analyze Bj background
		analyze("Bj-4p-0-300-v1510_14TEV", 34409.92339*1000, Background, Category);
		analyze("Bj-4p-300-600-v1510_14TEV", 2642.85309*1000, Background, Category);
		analyze("Bj-4p-600-1100-v1510_14TEV", 294.12311*1000, Background, Category);
		analyze("Bj-4p-1100-1800-v1510_14TEV", 25.95000*1000, Background, Category);
		analyze("Bj-4p-1800-2700-v1510_14TEV", 2.42111*1000, Background, Category);
		analyze("Bj-4p-2700-3700-v1510_14TEV", 0.22690*1000, Background, Category);
		analyze("Bj-4p-3700-100000-v1510_14TEV", 0.02767*1000, Background, Category);

		// Analyze Bjj-vbf background
		analyze("Bjj-vbf-4p-0-700-v1510_14TEV", 86.45604*1000, Background, Category);
		analyze("Bjj-vbf-4p-700-1400-v1510_14TEV", 4.34869*1000, Background, Category);
		analyze("Bjj-vbf-4p-1400-2300-v1510_14TEV", 0.32465*1000, Background, Category);
		analyze("Bjj-vbf-4p-2300-3400-v1510_14TEV", 0.03032*1000, Background, Category);
		//analyze("Bjj-vbf-4p-3400-100000-v1510_14TEV", 0.00313*1000, Background, Category);

		// Analyze H background
		analyze("H-4p-0-300-v1510_14TEV", 21.55990*1000, Background, Category);
		analyze("H-4p-300-800-v1510_14TEV", 1.11282*1000, Background, Category);
		analyze("H-4p-800-1500-v1510_14TEV", 0.09188*1000, Background, Category);
		analyze("H-4p-1500-100000-v1510_14TEV", 0.01009*1000, Background, Category);

		// Analyze LL background
		analyze("LL-4p-0-100-v1510_14TEV", 1341.36923*1000, Background, Category);
		analyze("LL-4p-100-200-v1510_14TEV", 156.29534*1000, Background, Category);
		analyze("LL-4p-200-500-v1510_14TEV", 42.40132*1000, Background, Category);
		analyze("LL-4p-500-900-v1510_14TEV", 2.84373*1000, Background, Category);
		analyze("LL-4p-900-1400-v1510_14TEV", 0.20914*1000, Background, Category);
		analyze("LL-4p-1400-100000-v1510_14TEV", 0.02891*1000, Background, Category);

		// Analyze LLB background
		analyze("LLB-4p-0-400-v1510_14TEV", 2.97380*1000, Background, Category);
		//analyze("LLB-4p-400-900-v1510_14TEV", 0.22854*1000, Background, Category);
		analyze("LLB-4p-900-100000-v1510_14TEV", 0.02080*1000, Background, Category);

		// Analyze tB background
		analyze("tB-4p-0-500-v1510_14TEV", 63.88923*1000, Background, Category);
		analyze("tB-4p-500-900-v1510_14TEV", 7.12172*1000, Background, Category);
		analyze("tB-4p-900-1500-v1510_14TEV", 0.98030*1000, Background, Category);
		analyze("tB-4p-1500-2200-v1510_14TEV", 0.08391*1000, Background, Category);
		analyze("tB-4p-2200-100000-v1510_14TEV", 0.00953*1000, Background, Category);

		// Analyze tj background
		analyze("tj-4p-0-500-v1510_14TEV", 109.73602*1000, Background, Category);
		analyze("tj-4p-500-1000-v1510_14TEV", 5.99325*1000, Background, Category);
		analyze("tj-4p-1000-1600-v1510_14TEV", 0.37680*1000, Background, Category);
		analyze("tj-4p-1600-2400-v1510_14TEV", 0.03462*1000, Background, Category);
		analyze("tj-4p-2400-100000-v1510_14TEV", 0.00312*1000, Background, Category);

		// Analyse tt background
		analyze("tt-4p-0-600-v1510_14TEV", 530.89358*1000, Background, Category);
		analyze("tt-4p-600-1100-v1510_14TEV", 42.55351*1000, Background, Category);
		analyze("tt-4p-1100-1700-v1510_14TEV", 4.48209*1000, Background, Category);
		analyze("tt-4p-1700-2500-v1510_14TEV", 0.52795*1000, Background, Category);
		analyze("tt-4p-2500-100000-v1510_14TEV", 0.05449*1000, Background, Category);

		// Analyze ttB background
		analyze("ttB-4p-0-900-v1510_14TEV", 2.6673*1000, Background, Category);
		analyze("ttB-4p-900-1600-v1510_14TEV", 0.250469*1000, Background, Category);
		analyze("ttB-4p-1600-2500-v1510_14TEV", 0.0237441*1000, Background, Category);
		analyze("ttB-4p-2500-100000-v1510_14TEV", 0.00208816*1000, Background, Category);
		
		// Analyze QCD backgrounds
		analyzeQCD("bbgg", 136, Category);
		analyzeQCD("bbjj", 6.1e+8, Category);
		analyzeQCD("ccgg", 1294.2, Category);
		analyzeQCD("ccjg", 2.5e+6, Category);
		analyzeQCD("ccjj", 6.5e+8, Category);
		analyzeQCD("jjgg", 2.2e+4, Category);
		analyzeQCD("jjjj", 2e+10, Category);
		
	}
	
	else {
		cout << "Background sample not found" << endl;
		return;
	}
	
	// Save results
	if(BackgroundSample == "Signal") saveResultsS();
	else saveResults();
	
}


/*
 * ANALYSIS
 */

void analyze(TString inputfile, Double_t crossSection, bool SorB, TString category)
{	
	
	// Set up temporal variables
	Double_t tempSelection=0, tempKinCut=0, tempMassCut=0;
	
	Double_t tempErrorSelection=0, tempErrorKinCut=0, tempErrorMassCut=0;
	
	const TString inputFile = "/afs/cern.ch/work/a/ariostas/public/bbgg/" + inputfile + ".root";
	
	inputfile = "Reading " + inputfile + " events... ";
	
	inputfile.Resize(60);
	
	cout << inputfile;
	
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
		
		(SorB ? hPhoton1PtS : hPhoton1PtB)->Fill(photon1->PT, weight);
		(SorB ? hPhoton2PtS : hPhoton2PtB)->Fill(photon2->PT, weight);
		(SorB ? hBJet1PtS : hBJet1PtB)->Fill(bJet1->PT, weight);
		(SorB ? hBJet2PtS : hBJet2PtB)->Fill(bJet2->PT, weight);

		(SorB ? hnLeptonsS : hnLeptonsB)->Fill(nLeptons, weight);
		(SorB ? hnJetsS : hnJetsB)->Fill(nJets, weight);
		
		if(category == "HighPt"){
			if((photon1->PT < 80) || (photon2->PT < 50) || (bJet1->PT < 80) || (bJet2->PT < 50) 
			|| (nLeptons != 0) || (nJets >= 4)) continue;
		}
		else if(category == "LowPt"){
			if((photon1->PT < 40) || (photon2->PT < 25) || (bJet1->PT < 30) || (bJet2->PT < 30) 
			|| ((photon1->PT > 80) && (photon2->PT > 50) && (bJet1->PT > 80) && (bJet2->PT > 50)) || (nLeptons != 0) | (nJets >= 4)) continue;
		}
			
		// Set up four-vectors
		vPhoton1.SetPtEtaPhiE(photon1->PT, photon1->Eta, photon1->Phi, photon1->E);
		vPhoton2.SetPtEtaPhiE(photon2->PT, photon2->Eta, photon2->Phi, photon2->E);
		vBJet1.SetPtEtaPhiM(bJet1->PT, bJet1->Eta, bJet1->Phi, bJet1->Mass);
		vBJet2.SetPtEtaPhiM(bJet2->PT, bJet2->Eta, bJet2->Phi, bJet2->Mass);
	
		photonSystem = vPhoton1 + vPhoton2;
		bJetSystem = vBJet1 + vBJet2;
		
		diHiggsSystem = vPhoton1 + vPhoton2 + vBJet1 + vBJet2;
					
		tempSelection += weight;
		tempErrorSelection++;
		
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
		(SorB ? hdRBJetsS : hdRBJetsB)->Fill(dRBJets, weight);
		
		// Check for kinematic requirements
		if(category == "HighPt"){
			if((dRPhotons < 0.5) || (dRPhotons > 1.25) || (dRBJets < 0.5) || (dRBJets > 1.25) || (dRPhotonBJet < 1.75) || (dRPhotonBJet > 3)) continue;
		}
		else if(category == "LowPt"){
			if((dRPhotons < 1) || (dRPhotons > 1.5) || (dRBJets < 0.75) || (dRBJets > 2) || (dRPhotonBJet < 0.5) || (dRPhotonBJet > 2.25)) continue;
		}
					
		tempKinCut += weight;
		tempErrorKinCut++;
		
		(SorB ? hmPhotonsS : hmPhotonsB)->Fill(photonSystem.M(), weight);
		(SorB ? hmBJetsS : hmBJetsB)->Fill(bJetSystem.M(), weight);
		(SorB ? hmDiHiggsS : hmDiHiggsB)->Fill(diHiggsSystem.M(), weight);
	
		// Check for mass cut requirements
		if(category == "HighPt"){
			if((photonSystem.M() < 120) || (photonSystem.M() > 130) || (bJetSystem.M() < 100) || (bJetSystem.M() > 135) 
			|| (diHiggsSystem.M() < 400) || (diHiggsSystem.M() > 1000)) continue;
		}
		else if(category == "LowPt"){
			if((photonSystem.M() < 122.5) || (photonSystem.M() > 127.5) || (bJetSystem.M() < 80) || (bJetSystem.M() > 130) 
			|| (diHiggsSystem.M() < 350) || (diHiggsSystem.M() > 600)) continue;
		}
		
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
	
	TString out = "";
	out += tempMassCut;
	out.Resize(8);
	
	cout << out << " events passed all cuts" << endl;

}

/* 
 * QCD ANALYSIS
 */
 
void analyzeQCD(TString inputfile, const Double_t crossSection, TString category){
	
	const TString inputFile = "/afs/cern.ch/work/a/ariostas/public/QCD_backgrounds/" + inputfile + ".root";
	
	inputfile = "Reading QCD " + inputfile + " events... ";
	
	inputfile.Resize(60);
	
	cout << inputfile;
	
	// Set up storage variables
	Int_t n_particles=0;
	vector<Int_t> *PID_v=0;
	vector<Double_t> *P_X_v=0;
	vector<Double_t> *P_Y_v=0;
	vector<Double_t> *P_Z_v=0;
	vector<Double_t> *E_v=0;
	vector<Double_t> *M_v=0;
	TLorentzVector vTemp, vPhoton1, vPhoton2, vBJet1, vBJet2, photonSystem, bJetSystem, diHiggsSystem;
	Int_t p1ID=0, p2ID=0, p3ID=0, p4ID=0;
	Double_t dRPhotons, dRBJets, dRPhotonBJet;

	TFile* infile = new TFile(inputFile); assert(infile);
	TTree* intree = (TTree*) infile->Get("Events"); assert(intree);
	
	intree->SetBranchAddress("n_particles",		&n_particles);
	intree->SetBranchAddress("PID",		&PID_v);
	intree->SetBranchAddress("P_X",		&P_X_v);
	intree->SetBranchAddress("P_Y",		&P_Y_v);
	intree->SetBranchAddress("P_Z",		&P_Z_v);
	intree->SetBranchAddress("E",		&E_v);
	intree->SetBranchAddress("M",       &M_v);
	
	Double_t selection=0, kinCut=0, massCut=0;
	Double_t tempErrorSelection=0, tempErrorKinCut=0, tempErrorMassCut=0;

			
	for (Long64_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // Event loop
		intree->GetEntry(iEntry);

		vTemp.SetPxPyPzE(0,0,0,0);
		vPhoton1=vPhoton2=vBJet1=vBJet2=vTemp;
		
		Double_t tempSelection=0, tempKinCut=0, tempMassCut=0, tempWeight=0;
		
		TLorentzVector v1, v2, v3, v4;

		for(Int_t i=0; i<n_particles; i++){
			
			vTemp.SetPxPyPzE(P_X_v->at(i), P_Y_v->at(i), P_Z_v->at(i), E_v->at(i));
			
			if((abs(PID_v->at(i)) >= 1 && abs(PID_v->at(i)) <= 5) || abs(PID_v->at(i)) == 21 || abs(PID_v->at(i)) == 22){
			
				if(vTemp.Pt() > v1.Pt()){
					p4ID = p3ID;
					v4 = v3;
					p3ID = p2ID;
					v3 = v2;
					p2ID = p1ID;
					v2 = v1;
					p1ID = abs(PID_v->at(i));
					v1 = vTemp;
				}
				else if(vTemp.Pt() > v2.Pt()){
					p4ID = p3ID;
					v4 = v3;
					p3ID = p2ID;
					v3 = v2;
					p2ID = abs(PID_v->at(i));
					v2 = vTemp;
				}
				else if(vTemp.Pt() > v3.Pt()){
					p4ID = p3ID;
					v4 = v3;
					p3ID = abs(PID_v->at(i));
					v3 = vTemp;
				}
				else if(vTemp.Pt() > v4.Pt()){
					p4ID = abs(PID_v->at(i));
					v4 = vTemp;
				}
			}
			
		}
		
		for(Int_t x=0; x<6; x++){
			
			if(x==0){
				tempWeight=3000*crossSection/intree->GetEntries();
				vBJet1 = v1;
				tempWeight *= (p1ID == 22 ? 0 : (p1ID == 5 ? 0.6 : (p1ID == 4 ? 0.15 : 0.01)));
				vBJet2 = v2;
				tempWeight *= (p2ID == 22 ? 0 : (p2ID == 5 ? 0.6 : (p2ID == 4 ? 0.15 : 0.01)));
				vPhoton1 = v3;
				tempWeight *= (p3ID == 22 ? 0.8 : (p3ID == 21 ? 5e-4 : 2e-3));
				vPhoton2 = v4;
				tempWeight *= (p4ID == 22 ? 0.8 : (p4ID == 21 ? 5e-4 : 2e-3));
			}
			else if(x==1){
				tempWeight=3000*crossSection/intree->GetEntries();
				vBJet1 = v3;
				tempWeight *= (p3ID == 22 ? 0 : (p3ID == 5 ? 0.6 : (p3ID == 4 ? 0.15 : 0.01)));
				vBJet2 = v4;
				tempWeight *= (p4ID == 22 ? 0 : (p4ID == 5 ? 0.6 : (p4ID == 4 ? 0.15 : 0.01)));
				vPhoton1 = v1;
				tempWeight *= (p1ID == 22 ? 0.8 : (p1ID == 21 ? 5e-4 : 2e-3));
				vPhoton2 = v2;
				tempWeight *= (p2ID == 22 ? 0.8 : (p2ID == 21 ? 5e-4 : 2e-3));
			}
			else if(x==2){
				tempWeight=3000*crossSection/intree->GetEntries();
				vBJet1 = v1;
				tempWeight *= (p1ID == 22 ? 0 : (p1ID == 5 ? 0.6 : (p1ID == 4 ? 0.15 : 0.01)));
				vBJet2 = v3;
				tempWeight *= (p3ID == 22 ? 0 : (p3ID == 5 ? 0.6 : (p3ID == 4 ? 0.15 : 0.01)));
				vPhoton1 = v2;
				tempWeight *= (p2ID == 22 ? 0.8 : (p2ID == 21 ? 5e-4 : 2e-3));
				vPhoton2 = v4;
				tempWeight *= (p4ID == 22 ? 0.8 : (p4ID == 21 ? 5e-4 : 2e-3));
			}
			else if(x==3){
				tempWeight=3000*crossSection/intree->GetEntries();
				vBJet1 = v2;
				tempWeight *= (p2ID == 22 ? 0 : (p2ID == 5 ? 0.6 : (p2ID == 4 ? 0.15 : 0.01)));
				vBJet2 = v4;
				tempWeight *= (p4ID == 22 ? 0 : (p4ID == 5 ? 0.6 : (p4ID == 4 ? 0.15 : 0.01)));
				vPhoton1 = v1;
				tempWeight *= (p1ID == 22 ? 0.8 : (p1ID == 21 ? 5e-4 : 2e-3));
				vPhoton2 = v3;
				tempWeight *= (p3ID == 22 ? 0.8 : (p3ID == 21 ? 5e-4 : 2e-3));
			}
			else if(x==4){
				tempWeight=3000*crossSection/intree->GetEntries();
				vBJet1 = v1;
				tempWeight *= (p1ID == 22 ? 0 : (p1ID == 5 ? 0.6 : (p1ID == 4 ? 0.15 : 0.01)));
				vBJet2 = v4;
				tempWeight *= (p4ID == 22 ? 0 : (p4ID == 5 ? 0.6 : (p4ID == 4 ? 0.15 : 0.01)));
				vPhoton1 = v2;
				tempWeight *= (p2ID == 22 ? 0.8 : (p2ID == 21 ? 5e-4 : 2e-3));
				vPhoton2 = v3;
				tempWeight *= (p3ID == 22 ? 0.8 : (p3ID == 21 ? 5e-4 : 2e-3));
			}
			else if(x==5){
				tempWeight=3000*crossSection/intree->GetEntries();
				vBJet1 = v2;
				tempWeight *= (p2ID == 22 ? 0 : (p2ID == 5 ? 0.6 : (p2ID == 4 ? 0.15 : 0.01)));
				vBJet2 = v3;
				tempWeight *= (p3ID == 22 ? 0 : (p3ID == 5 ? 0.6 : (p3ID == 4 ? 0.15 : 0.01)));
				vPhoton1 = v1;
				tempWeight *= (p1ID == 22 ? 0.8 : (p1ID == 21 ? 5e-4 : 2e-3));
				vPhoton2 = v4;
				tempWeight *= (p4ID == 22 ? 0.8 : (p4ID == 21 ? 5e-4 : 2e-3));
			}
		
			//Check if there are selected objects
			if((vPhoton1.Pt() < 25) || (vPhoton2.Pt() < 25) || (vBJet1.Pt() < 30) || (vBJet2.Pt() < 30)
				|| (vPhoton1.Eta() > 2.5) || (vPhoton2.Eta() > 2.5) || (vBJet1.Eta() > 2.4) || (vBJet2.Eta() > 2.4)) continue;
			
			hPhoton1PtB->Fill(vPhoton1.Pt(), tempWeight);
			hPhoton2PtB->Fill(vPhoton2.Pt(), tempWeight);
			hBJet1PtB->Fill(vBJet1.Pt(), tempWeight);
			hBJet2PtB->Fill(vBJet2.Pt(), tempWeight);
		
			//Check if objects satisfy selection requirements
			if(category == "HighPt"){
				if((vPhoton1.Pt() < 80) || (vPhoton2.Pt() < 50) || (vBJet1.Pt() < 80) || (vBJet2.Pt() < 50)) continue;
			}
			else if(category == "LowPt"){
				if((vPhoton1.Pt() < 40) || ((vPhoton1.Pt() > 80) && (vPhoton2.Pt() > 50) && (vBJet1.Pt() > 80) && (vBJet2.Pt() > 50))) continue;
			}
		
			photonSystem = vPhoton1 + vPhoton2;
			bJetSystem = vBJet1 + vBJet2;
			diHiggsSystem = vPhoton1 + vPhoton2 + vBJet1 + vBJet2;
		
			if(tempWeight > tempSelection) tempSelection = tempWeight;
			
			//Calculate required deltaRs
			dRPhotons = deltaR(vPhoton1.Eta(), vPhoton2.Eta(), vPhoton1.Phi(), vPhoton2.Phi());
			dRBJets = deltaR(vBJet1.Eta(), vBJet2.Eta(), vBJet1.Phi(), vBJet2.Phi());
			Float_t dR1 = deltaR(vBJet1.Eta(), vPhoton1.Eta(), vBJet1.Phi(), vPhoton1.Phi());
			Float_t dR2 = deltaR(vBJet1.Eta(), vPhoton2.Eta(), vBJet1.Phi(), vPhoton2.Phi());
			Float_t dR3 = deltaR(vBJet2.Eta(), vPhoton1.Eta(), vBJet2.Phi(), vPhoton1.Phi());
			Float_t dR4 = deltaR(vBJet2.Eta(), vPhoton2.Eta(), vBJet2.Phi(), vPhoton2.Phi());
		
			dRPhotonBJet = dR1;
			if (dRPhotonBJet > dR2) dRPhotonBJet = dR2;
			if (dRPhotonBJet > dR3) dRPhotonBJet = dR3;
			if (dRPhotonBJet > dR4) dRPhotonBJet = dR4;
		
			hdRPhotonsB->Fill(dRPhotons, tempWeight);
			hdRPhotonBJetB->Fill(dRPhotonBJet, tempWeight);
			hdRBJetsB->Fill(dRBJets, tempWeight);
		
			//Check for kinematic requirements
			if(category == "HighPt"){
				if((dRPhotons < 0.5) || (dRPhotons > 1.25) || (dRBJets < 0.5) || (dRBJets > 1.25) || (dRPhotonBJet < 1.75) || (dRPhotonBJet > 3)) continue;
			}
			else if(category == "LowPt"){
				if((dRPhotons < 1) || (dRPhotons > 1.5) || (dRBJets < 0.75) || (dRBJets > 2) || (dRPhotonBJet < 0.5) || (dRPhotonBJet > 2.25)) continue;
			}
		
			if(tempWeight > tempKinCut) tempKinCut = tempWeight;
			
			hmPhotonsB->Fill(photonSystem.M(), tempWeight);
			hmBJetsB->Fill(bJetSystem.M(), tempWeight);
			hmDiHiggsB->Fill(diHiggsSystem.M(), tempWeight);
		
			//Check for mass cut requirements
			if(category == "HighPt"){
				if((photonSystem.M() < 120) || (photonSystem.M() > 130) || (bJetSystem.M() < 100) || (bJetSystem.M() > 135) 
				|| (diHiggsSystem.M() < 400) || (diHiggsSystem.M() > 1000)) continue;
			}
			else if(category == "LowPt"){
				if((photonSystem.M() < 122.5) || (photonSystem.M() > 127.5) || (bJetSystem.M() < 80) || (bJetSystem.M() > 130) 
				|| (diHiggsSystem.M() < 350) || (diHiggsSystem.M() > 600)) continue;
			}
		
			if(tempWeight > tempMassCut) tempMassCut = tempWeight;
		
		}
		
		selection+=tempSelection;
		kinCut+=tempKinCut;
		massCut+=tempMassCut;
	
		if(tempSelection > 0) tempErrorSelection++;
		if(tempKinCut > 0) tempErrorKinCut++;
		if(tempMassCut > 0) tempErrorMassCut++;
		
	}
	
	if(tempErrorSelection != 0) errorSelectionBackground += sqrtf(tempErrorSelection)*selection/tempErrorSelection;
	if(tempErrorKinCut != 0) errorKinCutBackground += sqrtf(tempErrorKinCut)*kinCut/tempErrorKinCut;
	if(tempErrorMassCut != 0) errorMassCutBackground += sqrt(tempErrorMassCut)*massCut/tempErrorMassCut;
	
	totalBackground += 3000*crossSection;
	selectionBackground += selection;
	kinCutBackground += kinCut;
	massCutBackground += massCut;
	
	TString out = "";
	out += massCut;
	out.Resize(8);
	
	cout << out << " events passed all cuts" << endl;
	
}

/*
 * SAVE RESULTS
 */
 
void saveResults(){
	
	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);
	
	gStyle->SetOptStat(kFALSE);
	
	histogram(hnLeptonsS, hnLeptonsB, c1, "Number of Leptons", "Count", "./Histograms/histogramNL_" + BackgroundSample + "_" + Category + ".jpg");
	histogram(hnJetsS, hnJetsB, c1, "Number of Jets with |Eta|<2.5", "Count", "./Histograms/histogramNJ_" + BackgroundSample + "_" + Category + ".jpg");
	histogram(hdRPhotonsS, hdRPhotonsB, c1, "delta R of Photons", "Count", "./Histograms/histogramDRGG_" + BackgroundSample + "_" + Category + ".jpg");
	histogram(hdRBJetsS, hdRBJetsB, c1, "delta R of BJets", "Count", "./Histograms/histogramdRBB_" + BackgroundSample + "_" + Category + ".jpg");
	histogram(hdRPhotonBJetS, hdRPhotonBJetB, c1, "min delta R of Photon and BJet", "Count", "./Histograms/histogramdRPGB_" + BackgroundSample + "_" + Category + ".jpg");
	histogram(hmPhotonsS, hmPhotonsB, c1, "Reconstructed mass from Photons", "Count", "./Histograms/histogramMP_" + BackgroundSample + "_" + Category + ".jpg");
	histogram(hmBJetsS, hmBJetsB, c1, "Reconstructed mass from BJets", "Count", "./Histograms/histogramMBJ_" + BackgroundSample + "_" + Category + ".jpg");
	histogram(hmDiHiggsS, hmDiHiggsB, c1, "Reconstructed mass from diHiggs", "Count", "./Histograms/histogramMDH_" + BackgroundSample + "_" + Category + ".jpg");
	
	histogram(hPhoton1PtS, hPhoton1PtB, c1, "Pt of leading photon", "Count", "./Histograms/histogramPhoton1Pt_" + BackgroundSample + "_" + Category + ".jpg");
	histogram(hPhoton2PtS, hPhoton2PtB, c1, "Pt of second leading photon", "Count", "./Histograms/histogramPhoton2Pt_" + BackgroundSample + "_" + Category + ".jpg");
	histogram(hBJet1PtS, hBJet1PtB, c1, "Pt of leading bjet", "Count", "./Histograms/histogramBJet1Pt_" + BackgroundSample + "_" + Category + ".jpg");
	histogram(hBJet2PtS, hBJet2PtB, c1, "Pt of second leading bjet", "Count", "./Histograms/histogramBJet2Pt_" + BackgroundSample + "_" + Category + ".jpg");
	
	cout << "\nSignal\n" << endl;
	cout << "Total Events " << totalSignal << endl;
	cout << "Events after selection: " << selectionSignal << " +- " << errorSelectionSignal << endl;
	cout << "Events after kinematic cut: " << kinCutSignal << " +- " << errorKinCutSignal << endl;
	cout << "Events after mass cut: " << massCutSignal << " +- " << errorMassCutSignal << endl;
	
	cout << "\n\n" << BackgroundSample << " Background\n" << endl;
	cout << "Total Events " << totalBackground << endl;
	cout << "Events after selection: " << selectionBackground << " +- " << errorSelectionBackground << endl;
	cout << "Events after kinematic cut: " << kinCutBackground << " +- " << errorKinCutBackground << endl;
	cout << "Events after mass cut: " << massCutBackground << " +- " << errorMassCutBackground << endl;
	
}

void saveResultsS(){
	
	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);
	
	gStyle->SetOptStat(kFALSE);
	
	histogram(hnLeptonsS, c1, "Number of Leptons", "Count", "./Histograms/histogramnL_Signal" + "_" + Category + ".jpg");
	histogram(hnJetsS, c1, "Number of Jets with |Eta|<2.5", "Count", "./Histograms/histogramnJ_Signal" + "_" + Category + ".jpg");
	histogram(hdRPhotonsS, c1, "delta R of Photons", "Count", "./Histograms/histogramdRP_Signal" + "_" + Category + ".jpg");
	histogram(hdRBJetsS, hdRBJetsB, c1, "delta R of BJets", "Count", "./Histograms/histogramdRBJ_Signal" + "_" + Category + ".jpg");
	histogram(hdRPhotonBJetS, c1, "min delta R of Photon and BJet", "Count", "./Histograms/histogramdRPBJ_Signal" + "_" + Category + ".jpg");
	histogram(hmPhotonsS, c1, "Reconstructed mass from Photons", "Count", "./Histograms/histogramdmP_Signal" + "_" + Category + ".jpg");
	histogram(hmBJetsS, c1, "Reconstructed mass from BJets", "Count", "./Histograms/histogramdmBJ_Signal" + "_" + Category + ".jpg");
	
	histogram(hPhoton1PtS, c1, "Pt of leading photon", "Count", "./Histograms/histogramPhoton1Pt_Signal_" + Category + ".jpg");
	histogram(hPhoton2PtS, c1, "Pt of second leading photon", "Count", "./Histograms/histogramPhoton2Pt_Signal_" + Category + ".jpg");
	histogram(hBJet1PtS, c1, "Pt of leading bjet", "Count", "./Histograms/histogramBJet1Pt_Signal_" + Category + ".jpg");
	histogram(hBJet2PtS, c1, "Pt of second leading bjet", "Count", "./Histograms/histogramBJet2Pt_Signal_" + Category + ".jpg");
	
	cout << "\nSignal\n" << endl;
	cout << "Total Events " << totalSignal << endl;
	cout << "Events after selection: " << selectionSignal << " +- " << errorSelectionSignal << endl;
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
