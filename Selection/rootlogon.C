{
	TString addedLibs(gSystem->GetLibraries());
	if(!addedLibs.Contains("setRootEnv_C.so")) {
      gROOT->Macro("/afs/cern.ch/user/a/ariostas/bbgg/Selection/setRootEnv.C+");
    }
	gSystem->Load("/afs/cern.ch/user/a/ariostas/Delphes/libDelphes.so");
	loadLibraries("libTauAnalysisSVFitHelper.so");
    loadLibraries("libTauAnalysisCandidateTools.so");
    loadLibraries("libJECJECHelper.so");
	gROOT->ProcessLine(".include /afs/cern.ch/user/a/ariostas/Delphes");
	gROOT->ProcessLine(".include /afs/cern.ch/user/a/ariostas/Delphes/external");

}

