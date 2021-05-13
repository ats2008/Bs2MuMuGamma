#include "Bs2MuMuGamma/Ntuplizer/interface/TreeContent.h"
#include <iostream>

TreeContent::TreeContent ()
{
  ClearScalars();

  // ### Trigger ###
  TrigTable     = nullptr;
  TrigPrescales = nullptr;
  L1Table       = nullptr;
  L1Prescales   = nullptr;
 // hltObjs       = nullptr;
  
    nMuons=0;
    muon_pt		=nullptr	;
    muon_eta		=nullptr	;
    muon_phi		=nullptr	;
    mum_dz		=nullptr	;
    muon_dxy		=nullptr	;
    mum_dz_error	=nullptr	;
    muon_dxy_error	=nullptr	;
    muon_vx		=nullptr	;
    muon_vy		=nullptr	;
    muon_vz		=nullptr	;
    muon_vertexChi2	=nullptr	;
    muon_vertexNDoF	=nullptr	;
    muon_charge	        =nullptr	;
    muon_isGlobalMuon	=nullptr	;
    muon_isTrackerMuon	=nullptr	;
    muon_StandAloneMuon =nullptr	;
    muon_isCaloMuon	=nullptr	;
    muon_isPFMuon	=nullptr	;

    muon_selector         	= nullptr; 
    muon_isIsolationValid 	= nullptr;
    muon_isolationR03      	= nullptr;
    muon_isolationR05      	= nullptr;
    muon_isPFIsolationValid 	= nullptr;
    muon_pfIsolationR03      	= nullptr;
    muon_pfIsolationR04      	= nullptr;

  // ### mu- ###
  //mumHighPurity    = nullptr;
  //mumCL            = nullptr;
  //mumNormChi2      = nullptr;
  //mumPx            = nullptr;
  //mumPy            = nullptr;
  //mumPz            = nullptr;
  //mumDCAVtx        = nullptr;
  //mumDCAVtxE       = nullptr;
  //mumDCABS         = nullptr;
  //mumDCABSE        = nullptr;
  //mumKinkChi2      = nullptr;
  //mumFracHits      = nullptr;
  //mumdxyBS         = nullptr;
  //mumdzBS          = nullptr;
  //mumMinIP2D       = nullptr;
  //mumMinIP2DE      = nullptr;
  //mumMinIP         = nullptr;
  //mumMinIPS        = nullptr;
  //mumDeltaRwithMC  = nullptr;
  //mumCat           = nullptr;
  //mumNPixHits      = nullptr;
  //mumNPixLayers    = nullptr;
  //mumNTrkHits      = nullptr;
  //mumNTrkLayers    = nullptr;
  //mumNMuonHits     = nullptr;
  //mumNMatchStation = nullptr;
  //mumIso           = nullptr;
  //mumIsoPt         = nullptr;
  //mumIsodR         = nullptr;

  }

void TreeContent::Init ()
{
  // ### Trigger ###
  TrigTable     = new std::vector<std::string>;
  TrigPrescales = new std::vector<int>;
  L1Table       = new std::vector<std::string>;
  L1Prescales   = new std::vector<int>;
  //hltObjs       = new std::vector<miniHLTObj>;
  
    nMuons= new int;
    muon_pt		= new std::vector<double>	;
    muon_eta		= new std::vector<double>	;
    muon_phi		= new std::vector<double>	;
    mum_dz		= new std::vector<double>	;
    muon_dxy		= new std::vector<double>	;
    mum_dz_error	= new std::vector<double>	;
    muon_dxy_error	= new std::vector<double>	;
    muon_vx		= new std::vector<double>	;
    muon_vy		= new std::vector<double>	;
    muon_vz		= new std::vector<double>	;
    muon_vertexChi2	= new std::vector<double>	;
    muon_vertexNDoF	= new std::vector<double>	;
    muon_charge	        = new std::vector<int>	;
    muon_isGlobalMuon	= new std::vector<bool>	;
    muon_isTrackerMuon	= new std::vector<bool>	;
    muon_StandAloneMuon = new std::vector<bool>	;
    muon_isCaloMuon	= new std::vector<bool>	;
    muon_isPFMuon	= new std::vector<bool>	;

    muon_selector         =  new std::vector<uint64_t>; 
 
    muon_isIsolationValid =  new std::vector<bool>;;
    muon_isolationR03      =  new std::vector<reco::MuonIsolation>;
    muon_isolationR05      =  new std::vector<reco::MuonIsolation>;
    muon_isPFIsolationValid =  new std::vector<bool>;;
    muon_pfIsolationR03      =  new std::vector<reco::MuonPFIsolation>;
    muon_pfIsolationR04      =  new std::vector<reco::MuonPFIsolation>;

}

TreeContent::~TreeContent ()
{
   delete  TrigTable     	;	
   delete  TrigPrescales 	;
   delete  L1Table       	;
   delete  L1Prescales   	;
   //delete  hltObjs       	;
  
   delete  nMuons		;
   delete  muon_pt		;
   delete  muon_eta		;
   delete  muon_phi		;
   delete  mum_dz		;
   delete  muon_dxy		;
   delete  mum_dz_error	;
   delete  muon_dxy_error	;
   delete  muon_vx		;
   delete  muon_vy		;
   delete  muon_vz		;
   delete  muon_vertexChi2	;
   delete  muon_vertexNDoF	;
   delete  muon_charge	        ;
   delete  muon_isGlobalMuon	;
   delete  muon_isTrackerMuon	;
   delete  muon_StandAloneMuon ;
   delete  muon_isCaloMuon	;
   delete  muon_isPFMuon	;

   delete  muon_selector       ;
 
   delete  muon_isIsolationValid;
   delete  muon_isolationR03   ;
   delete  muon_isolationR05   ;
   delete  muon_isPFIsolationValid;
   delete  muon_pfIsolationR03 ;
   delete  muon_pfIsolationR04 ;


}

void TreeContent::ClearScalars ()
{
  ClearScalarsMonteCarlo();
}

void TreeContent::ClearScalarsMonteCarlo ()
{

}

void TreeContent::ClearVectors ()
{
  // ### Trigger ###
  TrigTable->clear();
  TrigPrescales->clear();
  L1Table->clear();
  L1Prescales->clear();
   //hltObjs -> clear();
  
   *nMuons=0;
   muon_pt		->clear();
   muon_eta		->clear();
   muon_phi		->clear();
   mum_dz		->clear();
   muon_dxy		->clear();
   mum_dz_error	->clear();
   muon_dxy_error	->clear();
   muon_vx		->clear();
   muon_vy		->clear();
   muon_vz		->clear();
   muon_vertexChi2	->clear();
   muon_vertexNDoF	->clear();
   muon_charge	        ->clear();
   muon_isGlobalMuon	->clear();
   muon_isTrackerMuon	->clear();
   muon_StandAloneMuon ->clear();
   muon_isCaloMuon	->clear();
   muon_isPFMuon	->clear();

   muon_selector       ->clear();
 
   muon_isIsolationValid->clear();
   muon_isolationR03   ->clear();
   muon_isolationR05   ->clear();
   muon_isPFIsolationValid->clear();
   muon_pfIsolationR03 ->clear();
   muon_pfIsolationR04 ->clear();


  ClearVectorsMonteCarlo();
}

void TreeContent::ClearVectorsMonteCarlo ()
{
  // ### Matching Between Reconstructed and Generated ###
}

void TreeContent::ClearNTuple ()
{
  ClearScalars();
  ClearVectors();
}

void TreeContent::ClearMonteCarlo ()
{
  ClearScalarsMonteCarlo();
  ClearVectorsMonteCarlo();
}

void TreeContent::MakeTreeBranches (TTree* theTree)
{
  // ########################################################
  // # Run Number, event number, #reco vtx and event weight #
  // ########################################################
  // theTree->Branch("runN",            &runN,            "runN/i");
  // theTree->Branch("eventN",          &eventN,          "eventN/i");
  // theTree->Branch("recoVtxN",        &recoVtxN,        "recoVtxN/i");
  // theTree->Branch("evWeight",        &evWeight,        "evWeight/D");
  // theTree->Branch("evWeightE2",      &evWeightE2,      "evWeightE2/D");
  // theTree->Branch("numEventsTried",  &numEventsTried,  "numEventsTried/i");
  // theTree->Branch("numEventsPassed", &numEventsPassed, "numEventsPassed/i");

  // ### Trigger ###
  theTree->Branch("TrigTable",     &TrigTable);
  theTree->Branch("TrigPrescales", &TrigPrescales);
  theTree->Branch("L1Table",       &L1Table);
  theTree->Branch("L1Prescales",   &L1Prescales);
  //theTree->Branch("hltObjs"    ,   &hltObjs);

   theTree->Branch("nMuons"			,  nMuons			);
   theTree->Branch("muon_pt"			,  &muon_pt			);
   theTree->Branch("muon_eta"		   	,  &muon_eta			);
   theTree->Branch("muon_phi"		   	,  &muon_phi			);
   theTree->Branch("mum_dz"		   	,  &mum_dz			);
   theTree->Branch("muon_dxy"		   	,  &muon_dxy			);
   theTree->Branch("mum_dz_error"	   	,  &mum_dz_error   		);
   theTree->Branch("muon_dxy_error"	   	,  &muon_dxy_error		);
   theTree->Branch("muon_vx"		   	,  &muon_vx			);
   theTree->Branch("muon_vy"		   	,  &muon_vy			);
   theTree->Branch("muon_vz"		   	,  &muon_vz			);
   theTree->Branch("muon_vertexChi2"	   	,  &muon_vertexChi2		);
   theTree->Branch("muon_vertexNDoF"	   	,  &muon_vertexNDoF		);
   theTree->Branch("muon_charge"	        ,  &muon_charge	        	);
   theTree->Branch("muon_isGlobalMuon"	   	,  &muon_isGlobalMuon		);
   theTree->Branch("muon_isTrackerMuon"	   	,  &muon_isTrackerMuon		);
   theTree->Branch("muon_StandAloneMuon"    	,  &muon_StandAloneMuon 	);
   theTree->Branch("muon_isCaloMuon"	   	,  &muon_isCaloMuon		);
   theTree->Branch("muon_isPFMuon"	   	,  &muon_isPFMuon		);
   theTree->Branch("muon_selector"          	,  &muon_selector         	); 
   theTree->Branch("muon_isIsolationValid"  	,  &muon_isIsolationValid 	);
   //theTree->Branch("muon_isolationR03"      	,  &muon_isolationR03      	);
   //theTree->Branch("muon_isolationR05"      	,  &muon_isolationR05      	);
   theTree->Branch("muon_isPFIsolationValid"	,  &muon_isPFIsolationValid 	);
   //theTree->Branch("muon_pfIsolationR03"    	,  &muon_pfIsolationR03      	);
   //theTree->Branch("muon_pfIsolationR04"    	,  &muon_pfIsolationR04      	);


}




