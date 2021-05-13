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

    muon_selector         = nullptr; 
 
    muon_isIsolationValid = nullptr;;
    muon_isolationR03      = nullptr;
    muon_isolationR05      = nullptr;
    muon_isPFIsolationValid = nullptr;;
    muon_pfIsolationR03      = nullptr;
    muon_pfIsolationR04      = nullptr;

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

}

TreeContent::~TreeContent ()
{

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
  theTree->Branch("runN",            &runN,            "runN/i");
  theTree->Branch("eventN",          &eventN,          "eventN/i");
  theTree->Branch("recoVtxN",        &recoVtxN,        "recoVtxN/i");
  theTree->Branch("evWeight",        &evWeight,        "evWeight/D");
  theTree->Branch("evWeightE2",      &evWeightE2,      "evWeightE2/D");
  theTree->Branch("numEventsTried",  &numEventsTried,  "numEventsTried/i");
  theTree->Branch("numEventsPassed", &numEventsPassed, "numEventsPassed/i");

  // ### Trigger ###
  theTree->Branch("TrigTable",     &TrigTable);
  theTree->Branch("TrigPrescales", &TrigPrescales);
  theTree->Branch("L1Table",       &L1Table);
  theTree->Branch("L1Prescales",   &L1Prescales);
  //theTree->Branch("hltObjs"    ,   &hltObjs);


}




