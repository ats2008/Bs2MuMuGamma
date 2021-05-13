// -*- C++ -*-
//
// Package:    Bs2MuMuGamma/NTuplizer
// Class:      NTuplizer
//
/**\class NTuplizer NTuplizer.cc Bs2MuMuGamma/NTuplizer/plugins/NTuplizer.cc

 Description: Takes in AOD and makes NTuples of Muon/DiMuons/Photons

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed, 12 May 2021 18:51:51 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"



// user include files
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"



#include "Bs2MuMuGamma/Ntuplizer/interface/TreeContent.h"
#include "Bs2MuMuGamma/Ntuplizer/interface/Utils.h"




//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class NTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit NTuplizer(const edm::ParameterSet&);
      ~NTuplizer();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      TTree* theTree;
      TreeContent* NTuple;
      Utils* Utility;
      //l1t::L1TGlobalUtil *fGtUtil;

      // selection cuts;
      double pTMinMuons;

      edm::EDGetTokenT<reco::GenParticleCollection>         prunedGenToken_;
      
      //edm::EDGetTokenT<edm::TriggerResults>                    triggerBits_;
      //edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      //edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescales_;
      //std::vector<std::string> trigTable_;
      //std::vector<std::string> l1Table_;
      
      edm::EDGetTokenT<reco::BeamSpot>                 beamSpotToken_;
      edm::EDGetTokenT<reco::VertexCollection>         vtxToken_;
      edm::EDGetTokenT<std::vector<reco::Muon>>         muonToken_;



     // Dimuon Reco vars
         TrajectoryStateClosestToPoint theDCAXBS;
	 ClosestApproachInRPhi ClosestApp;
	 GlobalPoint XingPoint;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
NTuplizer::NTuplizer(const edm::ParameterSet& iConfig)
 :
  muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))

{
   //now do what ever initialization is needed
 
    beamSpotToken_  =consumes<reco::BeamSpot>                    (iConfig.getParameter<edm::InputTag>("beamSpot"));
    vtxToken_       =consumes<reco::VertexCollection>            (iConfig.getParameter<edm::InputTag>("vertices"));
    //trackToken_     =consumes<reco::TrackCollection>             (iConfig.getParameter<edm::InputTag>("tracks"));
    //triggerBits_     =consumes<edm::TriggerResults>                    (iConfig.getParameter<edm::InputTag>("bits"));
    //triggerPrescales_ =consumes<pat::PackedTriggerPrescales>            (iConfig.getParameter<edm::InputTag>("prescales"));
    //trigTable_        =iConfig.getParameter<std::vector<std::string> >("TriggerNames");   
    //triggerObjects_   (consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("objects"));
    
    //l1Table_           = iConfig.getParameter<std::vector<std::string> >("L1Names");   
    //mcGenToken_        = consumes<reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>(""));

    pTMinMuons = 3.50;
    
    NTuple = new TreeContent;
    NTuple->Init();

}


NTuplizer::~NTuplizer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  
  //delete NTuple;
  //delete Utility;


}


//
// member functions
//

// ------------ method called for each event  ------------
void
NTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
 
     // Get magnetic field
    //edm::ESHandle<MagneticField> bFieldHandle;
    //iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);      
    // Get BeamSpot
    edm::Handle<reco::BeamSpot> beamSpotH;
    iEvent.getByToken(beamSpotToken_, beamSpotH);
    //reco::BeamSpot beamSpot = *beamSpotH;

    edm::Handle<std::vector<reco::Muon>> muons;
    iEvent.getByToken(muonToken_, muons);
 
    //  Muon Ntuplizing
    //  TODO : Add details to closest PV
    //         Add details with BS
    //
    int x=0;
    reco::TrackRef aMuonInnerTrack;
    for ( uint32_t i=0;i<muons->size();i++) 
    {
    	 auto &aMuon=muons->at(i);
         
         if(aMuon.pt()  < pTMinMuons) continue;
	 aMuonInnerTrack= aMuon.innerTrack();
         
	 
	 NTuple->muon_pt		      ->push_back(aMuon.pt());
         NTuple->muon_eta		      ->push_back(aMuon.eta());
         NTuple->muon_phi		      ->push_back(aMuon.phi());
	 if(not aMuonInnerTrack.isNull())
	 {
         	NTuple->mum_dz		              ->push_back(aMuonInnerTrack->dz());
         	NTuple->muon_dxy		      ->push_back(aMuonInnerTrack->dxy());
         	NTuple->mum_dz_error	              ->push_back(aMuonInnerTrack->dzError() );
         	NTuple->muon_dxy_error	              ->push_back(aMuonInnerTrack->dxyError());
         }
	 else
	 {
	         NTuple->mum_dz		              ->push_back(1e5);
	         NTuple->muon_dxy		      ->push_back(1e5);
	         NTuple->mum_dz_error	              ->push_back(1e9);
	         NTuple->muon_dxy_error	              ->push_back(1e9);
	 }

	 NTuple->muon_vx		      ->push_back(aMuon.vx() );
         NTuple->muon_vy		      ->push_back(aMuon.vy());
         NTuple->muon_vz		      ->push_back(aMuon.vz());
         NTuple->muon_vertexChi2	      ->push_back(aMuon.vertexChi2());
         NTuple->muon_vertexNDoF	      ->push_back(aMuon.vertexNdof());
         NTuple->muon_charge	              ->push_back(aMuon.charge());
         NTuple->muon_isGlobalMuon	      ->push_back(aMuon.isGlobalMuon());
         NTuple->muon_isTrackerMuon	      ->push_back(aMuon.isTrackerMuon());
         NTuple->muon_StandAloneMuon          ->push_back(aMuon.isStandAloneMuon());
         NTuple->muon_isCaloMuon	      ->push_back(aMuon.isCaloMuon());
         NTuple->muon_isPFMuon	              ->push_back(aMuon.isPFMuon());

         NTuple->muon_selector            ->push_back(aMuon.selectors()); 
         NTuple->muon_isIsolationValid    ->push_back(aMuon.isIsolationValid());
         NTuple->muon_isolationR03         ->push_back(aMuon.isolationR03());
         NTuple->muon_isolationR05         ->push_back(aMuon.isolationR05());
         NTuple->muon_isPFIsolationValid  ->push_back(aMuon.isIsolationValid());
         NTuple->muon_pfIsolationR03       ->push_back(aMuon.pfIsolationR03());
         NTuple->muon_pfIsolationR04       ->push_back(aMuon.pfIsolationR04());
	 *(NTuple->nMuons)+=1; 
	 /*
	 if ((aMuonInnerTrack.isNull() == true) || (aMuonInnerTrack->charge() != -1)) continue;

	 const reco::TransientTrack muTrackpTT( aMuonInnerTrack, &(*bFieldHandle));
  	 if (!muTrackmTT.isValid()) continue;
	
	 // # Compute mu- DCA to BeamSpot #
	 theDCAXBS = muTrackmTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
	 if (theDCAXBS.isValid() == false)
	 {
	      if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu-" << std::endl;
	            continue;
	 }

        double DCAmumBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
        double DCAmumBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
        if (fabs(DCAmumBS) > DCAMUBS)
        {
          if (printMsg) std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu- : " << DCAmumBS << std::endl;
          continue;
        }
    	for ( uint32_t j=0;j<i;j++) 
    	{
    	 	auto &bMuon=muons->at(j);
               
  	       if(bMuon.pt()  < pTMinMuons) continue;
	       
	       bMuonInnerTrack = bMuon.innerTrack();
	       
               if ((bMuonInnerTrack.isNull() == true) || (bMuonInnerTrack->charge() != 1)) continue;
	       
		const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle));
                if (!muTrackpTT.isValid()) continue;
	       
	       // # Compute mu+ DCA to BeamSpot #
               theDCAXBS = muTrackpTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
               if (theDCAXBS.isValid() == false)
               {
                 if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu+" << std::endl;
                 continue;
               }
               double DCAmupBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
               double DCAmupBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
               if (fabs(DCAmupBS) > DCAMUBS)
               {
                 if (printMsg) std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu+: " << DCAmupBS << std::endl;
                 continue;
               }
	     
	    // # Check goodness of muons closest approach #
               ClosestApp.calculate(muTrackpTT.initialFreeState(),muTrackmTT.initialFreeState());
	       XingPoint = ClosestApp.crossingPoint();
            if (  (sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) > TRKMAXR) || (fabs(XingPoint.z()) > TRKMAXZ) )
             {
              if (printMsg) std::cout << __LINE__ << " : continue --> closest approach crossing point outside the tracker volume" << std::endl;
              continue ;
             }
	    double mumuDCA = ClosestApp.distance();
            if (mumuDCA > DCAMUMU)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad 3D-DCA of mu+(-) with respect to mu-(+): " << mumuDCA << std::endl;
              continue;
            }
            
	    // # Cut on the dimuon inviariant mass and pT #
            jpsi_lv.SetPxPyPzE( muTrackmTT.track().px() + muTrackpTT.track().px(), 
                                muTrackmTT.track().py() + muTrackpTT.track().py(),
                                muTrackmTT.track().pz() + muTrackpTT.track().pz(),
                                sqrt( pow(muTrackmTT.track().p(),2) + pow(Utility->muonMass,2) ) + sqrt( pow(muTrackpTT.track().p(),2) + pow(Utility->muonMass,2) )
                              );

	    if ((jpsi_lv.Pt() < MINMUMUPT)  || (jpsi_lv.M() < MINMUMUINVMASS) || (jpsi_lv.M() > MAXMUMUINVMASS))
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << jpsi_lv.Pt() << "\tinv. mass: " << jpsi_lv.M() << std::endl;
              continue;
            }

  	    chi2 = 0.;
            ndf  = 0.;
            // ####################################################
            // # Try to vertex the two muons to get dimuon vertex #
            // ####################################################
            std::vector<RefCountedKinematicParticle> muonParticles;
            muonParticles.push_back(partFactory.particle(muTrackmTT, muonMass,chi,ndf, muonMassErr));
            muonParticles.push_back(partFactory.particle(muTrackpTT, muonMass,chi,ndf, muonMassErr));
        
            RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles); 
            
	    if (mumuVertexFitTree->isValid() == false)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
              continue; 
            }
        
            mumuVertexFitTree->movePointerToTheTop();
            RefCountedKinematicVertex mumu_KV   = mumuVertexFitTree->currentDecayVertex();
            if (mumu_KV->vertexIsValid() == false)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
              continue;
            }
            if (TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))) < CLMUMUVTX)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad vtx CL from mu+ mu- fit: " << TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))) << std::endl;
              continue;
            }

            RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();

            // ######################################################
            // # Compute the distance between mumu vtx and BeamSpot #
            // ######################################################
            
            double MuMuLSBS;
            double MuMuLSBSErr;
            Utility->computeLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
                                beamSpot.position().x(),beamSpot.position().y(),0.0,
                                mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
                                mumu_KV->error().matrix()(0,1),0.0,0.0,
                                beamSpot.covariance()(0,0),beamSpot.covariance()(1,1),0.0,
                                beamSpot.covariance()(0,1),0.0,0.0,
                                &MuMuLSBS,&MuMuLSBSErr);
            if (MuMuLSBS/MuMuLSBSErr < LSMUMUBS)     
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad mumu L/sigma with respect to BeamSpot: " << MuMuLSBS << "+/-" << MuMuLSBSErr << std::endl;
              continue;
            }
            // ###################################################################
            // # Compute cos(alpha) between mumu momentum and mumuVtx - BeamSpot #
            // ###################################################################
            double MuMuCosAlphaBS;
            double MuMuCosAlphaBSErr;
            Utility->computeCosAlpha (mumu_KP->currentState().globalMomentum().x(),mumu_KP->currentState().globalMomentum().y(),0.0,
                                      mumu_KV->position().x() - beamSpot.position().x(),mumu_KV->position().y() - beamSpot.position().y(),0.0,
                                      mumu_KP->currentState().kinematicParametersError().matrix()(3,3),mumu_KP->currentState().kinematicParametersError().matrix()(4,4),0.0,
                                      mumu_KP->currentState().kinematicParametersError().matrix()(3,4),0.0,0.0,
                                      mumu_KV->error().cxx() + beamSpot.covariance()(0,0),mumu_KV->error().cyy() + beamSpot.covariance()(1,1),0.0,
                                      mumu_KV->error().matrix()(0,1) + beamSpot.covariance()(0,1),0.0,0.0,
                                      &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);
            if (MuMuCosAlphaBS < COSALPHAMUMUBS)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad mumu cos(alpha) with respect to BeamSpot: " << MuMuCosAlphaBS << "+/-" << MuMuCosAlphaBSErr << std::endl;
              continue;
            }	     
	 */


	//}



   }
 
    theTree->Fill();
    NTuple->ClearNTuple();

}


// ------------ method called once each job just before starting event loop  ------------
void
NTuplizer::beginJob()
{
    edm::Service<TFileService> outfile_;
    theTree = outfile_->make<TTree>("B0KstMuMuNTuple","B0KstMuMuNTuple");
    NTuple->MakeTreeBranches(theTree);
}

// ------------ method called once each job just after ending the event loop  ------------
void
NTuplizer::endJob()
{

    theTree->GetDirectory()->cd();
    theTree->Write();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NTuplizer);
