// -*- C++ -*-
//
// Package:    DiffractiveForwardAnalysis/SinglePhotonTreeProducer
// Class:      SinglePhotonTreeProducer
//
/**\class SinglePhotonTreeProducer SinglePhotonTreeProducer.cc DiffractiveForwardAnalysis/SinglePhotonTreeProducer/plugins/SinglePhotonTreeProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme
//         Created:  Tue, 24 Jul 2018 12:24:16 GMT
//
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/SinglePhotonEvent.h"

class SinglePhotonTreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
  public:
    explicit SinglePhotonTreeProducer( const edm::ParameterSet& );
    ~SinglePhotonTreeProducer() {}

    static void fillDescriptions( edm::ConfigurationDescriptions& descriptions );

  private:
    virtual void beginRun( const edm::Run&, const edm::EventSetup& );

    virtual void beginJob() override;
    virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
    virtual void endJob() override {}

    TTree* tree_;
    gggx::SinglePhotonEvent evt_;

    std::string hltMenuLabel_;
    std::vector<std::string> triggersList_;
    bool runOnMC_;
    unsigned int minPhotonMult_, minFwdTrks_;

    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
    edm::EDGetTokenT<edm::View<pat::Photon> > photonToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > phoWP80IdMapToken_, phoWP90IdMapToken_;
    edm::EDGetTokenT<edm::View<CTPPSLocalTrackLite> > ppsTracksToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> > verticesToken_;
    edm::EDGetTokenT<edm::View<pat::Jet> > jetsToken_;
    edm::EDGetTokenT<edm::View<pat::MET> > metsToken_;

    HLTConfigProvider hltConfig_;
    HLTPrescaleProvider hltPrescale_;
};

SinglePhotonTreeProducer::SinglePhotonTreeProducer( const edm::ParameterSet& iConfig ) :
  hltMenuLabel_  ( iConfig.getParameter<std::string>( "hltMenuTag" ) ),
  triggersList_  ( iConfig.getParameter<std::vector<std::string> >( "triggersList" ) ),
  runOnMC_       ( iConfig.getParameter<bool>( "runOnMC" ) ),
  minPhotonMult_ ( iConfig.getParameter<unsigned int>( "minPhotonMult" ) ),
  minFwdTrks_    ( iConfig.getParameter<unsigned int>( "minFwdTrks" ) ),
  triggerResultsToken_( consumes<edm::TriggerResults>            ( iConfig.getParameter<edm::InputTag>( "triggerResults" ) ) ),
  photonToken_        ( consumes<edm::View<pat::Photon> >        ( iConfig.getParameter<edm::InputTag>( "photonsTag" ) ) ),
  phoWP80IdMapToken_  ( consumes<edm::ValueMap<bool> >           ( iConfig.getParameter<edm::InputTag>( "phoWP80IdMap" ) ) ),
  phoWP90IdMapToken_  ( consumes<edm::ValueMap<bool> >           ( iConfig.getParameter<edm::InputTag>( "phoWP90IdMap" ) ) ),
  ppsTracksToken_     ( consumes<edm::View<CTPPSLocalTrackLite> >( iConfig.getParameter<edm::InputTag>( "ppsTracksTag" ) ) ),
  verticesToken_      ( consumes<edm::View<reco::Vertex> >       ( iConfig.getParameter<edm::InputTag>( "verticesTag" ) ) ),
  jetsToken_          ( consumes<edm::View<pat::Jet> >           ( iConfig.getParameter<edm::InputTag>( "jetsTag" ) ) ),
  metsToken_          ( consumes<edm::View<pat::MET> >           ( iConfig.getParameter<edm::InputTag>( "metsTag" ) ) ),
  hltPrescale_        ( iConfig, consumesCollector(), *this )
{
  usesResource( "TFileService" );
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>( "ntp1", "ntp1" );
}


void
SinglePhotonTreeProducer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  evt_.clear();

  // Run and BX information
  evt_.BX = iEvent.bunchCrossing();
  evt_.Run = iEvent.id().run();
  evt_.LumiSection = iEvent.luminosityBlock();
  evt_.EventNum = iEvent.id().event();

  //--- trigger information
  edm::Handle<edm::TriggerResults> hltResults;
  iEvent.getByToken( triggerResultsToken_, hltResults );
  const edm::TriggerNames& trigNames = iEvent.triggerNames( *hltResults );

  evt_.HLT_Name.reserve( trigNames.size() );
  for ( const auto& sel : triggersList_ ) {
    short trig_id = -1;
    const std::string trig_name = ( sel.size() > 1 && sel[sel.size()-1] == '*' )
      ? sel.substr( 0, sel.size()-1 ) // remove trailing '*'
      : sel;
    for ( unsigned int i = 0; i < trigNames.size(); ++i )
      if ( trigNames.triggerNames().at( i ).find( trig_name ) != std::string::npos ) {
        trig_id = i;
        break;
      }
    if ( trig_id < 0 )
      continue;
    evt_.HLT_Name[evt_.nHLT] = trigNames.triggerNames().at( trig_id );
    evt_.HLT_Accept[evt_.nHLT] = hltResults->accept( trig_id );
    evt_.HLT_Prescl[evt_.nHLT] = 0.;
    // extract prescale value for this path
    if ( !runOnMC_ ) {
      int prescale_set = hltPrescale_.prescaleSet( iEvent, iSetup );
      evt_.HLT_Prescl[evt_.nHLT] = ( prescale_set < 0 )
        ? 0.
        : hltConfig_.prescaleValue(prescale_set, trigNames.triggerNames().at( trig_id ) );
    }
    evt_.nHLT++;
  }

  //--- photons collection retrieval
  edm::Handle<edm::View<pat::Photon> > photonColl;
  iEvent.getByToken( photonToken_, photonColl );

  //--- photon identification
  edm::Handle<edm::ValueMap<bool> > wp80_id_decisions, wp90_id_decisions;
  iEvent.getByToken( phoWP80IdMapToken_, wp80_id_decisions );
  iEvent.getByToken( phoWP90IdMapToken_, wp90_id_decisions );

  for ( unsigned int i = 0; i < photonColl->size(); ++i ) {
    const edm::Ptr<pat::Photon> photon = photonColl->ptrAt( i );

    evt_.PhotonCand_pt[evt_.nPhotonCand] = photon->pt();
    evt_.PhotonCand_eta[evt_.nPhotonCand] = photon->eta();
    evt_.PhotonCand_phi[evt_.nPhotonCand] = photon->phi();
    evt_.PhotonCand_e[evt_.nPhotonCand] = photon->energy();
    evt_.PhotonCand_r9[evt_.nPhotonCand] = photon->r9();

    evt_.PhotonCand_wp80id[evt_.nPhotonCand] = wp80_id_decisions->operator[]( photon );
    evt_.PhotonCand_wp90id[evt_.nPhotonCand] = wp90_id_decisions->operator[]( photon );

    evt_.PhotonCand_electronveto[evt_.nPhotonCand] = photon->passElectronVeto();
    evt_.PhotonCand_pixelseed[evt_.nPhotonCand] = photon->hasPixelSeed();

    evt_.PhotonCand_drtrue[evt_.nPhotonCand] = -999.;
    evt_.PhotonCand_detatrue[evt_.nPhotonCand] = -999.;
    evt_.PhotonCand_dphitrue[evt_.nPhotonCand] = -999.;

    if ( runOnMC_ ) {
      double photdr = 999., photdeta = 999., photdphi = 999.;
      double endphotdr = 999., endphotdeta = 999., endphotdphi = 999.;
      for ( unsigned int j = 0; j < evt_.nGenPhotCand; ++j ) { // matching with the 'true' photon object from MC
        photdeta = ( evt_.PhotonCand_eta[evt_.nPhotonCand]-evt_.GenPhotCand_eta[j] );
        photdphi = ( evt_.PhotonCand_phi[evt_.nPhotonCand]-evt_.GenPhotCand_phi[j] );
        photdr = sqrt( photdeta*photdeta + photdphi*photdphi );
        if ( photdr < endphotdr ) {
          endphotdr = photdr;
          endphotdeta = photdeta;
          endphotdphi = photdphi;
        }
      }
      evt_.PhotonCand_detatrue[evt_.nPhotonCand] = endphotdeta;
      evt_.PhotonCand_dphitrue[evt_.nPhotonCand] = endphotdphi;
      evt_.PhotonCand_drtrue[evt_.nPhotonCand] = endphotdr;
    }

    evt_.nPhotonCand++;
  }
  //--- do not store if minimal number of photons is not reached
  if ( evt_.nPhotonCand < minPhotonMult_ )
    return;

  //--- PPS local tracks
  edm::Handle<edm::View<CTPPSLocalTrackLite> > ppsTracks;
  iEvent.getByToken( ppsTracksToken_, ppsTracks );
  for ( const auto& trk : *ppsTracks ) {
    if ( evt_.nFwdTrkCand >= gggx::SinglePhotonEvent::MAX_FWDTRKCAND-1 ) {
      edm::LogWarning( "SinglePhotonTreeProducer" )
        << "maximum number of local tracks in RPs is reached! increase MAX_FWDTRKCAND="
        << gggx::SinglePhotonEvent::MAX_FWDTRKCAND;
      break;
    }
    const CTPPSDetId detid( trk.getRPId() );
    evt_.FwdTrkCand_arm[evt_.nFwdTrkCand] = detid.arm();
    evt_.FwdTrkCand_station[evt_.nFwdTrkCand] = detid.station();
    evt_.FwdTrkCand_pot[evt_.nFwdTrkCand] = detid.rp();
    evt_.FwdTrkCand_x[evt_.nFwdTrkCand] = trk.getX()/1.e3;
    evt_.FwdTrkCand_y[evt_.nFwdTrkCand] = trk.getY()/1.e3;
    evt_.FwdTrkCand_xSigma[evt_.nFwdTrkCand] = trk.getXUnc()/1.e3;
    evt_.FwdTrkCand_ySigma[evt_.nFwdTrkCand] = trk.getYUnc()/1.e3;
    evt_.nFwdTrkCand++;
    LogDebug( "SinglePhotonTreeProducer" )
      << "Proton track candidate with origin: ( " << trk.getX() << ", " << trk.getY() << " ) extracted!";
  }
  //--- do not store if minimal number of PPS tracks is not reached
  if ( evt_.nFwdTrkCand < minFwdTrks_ )
    return;

  //--- vertex collection
  edm::Handle<edm::View<reco::Vertex> > verticesColl;
  iEvent.getByToken( verticesToken_, verticesColl );

  for ( unsigned int i = 0; i < verticesColl->size() && evt_.nPrimVertexCand < gggx::SinglePhotonEvent::MAX_VTX; ++i ) {
    const edm::Ptr<reco::Vertex> vertex = verticesColl->ptrAt( i );

    evt_.PrimVertexCand_x[evt_.nPrimVertexCand] = vertex->x();
    evt_.PrimVertexCand_y[evt_.nPrimVertexCand] = vertex->y();
    evt_.PrimVertexCand_z[evt_.nPrimVertexCand] = vertex->z();
    evt_.PrimVertexCand_tracks[evt_.nPrimVertexCand] = vertex->nTracks();
    evt_.PrimVertexCand_chi2[evt_.nPrimVertexCand] = vertex->chi2();
    evt_.PrimVertexCand_ndof[evt_.nPrimVertexCand] = vertex->ndof();
    evt_.nPrimVertexCand++;
  }

  //--- jets collection
  edm::Handle<edm::View<pat::Jet> > jetColl;
  iEvent.getByToken( jetsToken_, jetColl );

  double totalJetEnergy = 0., HEJet_pt = 0., HEJet_eta = 0., HEJet_phi = 0., HEJet_e = 0.;

  for ( unsigned int i = 0; i < jetColl->size() && evt_.nJetCand < gggx::SinglePhotonEvent::MAX_JETS; ++i ) {
    const edm::Ptr<pat::Jet> jet = jetColl->ptrAt( i );

    evt_.JetCand_e[evt_.nJetCand] = jet->energy();
    evt_.JetCand_pt[evt_.nJetCand] = jet->pt();
    evt_.JetCand_eta[evt_.nJetCand] = jet->eta();
    evt_.JetCand_phi[evt_.nJetCand] = jet->phi();
    totalJetEnergy += jet->energy();
    //--- find the highest energy jet
    if ( evt_.JetCand_e[evt_.nJetCand] > HEJet_e ) {
      HEJet_e = evt_.JetCand_e[evt_.nJetCand];
      HEJet_pt = evt_.JetCand_pt[evt_.nJetCand];
      HEJet_eta = evt_.JetCand_eta[evt_.nJetCand];
      HEJet_phi = evt_.JetCand_phi[evt_.nJetCand];
    }
    evt_.nJetCand++;
  }
  evt_.HighestJet_pt = HEJet_pt;
  evt_.HighestJet_eta = HEJet_eta;
  evt_.HighestJet_phi = HEJet_phi;
  evt_.HighestJet_e = HEJet_e;
  evt_.SumJet_e = totalJetEnergy;

  //--- missing ET
  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByToken( metsToken_, mets );
  const edm::View<pat::MET>* metColl = mets.product();
  edm::View<pat::MET>::const_iterator met = metColl->begin();

  evt_.Etmiss = met->et();
  evt_.Etmiss_phi = met->phi();
  evt_.Etmiss_significance = met->significance();


  tree_->Fill();
}

void
SinglePhotonTreeProducer::beginRun( const edm::Run& iRun, const edm::EventSetup& iSetup )
{
  //--- HLT part
  bool changed = true;
  if ( !hltPrescale_.init( iRun, iSetup, hltMenuLabel_, changed ) )
    throw cms::Exception( "GammaGammaLL" ) << " prescales extraction failure with process name " << hltMenuLabel_;
  // Initialise HLTConfigProvider
  hltConfig_ = hltPrescale_.hltConfigProvider();
  if ( !hltConfig_.init( iRun, iSetup, hltMenuLabel_, changed ) )
    throw cms::Exception( "GammaGammaLL" ) << " config extraction failure with process name " << hltMenuLabel_;
  else if ( hltConfig_.size() == 0 )
    edm::LogError( "GammaGammaLL" ) << "HLT config size error";
}

void
SinglePhotonTreeProducer::beginJob()
{
  evt_.attach( tree_, runOnMC_ );
}

void
SinglePhotonTreeProducer::fillDescriptions( edm::ConfigurationDescriptions& descriptions )
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault( desc );
}

//define this as a plug-in
DEFINE_FWK_MODULE( SinglePhotonTreeProducer );
