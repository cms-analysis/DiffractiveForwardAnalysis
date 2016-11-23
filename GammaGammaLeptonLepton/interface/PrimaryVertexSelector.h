#ifndef DiffractiveForwardAnalysis_PrimaryVertexSelector_h
#define DiffractiveForwardAnalysis_PrimaryVertexSelector_h

// system include files
#include <fstream>
#include <memory>
#include <vector>

// Muons collection
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// Electrons collection
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// Vertices collection
#include "DataFormats/Common/interface/RefToBase.h" 

#include "TLorentzVector.h"

//
// class declaration
//

class PrimaryVertexSelector {
  public:
    typedef std::vector< std::pair<int,reco::TrackRef> > MatchedLeptonsMap;

  public:
    explicit PrimaryVertexSelector(const std::map<int,TLorentzVector>&, const std::map<int,TLorentzVector>&);
    inline ~PrimaryVertexSelector() {;}

    void feedTracks(const reco::Vertex::trackRef_iterator&, const reco::Vertex::trackRef_iterator&);

    inline MatchedLeptonsMap matchedElectrons() const { return matchedElectrons_; }
    int matchedElectron(const reco::TrackRef&) const;

    inline MatchedLeptonsMap matchedMuons() const { return matchedMuons_; }
    int matchedMuon(const reco::TrackRef&) const;

  private:
    typedef std::map<int,TLorentzVector> LeptonsMap;

    LeptonsMap muonMomenta_, electronMomenta_;
    MatchedLeptonsMap matchedMuons_, matchedElectrons_;
};

#endif
