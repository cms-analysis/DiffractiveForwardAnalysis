// -*- C++ -*-
//
// Package:    GammaGammaLeptonLepton
// Class:      PrimaryVertexSelector
// 
/**\class PrimaryVertexSelector PrimaryVertexSelector.cc DiffractiveForwardAnalysis/GammaGammaLeptonLepton/src/PrimaryVertexSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme,40 4-B20,+41227671567,
//         Created:  Thu Sep 13 15:17:14 CET 2012
// $Id: PrimaryVertexSelector.cc,v 1.3 2013/04/28 08:40:45 lforthom Exp $
//
//

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/PrimaryVertexSelector.h"

//
// constructors and destructor
//
PrimaryVertexSelector::PrimaryVertexSelector(const std::map<int,TLorentzVector>& mu, const std::map<int,TLorentzVector>& ele) :
  muonMomenta_(mu), electronMomenta_(ele)
{
}

/**
 * \brief Matches a track arising from the vertex with a lepton track from the
 *  internal collection
 */
void
PrimaryVertexSelector::feedTracks(const reco::Vertex::trackRef_iterator& begin, const reco::Vertex::trackRef_iterator& end)
{
  matchedElectrons_.clear();
  matchedMuons_.clear();

  const float dr_max = 0.1;

  for (reco::Vertex::trackRef_iterator trk_it=begin; trk_it!=end; trk_it++) {
    const reco::TrackRef trk = trk_it->castTo<reco::TrackRef>();

    const TVector3 trk_vec(trk->px(), trk->py(), trk->pz());

    // look at electron matching
    for (LeptonsMap::const_iterator mu=muonMomenta_.begin(); mu!=muonMomenta_.end(); mu++) {
      if (trk_vec.DeltaR(mu->second.Vect())>dr_max) continue;
      matchedMuons_.push_back(std::pair<int,reco::TrackRef>(mu->first, trk));
    }

    // then look at muon matching
    for (LeptonsMap::const_iterator ele=electronMomenta_.begin(); ele!=electronMomenta_.end(); ele++) {
      if (trk_vec.DeltaR(ele->second.Vect())>dr_max) continue;
      matchedElectrons_.push_back(std::pair<int,reco::TrackRef>(ele->first, trk));
    }
  }
}

int
PrimaryVertexSelector::matchedElectron(const reco::TrackRef& trk) const
{
  for (MatchedLeptonsMap::const_iterator it=matchedElectrons_.begin(); it!=matchedElectrons_.end(); it++) {
    if (it->second==trk);
  }
  return -1;
}

int
PrimaryVertexSelector::matchedMuon(const reco::TrackRef& trk) const
{
  for (MatchedLeptonsMap::const_iterator it=matchedMuons_.begin(); it!=matchedMuons_.end(); it++) {
    if (it->second==trk);
  }
  return -1;
}
