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
PrimaryVertexSelector::PrimaryVertexSelector(std::vector<std::string>& _leptonsType, std::map<int,TLorentzVector>& _muonsMomenta, std::map<int,TLorentzVector>& _electronsMomenta) :
  nTracks(0),
  nMatchedTracks(0),
  nUnmatchedTracks(0),
  nMatchedMuons(0),
  nMatchedElectrons(0),
  FetchMuons(false),
  FetchElectrons(false)
{
  //LeptonsType = _leptonsType;
  for (i=0; i<_leptonsType.size(); i++) {
    if (_leptonsType[i]=="Muon") {
      FetchMuons = true;
    }
    else if (_leptonsType[i]=="Electron") {
      FetchElectrons = true;
    }
  }
  MuonMomenta = _muonsMomenta;
  ElectronMomenta = _electronsMomenta;
}

PrimaryVertexSelector::~PrimaryVertexSelector() {
}

void
PrimaryVertexSelector::SetPosition(double _x, double _y, double _z)
{
  Position.SetXYZ(_x, _y, _z);
#ifdef DEBUG
  std::cout << "[PrimaryVertexSelector] SetPosition :: Vertex located at (" << Position.x() << ", " << Position.y() << ", " << Position.z() << ")" << std::endl;
#endif
}

/**
 * \brief Matches a track arising from the vertex with a lepton track from the
 *  internal collection
 */
int
PrimaryVertexSelector::AddTrack(const reco::TrackRef& _track, TString& _leptonType)
{
  nTracks++; // total number of tracks matched with the vertex
  std::map<int,TLorentzVector>::iterator lep;
  for (lep=MuonMomenta.begin(); lep!=MuonMomenta.end(); lep++) {
    if (fabs(_track->p()-lep->second.P())>.01) continue;
    if (fabs(_track->pt()-lep->second.Pt())>.01) continue;
    if (fabs(_track->eta()-lep->second.Eta())>.01) continue;
    if (fabs(_track->phi()-lep->second.Phi())>.01) continue;
    _leptonType = "muon";
    MatchedMuons.push_back(lep->first);
    nMatchedMuons++;
    nMatchedTracks++;
    return lep->first;
  }
  for (lep=ElectronMomenta.begin(); lep!=ElectronMomenta.end(); lep++) {
    if (fabs(_track->p()-lep->second.P())>.01) continue;
    if (fabs(_track->pt()-lep->second.Pt())>.01) continue;
    if (fabs(_track->eta()-lep->second.Eta())>.01) continue;
    if (fabs(_track->phi()-lep->second.Phi())>.01) continue;
    _leptonType = "electron";
    MatchedElectrons.push_back(lep->first);
    nMatchedElectrons++;
    nMatchedTracks++;
    return lep->first;
  }
  nUnmatchedTracks++;
  return -1;
}

double
PrimaryVertexSelector::dZ(TVector3 _vmu, int _muind)
{
  TLorentzVector m(MuonMomenta[_muind]);
  return (_vmu.Z()-Position.Z())-((_vmu.X()-Position.X())*m.Px()+(_vmu.Y()-Position.Y())*m.Py())/m.Pt()*m.Pz()/m.Pt();
}
