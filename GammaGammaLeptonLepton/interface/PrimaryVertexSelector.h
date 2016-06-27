#ifndef DiffractiveForwardAnalysis_PrimaryVertexSelector_h
#define DiffractiveForwardAnalysis_PrimaryVertexSelector_h

// system include files
#include <fstream>
#include <memory>
#include <map>

// Muons collection
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// Electrons collection
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// Vertices collection
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h" 
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include <TVector3.h>
#include <TLorentzVector.h>

//
// class declaration
//

class PrimaryVertexSelector : public reco::Vertex {
  public:
    explicit PrimaryVertexSelector(std::vector<std::string>&, std::map<int,TLorentzVector>&, std::map<int,TLorentzVector>&);
    ~PrimaryVertexSelector();
    void SetPosition(double, double, double);
    int AddTrack(const reco::TrackRef&, TString&);
    inline int Electrons() { return nMatchedElectrons; }
    inline int Muons() { return nMatchedMuons; }
    double dZ(TVector3, int);
    TVector3 Position;
    int nTracks, nMatchedTracks, nUnmatchedTracks;
    std::vector<int> MatchedMuons, MatchedElectrons;
  private:
    unsigned int i;
    int nMatchedMuons, nMatchedElectrons;
    bool FetchMuons, FetchElectrons;
    std::map<int,TLorentzVector> MuonMomenta;
    std::map<int,TLorentzVector> ElectronMomenta;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

#endif
