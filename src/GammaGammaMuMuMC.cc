// -*- C++ -*-
//
// Package:    GammaGammaMuMuMC
// Class:      GammaGammaMuMuMC
// 
/**\class GammaGammaMuMuMC GammaGammaMuMuMC.cc GammaGammaLeptonLepton/GammaGammaMuMuMC/src/GammaGammaMuMuMC.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jonathan Hollar
//         Created:  Wed Sep 20 10:08:38 BST 2006
// $Id: GammaGammaMuMuMC.cc,v 1.5 2007/12/04 13:13:25 jjhollar Exp $
//
//


#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaMuMuMC.h"

#include <FWCore/MessageLogger/interface/MessageLogger.h> 
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// C++
#include <memory>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

using namespace edm;
using namespace std;
using namespace HepMC;
using namespace reco;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GammaGammaMuMuMC::GammaGammaMuMuMC(const edm::ParameterSet& pset)
{
   //now do what ever initialization is needed
  rootfilename       = pset.getUntrackedParameter<std::string>("outfilename","test.root");
  fillallmc          = pset.getParameter<bool>("FillAllMCParticles");

  //  nEvt=0;
  MUONMAX=10;
  PROTMAX=30;
  
  thefile = new TFile(rootfilename.c_str(),"recreate");
  thefile->cd();
  thetree= new TTree("ntp1","ntp1");

  thetree->Branch("nGenProtCand",&nGenProtCand,"nGenProtCand/I");
  thetree->Branch("GenProtCand_px",GenProtCand_px,"GenProtCand_px[nGenProtCand]/D");
  thetree->Branch("GenProtCand_py",GenProtCand_py,"GenProtCand_py[nGenProtCand]/D");
  thetree->Branch("GenProtCand_pz",GenProtCand_pz,"GenProtCand_pz[nGenProtCand]/D");
  thetree->Branch("GenProtCand_e",GenProtCand_e,"GenProtCand_e[nGenProtCand]/D");

  thetree->Branch("nGenMuonCand",&nGenMuonCand,"nGenMuonCand/I");
  thetree->Branch("GenMuonCand_px",GenMuonCand_px,"GenMuonCand_px[nGenMuonCand]/D");
  thetree->Branch("GenMuonCand_py",GenMuonCand_py,"GenMuonCand_py[nGenMuonCand]/D");
  thetree->Branch("GenMuonCand_pz",GenMuonCand_pz,"GenMuonCand_pz[nGenMuonCand]/D");
  thetree->Branch("GenMuonCand_p",GenMuonCand_p,"GenMuonCand_p[nGenMuonCand]/D");
  thetree->Branch("GenMuonCand_pt",GenMuonCand_pt,"GenMuonCand_pt[nGenMuonCand]/D");
  thetree->Branch("GenMuonCand_eta",GenMuonCand_eta,"GenMuonCand_eta[nGenMuonCand]/D");
  thetree->Branch("GenMuonCand_phi",GenMuonCand_phi,"GenMuonCand_phi[nGenMuonCand]/D");
  thetree->Branch("GenMuonCand_charge",GenMuonCand_charge,"GenMuonCand_charge[nGenMuonCand]/I");

  if(fillallmc == true)
    {
      thetree->Branch("nMCPar",&nMCPar,"nMCPar/I"); 
      thetree->Branch("MCPar_px",MCPar_px,"MCPar_px[nMCPar]/D"); 
      thetree->Branch("MCPar_py",MCPar_py,"MCPar_py[nMCPar]/D"); 
      thetree->Branch("MCPar_pz",MCPar_pz,"MCPar_pz[nMCPar]/D"); 
      thetree->Branch("MCPar_e",MCPar_e,"MCPar_e[nMCPar]/D"); 
      thetree->Branch("MCPar_eta",MCPar_eta,"MCPar_eta[nMCPar]/D"); 
      thetree->Branch("MCPar_phi",MCPar_phi,"MCPar_phi[nMCPar]/D"); 
      thetree->Branch("MCPar_pdgid",MCPar_pdgid,"MCPar_pdgid[nMCPar]/I");
      thetree->Branch("MCPar_status",MCPar_status,"MCPar_status[nMCPar]/I");
    }

  thetree->Branch("HitInZDC",&HitInZDC,"HitInZDC/I");
  thetree->Branch("HitInCastor",&HitInCastor,"HitInCastor/I");

  thetree->Branch("GenMuMu_mass",&GenMuMu_mass,"GenMuMu_mass/D");
  thetree->Branch("GenMuMu_dphi",&GenMuMu_dphi,"GenMuMu_dphi/D");

  nEvt = 0;
}


GammaGammaMuMuMC::~GammaGammaMuMuMC()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GammaGammaMuMuMC::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  nGenMuonCand=0;
  nGenProtCand=0;
  nMCPar=0;
  HitInZDC=0;
  HitInCastor=0;

  GenMuMu_mass = -1;
  GenMuMu_dphi = -1;

  bool passed = true;

  nEvt++;
  if((nEvt%10==0 && nEvt<=100)||(nEvt%100==0 && nEvt>100))
    std::cout << "reading event " << nEvt << std::endl;
  
  // step 1: fill some basic MC information into the root tree
  Handle<CandidateCollection> genParticles;
  event.getByLabel( "genParticleCandidates", genParticles );

  for ( size_t i = 0; i < genParticles->size() && i < 500; ++ i ) 
    {
      const Candidate & p = (*genParticles)[ i ];

       MCPar_status[nMCPar]= p.status();
       if(p.status() == 3)
	 continue;
       MCPar_px[nMCPar]=p.px();
       MCPar_py[nMCPar]=p.py();
       MCPar_pz[nMCPar]=p.pz();
       MCPar_phi[nMCPar]=atan2(MCPar_py[nMCPar],MCPar_px[nMCPar]);
       MCPar_eta[nMCPar] = p.eta();
       MCPar_pdgid[nMCPar]=p.pdgId();
       MCPar_mass[nMCPar]=p.mass();
       MCPar_e[nMCPar] = sqrt(MCPar_mass[nMCPar]*MCPar_mass[nMCPar] + (MCPar_px[nMCPar]*MCPar_px[nMCPar] + MCPar_py[nMCPar]*MCPar_py[nMCPar] + MCPar_pz[nMCPar]*MCPar_pz[nMCPar])); 

       if(abs(p.pdgId()) == 13)
	{
	  GenMuonCand_px[nGenMuonCand] = MCPar_px[nMCPar];
	  GenMuonCand_py[nGenMuonCand] = MCPar_py[nMCPar];
	  GenMuonCand_pz[nGenMuonCand] = MCPar_pz[nMCPar];
	  GenMuonCand_phi[nGenMuonCand] = MCPar_phi[nMCPar];
	  GenMuonCand_charge[nGenMuonCand] = (int)(MCPar_pdgid[nMCPar]/13);
	  GenMuonCand_eta[nGenMuonCand] = MCPar_eta[nMCPar];
	  GenMuonCand_p[nGenMuonCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+ 
					     (MCPar_py[nMCPar]*MCPar_py[nMCPar])+
					     (MCPar_pz[nMCPar]*MCPar_pz[nMCPar]));
	  GenMuonCand_pt[nGenMuonCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+ 
					      (MCPar_py[nMCPar]*MCPar_py[nMCPar]));
	  
	  nGenMuonCand++;
	}
      if (p.pdgId() == 2212)
	{
	  GenProtCand_px[nGenProtCand] = MCPar_px[nMCPar];
	   GenProtCand_py[nGenProtCand] = MCPar_py[nMCPar];
	   GenProtCand_pz[nGenProtCand] = MCPar_pz[nMCPar];
	   GenProtCand_e[nGenProtCand] = MCPar_e[nMCPar];
	   
	   nGenProtCand++; 
	}       

      if(MCPar_pdgid[nMCPar] == 22 && abs(MCPar_eta[nMCPar]) > 8.6 && MCPar_e[nMCPar] > 20.0) 
	HitInZDC++;
      if(MCPar_pdgid[nMCPar] == 2112 && abs(MCPar_eta[nMCPar]) > 8.6 && MCPar_e[nMCPar] > 50.0)
	HitInZDC++;
      if((MCPar_pdgid[nMCPar] != 22 && MCPar_pdgid[nMCPar] != 2112) && (abs(MCPar_eta[nMCPar]) > 5.2 && abs(MCPar_eta[nMCPar]) < 6.6))
	HitInCastor++;

      nMCPar++;
    }

  // Two opposite sign dimuons
  if(nGenMuonCand == 2 && (GenMuonCand_charge[0]*GenMuonCand_charge[1] < 0))
    {
      double mass = pow(GenMuonCand_p[0]+GenMuonCand_p[1],2);
      mass-=pow(GenMuonCand_px[0]+GenMuonCand_px[1],2);
      mass-=pow(GenMuonCand_py[0]+GenMuonCand_py[1],2);
      mass-=pow(GenMuonCand_pz[0]+GenMuonCand_pz[1],2);
      GenMuMu_mass = sqrt(mass);      
      
      double dphi = abs(GenMuonCand_phi[0]-GenMuonCand_phi[1]);
      if(dphi > 3.14159)
	dphi = (2.0*3.14159) - dphi;       
      GenMuMu_dphi = dphi;
    }
  
  if(passed == true)
    thetree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaMuMuMC::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaMuMuMC::endJob() {
  thefile->Write();
  thefile->Close();
}
  
