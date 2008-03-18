// -*- C++ -*-
//
// Package:    GammaGammaEEMC
// Class:      GammaGammaEEMC
// 
/**\class GammaGammaEEMC GammaGammaEEMC.cc GammaGammaLeptonLepton/GammaGammaEEMC/src/GammaGammaEEMC.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jonathan Hollar
//         Created:  Wed Sep 20 10:08:38 BST 2006
// $Id: GammaGammaEEMC.cc,v 1.5 2007/12/04 13:13:25 jjhollar Exp $
//
//


#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaEEMC.h"

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
GammaGammaEEMC::GammaGammaEEMC(const edm::ParameterSet& pset)
{
   //now do what ever initialization is needed
  rootfilename       = pset.getUntrackedParameter<std::string>("outfilename","test.root");
  fillallmc          = pset.getParameter<bool>("FillAllMCParticles");

  //  nEvt=0;
  ELEMAX=10;
  PROTMAX=30;
  
  thefile = new TFile(rootfilename.c_str(),"recreate");
  thefile->cd();
  thetree= new TTree("ntp1","ntp1");

  thetree->Branch("nGenEleCand",&nGenEleCand,"nGenEleCand/I");
  thetree->Branch("GenEleCand_px",GenEleCand_px,"GenEleCand_px[nGenEleCand]/D");
  thetree->Branch("GenEleCand_py",GenEleCand_py,"GenEleCand_py[nGenEleCand]/D");
  thetree->Branch("GenEleCand_pz",GenEleCand_pz,"GenEleCand_pz[nGenEleCand]/D");
  thetree->Branch("GenEleCand_p",GenEleCand_p,"GenEleCand_p[nGenEleCand]/D");
  thetree->Branch("GenEleCand_e",GenEleCand_e,"GenEleCand_e[nGenEleCand]/D");
  thetree->Branch("GenEleCand_pt",GenEleCand_pt,"GenEleCand_pt[nGenEleCand]/D");
  thetree->Branch("GenEleCand_eta",GenEleCand_eta,"GenEleCand_eta[nGenEleCand]/D");
  thetree->Branch("GenEleCand_phi",GenEleCand_phi,"GenEleCand_phi[nGenEleCand]/D");
  thetree->Branch("GenEleCand_charge",GenEleCand_charge,"GenEleCand_charge[nGenEleCand]/D");

  thetree->Branch("nGenProtCand",&nGenProtCand,"nGenProtCand/I");
  thetree->Branch("GenProtCand_px",GenProtCand_px,"GenProtCand_px[nGenProtCand]/D");
  thetree->Branch("GenProtCand_py",GenProtCand_py,"GenProtCand_py[nGenProtCand]/D");
  thetree->Branch("GenProtCand_pz",GenProtCand_pz,"GenProtCand_pz[nGenProtCand]/D");
  thetree->Branch("GenProtCand_e",GenProtCand_e,"GenProtCand_e[nGenProtCand]/D");

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
  
  thetree->Branch("GenElEl_mass",&GenElEl_mass,"GenElEl_mass/D");
  thetree->Branch("GenElEl_dphi",&GenElEl_dphi,"GenElEl_dphi/D");
  
}


GammaGammaEEMC::~GammaGammaEEMC()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GammaGammaEEMC::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  nGenEleCand=0;
  nGenProtCand=0;
  nMCPar=0;
  HitInZDC=0;
  HitInCastor=0;

  GenElEl_mass = -1;
  GenElEl_dphi = -1;

  bool passed = true;

  nEvt++;
  if((nEvt%10==0 && nEvt<=100)||(nEvt%100==0 && nEvt>100))
     std::cout << "reading event " << nEvt << std::endl;
  
  
  // step 1: fill some basic MC information into the root tree
  Handle<CandidateCollection> genParticles;
  event.getByLabel( "genParticleCandidates", genParticles );
  

   for( size_t i = 0; i < genParticles->size() && i < 500; ++ i )
    {
      const Candidate & p = (*genParticles)[ i ];

       MCPar_status[nMCPar]=p.status();
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


       if (abs(p.pdgId()) == 11)
	{
	  GenEleCand_px[nGenEleCand] = MCPar_px[nMCPar];
	  GenEleCand_py[nGenEleCand] = MCPar_py[nMCPar];
	  GenEleCand_pz[nGenEleCand] = MCPar_pz[nMCPar];
	  GenEleCand_phi[nGenEleCand] = MCPar_phi[nMCPar];
	  GenEleCand_charge[nGenEleCand] = (int)(MCPar_pdgid[nMCPar]/11);
	  GenEleCand_eta[nGenEleCand] = MCPar_eta[nMCPar];
	  GenEleCand_p[nGenEleCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+ 
					     (MCPar_py[nMCPar]*MCPar_py[nMCPar])+
					     (MCPar_pz[nMCPar]*MCPar_pz[nMCPar]));
	  GenEleCand_pt[nGenEleCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+ 
					      (MCPar_py[nMCPar]*MCPar_py[nMCPar]));
	  
	  nGenEleCand++;
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
  
  // Two opposite sign dielectrons
  if(nGenEleCand == 2 && (GenEleCand_charge[0]*GenEleCand_charge[1] < 0))
    {
      double mass = pow(GenEleCand_p[0]+GenEleCand_p[1],2);
      mass-=pow(GenEleCand_px[0]+GenEleCand_px[1],2);
      mass-=pow(GenEleCand_py[0]+GenEleCand_py[1],2);
      mass-=pow(GenEleCand_pz[0]+GenEleCand_pz[1],2);
      GenElEl_mass = sqrt(mass);      
      
      double dphi = abs(GenEleCand_phi[0]-GenEleCand_phi[1]);
      if(dphi > 3.14159)
	dphi = (2.0*3.14159) - dphi;       
      GenElEl_dphi = dphi;
    }

  if(passed == true)
    thetree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaEEMC::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaEEMC::endJob() {
  thefile->Write();
  thefile->Close();
}
  
