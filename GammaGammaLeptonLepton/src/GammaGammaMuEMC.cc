// -*- C++ -*-
//
// Package:    GammaGammaMuEMC
// Class:      GammaGammaMuEMC
// 
/**\class GammaGammaMuEMC GammaGammaMuEMC.cc GammaGammaLeptonLepton/GammaGammaMuEMC/src/GammaGammaMuEMC.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jonathan Hollar
//         Created:  Wed Sep 20 10:08:38 BST 2006
// $Id: GammaGammaMuEMC.cc,v 1.10 2010/02/10 13:49:06 jjhollar Exp $
//
//


#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaMuEMC.h"

#include <FWCore/MessageLogger/interface/MessageLogger.h> 
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h" 
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "FWCore/ParameterSet/interface/ParameterDescriptionNode.h" 

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
GammaGammaMuEMC::GammaGammaMuEMC(const edm::ParameterSet& pset)
{
   //now do what ever initialization is needed
  rootfilename       = pset.getUntrackedParameter<std::string>("outfilename","test.root");
  fillallmc          = pset.getParameter<bool>("FillAllMCParticles");

  //  nEvt=0;
  MUONMAX=10;
  ELEMAX=10;
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
  thetree->Branch("GenMuonCand_status",GenMuonCand_status,"GenMuonCand_status[nGenMuonCand]/I");


  thetree->Branch("nGenEleCand",&nGenEleCand,"nGenEleCand/I"); 
  thetree->Branch("GenEleCand_px",GenEleCand_px,"GenEleCand_px[nGenEleCand]/D"); 
  thetree->Branch("GenEleCand_py",GenEleCand_py,"GenEleCand_py[nGenEleCand]/D"); 
  thetree->Branch("GenEleCand_pz",GenEleCand_pz,"GenEleCand_pz[nGenEleCand]/D"); 
  thetree->Branch("GenEleCand_p",GenEleCand_p,"GenEleCand_p[nGenEleCand]/D"); 
  thetree->Branch("GenEleCand_pt",GenEleCand_pt,"GenEleCand_pt[nGenEleCand]/D"); 
  thetree->Branch("GenEleCand_eta",GenEleCand_eta,"GenEleCand_eta[nGenEleCand]/D"); 
  thetree->Branch("GenEleCand_phi",GenEleCand_phi,"GenEleCand_phi[nGenEleCand]/D"); 
  thetree->Branch("GenEleCand_charge",GenEleCand_charge,"GenEleCand_charge[nGenEleCand]/I"); 
  thetree->Branch("GenEleCand_status",GenEleCand_status,"GenEleCand_status[nGenEleCand]/I");

  thetree->Branch("nGenPhotCand",&nGenPhotCand,"nGenPhotCand/I");  
  thetree->Branch("GenPhotCand_px",GenPhotCand_px,"GenPhotCand_px[nGenPhotCand]/D");  
  thetree->Branch("GenPhotCand_py",GenPhotCand_py,"GenPhotCand_py[nGenPhotCand]/D");  
  thetree->Branch("GenPhotCand_pz",GenPhotCand_pz,"GenPhotCand_pz[nGenPhotCand]/D");  
  thetree->Branch("GenPhotCand_p",GenPhotCand_p,"GenPhotCand_p[nGenPhotCand]/D");  
  thetree->Branch("GenPhotCand_pt",GenPhotCand_pt,"GenPhotCand_pt[nGenPhotCand]/D");  
  thetree->Branch("GenPhotCand_eta",GenPhotCand_eta,"GenPhotCand_eta[nGenPhotCand]/D");  
  thetree->Branch("GenPhotCand_phi",GenPhotCand_phi,"GenPhotCand_phi[nGenPhotCand]/D");  
  thetree->Branch("GenPhotCand_charge",GenPhotCand_charge,"GenPhotCand_charge[nGenPhotCand]/I");  

  thetree->Branch("PassPartonLevel_JetVeto",&PassPartonLevel_JetVeto,"PassPartonLevel_JetVeto/I");

  thetree->Branch("GenMET",&GenMET,"GenMET/D");

  thetree->Branch("nGenJets",&nGenJets,"nGenJets/I");
  thetree->Branch("GenJet_pt",GenJet_pt,"GenJet_pt[nGenJets]/D");
  thetree->Branch("GenJet_eta",GenJet_eta,"GenJet_eta[nGenJets]/D"); 
  thetree->Branch("GenJet_phi",GenJet_phi,"GenJet_phi[nGenJets]/D"); 
  thetree->Branch("GenJet_dR",GenJet_dR,"GenJet_dR[nGenJets]/D");

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

  thetree->Branch("GenGamGam_mass",&GenGamGam_mass,"GenGamGam_mass/D"); 
  thetree->Branch("GenMuE_mass",&GenMuE_mass,"GenMuE_mass/D");
  thetree->Branch("GenMuE_dphi",&GenMuE_dphi,"GenMuE_dphi/D");
  thetree->Branch("GenMuE_pt",&GenMuE_pt,"GenMuE_pt/D"); 

  nEvt = 0;
  nWeird = 0;
}


void 
GammaGammaMuEMC::fillDescriptions(ConfigurationDescriptions & descriptions) { 
   
  descriptions.setComment("Exclusive dimuon GEN-level EDAnalyzer."); 
   
  edm::ParameterSetDescription iDesc;   

  iDesc.add<bool>("FillAllMCParticles", false)->setComment("Set false to fill MC truth information for only muons and protons"); 
  iDesc.addOptionalUntracked<std::string>("outfilename", ("mumu.genlevel.root"))->setComment("output flat ntuple file name");   

  descriptions.add("ParameterDescriptionsForGammaGammaMuEMC", iDesc); 
} 


GammaGammaMuEMC::~GammaGammaMuEMC()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GammaGammaMuEMC::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  nGenMuonCand=0;
  nGenEleCand=0;
  nGenPhotCand=0;
  nGenProtCand=0;
  nGenJets=0;
  nMCPar=0;
  HitInZDC=0;
  HitInCastor=0;
  PassPartonLevel_JetVeto = 1;

  GenGamGam_mass = -1;
  GenMuE_mass = -1;
  GenMuE_dphi = -1;
  GenMuE_pt = -1;
  GenMET = -1;

  bool passed = true;

  nEvt++;
  if((nEvt%10==0 && nEvt<=100)||(nEvt%100==0 && nEvt>100))
    std::cout << "reading event " << nEvt << std::endl;
  
  // step 1: fill some basic MC information into the root tree
  Handle<GenParticleCollection> genParticles; 
  event.getByLabel( "genParticles", genParticles ); 

  for ( size_t i = 0; i < genParticles->size() && i < 5000; ++ i ) 
    {
      const Candidate & p = (*genParticles)[ i ];

       MCPar_status[nMCPar]= p.status();
       //       if(p.status() == 3)
       //	 continue;
       MCPar_px[nMCPar]=p.px();
       MCPar_py[nMCPar]=p.py();
       MCPar_pz[nMCPar]=p.pz();
       MCPar_phi[nMCPar]=atan2(MCPar_py[nMCPar],MCPar_px[nMCPar]);
       MCPar_eta[nMCPar] = p.eta();
       MCPar_pdgid[nMCPar]=p.pdgId();
       MCPar_mass[nMCPar]=p.mass();
       MCPar_e[nMCPar] = sqrt(MCPar_mass[nMCPar]*MCPar_mass[nMCPar] + (MCPar_px[nMCPar]*MCPar_px[nMCPar] + MCPar_py[nMCPar]*MCPar_py[nMCPar] + MCPar_pz[nMCPar]*MCPar_pz[nMCPar])); 

       if(MCPar_status[nMCPar] == 1)
	 {
	   if(fabs(MCPar_pdgid[nMCPar]) < 5)
	     {
	       if(fabs(MCPar_eta[nMCPar]) < 4.7)
		 {
		   if(sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar]) + (MCPar_py[nMCPar]*MCPar_py[nMCPar])) > 30)
		     {
		       PassPartonLevel_JetVeto = 0;
		     }
		 }
	     }
	 }

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
	  GenMuonCand_status[nGenMuonCand] = MCPar_status[nMCPar];
	  
	  nGenMuonCand++;
	}
       if(abs(p.pdgId()) == 11) 
	 { 
	   GenEleCand_px[nGenEleCand] = MCPar_px[nMCPar]; 
	   GenEleCand_py[nGenEleCand] = MCPar_py[nMCPar]; 
	   GenEleCand_pz[nGenEleCand] = MCPar_pz[nMCPar]; 
	   GenEleCand_phi[nGenEleCand] = MCPar_phi[nMCPar]; 
	   GenEleCand_charge[nGenEleCand] = (int)(MCPar_pdgid[nMCPar]/13); 
	   GenEleCand_eta[nGenEleCand] = MCPar_eta[nMCPar]; 
	   GenEleCand_p[nGenEleCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+  
					      (MCPar_py[nMCPar]*MCPar_py[nMCPar])+ 
					      (MCPar_pz[nMCPar]*MCPar_pz[nMCPar])); 
	   GenEleCand_pt[nGenEleCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+  
					       (MCPar_py[nMCPar]*MCPar_py[nMCPar])); 
	   GenEleCand_status[nGenEleCand] = MCPar_status[nMCPar];
           
	   nGenEleCand++; 
	 } 
       if((abs(p.pdgId()) == 24) && (p.status() == 2))  
         {  
           GenPhotCand_px[nGenPhotCand] = MCPar_px[nMCPar];  
           GenPhotCand_py[nGenPhotCand] = MCPar_py[nMCPar];  
           GenPhotCand_pz[nGenPhotCand] = MCPar_pz[nMCPar];  
           GenPhotCand_phi[nGenPhotCand] = MCPar_phi[nMCPar];  
           GenPhotCand_charge[nGenPhotCand] = (int)(MCPar_pdgid[nMCPar]/13);  
           GenPhotCand_eta[nGenPhotCand] = MCPar_eta[nMCPar];  
           GenPhotCand_p[nGenPhotCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+   
					    (MCPar_py[nMCPar]*MCPar_py[nMCPar])+  
					    (MCPar_pz[nMCPar]*MCPar_pz[nMCPar]));  
           GenPhotCand_pt[nGenPhotCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+   
					     (MCPar_py[nMCPar]*MCPar_py[nMCPar]));  
            
           nGenPhotCand++;  
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

  edm::Handle<reco::GenMETCollection>  genmets;
  event.getByLabel( "genMetTrue", genmets );

  if (genmets.isValid()) {
    typedef reco::GenMETCollection::const_iterator gmiter;
    for ( gmiter i=genmets->begin(); i!=genmets->end(); i++) {
      //      mgenmet = i->pt();
      //      mgenphi = i->phi();
      //      mgensum = i->sumEt();
      GenMET = i->pt();
      cout << "MET = " << i->pt() << endl;
    }
  }
  
  edm::Handle<reco::GenJetCollection>  genjets; 
  event.getByLabel( "ak5GenJets", genjets ); 
  cout << "Event!" << endl;
  if (genjets.isValid()) 
    { 
      reco::GenJetCollection mygenjets; 
      mygenjets=*genjets; 
      //      std::sort(mygenjets.begin(),mygenjets.end(),PtGreater()); 
      typedef reco::GenJetCollection::const_iterator gjiter; 
      for ( gjiter i=mygenjets.begin(); i!=mygenjets.end() && nGenJets<1000; i++) 
	{ 
	  GenJet_pt[nGenJets] = i->pt(); 
	  GenJet_phi[nGenJets] = i->phi(); 
	  GenJet_eta[nGenJets] = i->eta(); 

	  double dR1 = sqrt((GenJet_eta[nGenJets]-GenMuonCand_eta[0])*(GenJet_eta[nGenJets]-GenMuonCand_eta[0]) + 
			    (GenJet_phi[nGenJets]-GenMuonCand_phi[0])*(GenJet_phi[nGenJets]-GenMuonCand_phi[0]));
          double dR2 = sqrt((GenJet_eta[nGenJets]-GenEleCand_eta[0])*(GenJet_eta[nGenJets]-GenEleCand_eta[0]) +
                            (GenJet_phi[nGenJets]-GenEleCand_phi[0])*(GenJet_phi[nGenJets]-GenEleCand_phi[0]));
	  double dR = dR1;
	  if(dR2 < dR1)
	    dR = dR2;
	  
	  GenJet_dR[nGenJets] = dR;

	  nGenJets++; 
	  cout << "\t" << nGenJets << endl;
	} 
    }

  // Two opposite sign dimuons
  if(nGenMuonCand >= 1 && nGenEleCand >= 1) 
    { 
      TLorentzVector muvec1; 
      TLorentzVector evec2; 
      TLorentzVector muevec; 
 
      muvec1.SetXYZM(GenMuonCand_px[0],GenMuonCand_py[0],GenMuonCand_pz[0],0.1057); 
      evec2.SetXYZM(GenEleCand_px[0],GenEleCand_py[0],GenEleCand_pz[0],0.000511);  
      muevec = muvec1 + evec2; 
      GenMuE_pt = muevec.Pt(); 
      GenMuE_mass = muevec.M();
      double dphi = abs(GenMuonCand_phi[0]-GenEleCand_phi[0]); 
      if(dphi > 3.14159) 
        dphi = (2.0*3.14159) - dphi;        
      GenMuE_dphi = dphi; 
    } 
  else
    {
      cout << "Found " << nGenMuonCand << " e- and " << nGenEleCand << " e+" << endl; 
      nWeird++;
    }

  TLorentzVector phovec1;  
  TLorentzVector phovec2;  
  TLorentzVector gamgamvec;  
  phovec1.SetXYZM(GenPhotCand_px[0],GenPhotCand_py[0],GenPhotCand_pz[0],80.4030568630);  
  phovec2.SetXYZM(GenPhotCand_px[1],GenPhotCand_py[1],GenPhotCand_pz[1],80.4030568630);   
  gamgamvec = phovec1 + phovec2;  
  GenGamGam_mass = gamgamvec.M(); 

  
  if(passed == true)
    thetree->Fill();

  cout << "N(weird) = " << nWeird << endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaMuEMC::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaMuEMC::endJob() {
  thefile->Write();
  thefile->Close();
}
  
