#define Dimuons2018Macro_cxx
#include "Dimuons2018Macro.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

bool Dimuons2018Macro::FiducalCuts(Float_t trackx210, Float_t tracky210, Float_t trackx220, Float_t tracky220, Int_t arm, Int_t run)
{
  float pixelX0_rotated = 0;
  float pixelY0_rotated = 0;
  float thex220 = 0.0;
  float they220 = 0.0;
  float thex210 = 0.0;
  float they210 = 0.0;


  float xmin_45_210, xmin_45_220, ymin_45_210, ymin_45_220, 
    xmax_45_210, xmax_45_220, ymax_45_210, ymax_45_220,
    xmin_56_210, xmin_56_220, ymin_56_210, ymin_56_220,
    xmax_56_210, xmax_56_220, ymax_56_210, ymax_56_220;
  
  if((run >= 315252) && (run <= 316995))
    {
      // 2018A                                                                                                                                  
      xmin_45_210 = 2.850; xmax_45_210 = 17.927; ymin_45_210 = -11.598; ymax_45_210 = 3.698;
      xmin_45_220 = 2.421; xmax_45_220 = 24.620; ymin_45_220 = -10.898; ymax_45_220 = 4.398;
      xmin_56_210 = 3.275; xmax_56_210 = 18.498; ymin_56_210 = -11.298; ymax_56_210 = 3.298;
      xmin_56_220 = 2.421; xmax_56_220 = 20.045; ymin_56_220 = -10.398; ymax_56_220 = 5.098;
    }
  if((run >= 316998) && (run <= 317696))
    {
      // 2018B1                                                                                                                                 
      xmin_45_210 = 2.850; xmax_45_210 = 17.927; ymin_45_210 = -11.598; ymax_45_210 = 3.698;
      xmin_45_220 = 2.421; xmax_45_220 = 24.620; ymin_45_220 = -10.898; ymax_45_220 = 4.198;
      xmin_56_210 = 3.275; xmax_56_210 = 18.070; ymin_56_210 = -11.198; ymax_56_210 = 4.098;
      xmin_56_220 = 2.564; xmax_56_220 = 20.045; ymin_56_220 = -10.398; ymax_56_220 = 5.098;
    }
  if((run >= 318622) && (run <= 319312))
    {
      // 2018B2                                                                                                                                 
      xmin_45_210 = 2.564; xmax_45_210 = 17.640; ymin_45_210 = -11.598; ymax_45_210 = 4.198;
      xmin_45_220 = 2.140; xmax_45_220 = 24.479; ymin_45_220 = -11.398; ymax_45_220 = 3.798;
      xmin_56_210 = 3.275; xmax_56_210 = 17.931; ymin_56_210 = -10.498; ymax_56_210 = 4.098;
      xmin_56_220 = 2.279; xmax_56_220 = 24.760; ymin_56_220 = -10.598; ymax_56_220 = 4.498;
    }
  if((run >= 319313) && (run <= 320393))
    {
      // 2018C                                                                                                                                  
      xmin_45_210 = 2.564; xmax_45_210 = 17.930; ymin_45_210 = -11.098; ymax_45_210 = 4.198;
      xmin_45_220 = 2.421; xmax_45_220 = 24.620; ymin_45_220 = -11.398; ymax_45_220 = 3.698;
      xmin_56_210 = 3.275; xmax_56_210 = 17.931; ymin_56_210 = -10.498; ymax_56_210 = 4.698;
      xmin_56_220 = 2.279; xmax_56_220 = 24.760; ymin_56_220 = -10.598; ymax_56_220 = 4.398;
    }
  if((run >= 320394) && (run <= 322633))
    {
      // 2018D1                                                                                                                                 
      xmin_45_210 = 2.850; xmax_45_210 = 17.931; ymin_45_210 = -11.098; ymax_45_210 = 4.098;
      xmin_45_220 = 2.421; xmax_45_220 = 24.620; ymin_45_220 = -11.398; ymax_45_220 = 3.698;
      xmin_56_210 = 3.275; xmax_56_210 = 17.931; ymin_56_210 = -10.498; ymax_56_210 = 4.698;
      xmin_56_220 = 2.279; xmax_56_220 = 24.760; ymin_56_220 = -10.598; ymax_56_220 = 4.398;
    }
  if((run >= 323363) && (run <= 325273))
    {
      // 2018D2                                                                                                                                 
      xmin_45_210 = 2.850; xmax_45_210 = 17.931; ymin_45_210 = -10.598; ymax_45_210 = 4.498;
      xmin_45_220 = 2.421; xmax_45_220 = 24.620; ymin_45_220 = -11.698; ymax_45_220 = 3.298;
      xmin_56_210 = 3.275; xmax_56_210 = 17.931; ymin_56_210 =  -9.998; ymax_56_210 = 4.698;
      xmin_56_220 = 2.279; xmax_56_220 = 24.760; ymin_56_220 = -10.598; ymax_56_220 = 3.898;
    }

  pixelX0_rotated = trackx210 * TMath::Cos((-8. / 180.) * TMath::Pi()) -
    tracky210 * TMath::Sin((-8. / 180.) * TMath::Pi());
  pixelY0_rotated = trackx210 * TMath::Sin((-8. / 180.) * TMath::Pi()) +
    tracky210 * TMath::Cos((-8. / 180.) * TMath::Pi());

  thex210 = pixelX0_rotated;
  they210 = pixelY0_rotated;

  thex220 = trackx220;
  they220 = tracky220;
  
  bool pass = false;

  if(arm == 0)
    {
      if(((thex210 >= xmin_45_210) && (thex210 <= xmax_45_210)) && 
	 ((they210 >= ymin_45_210) && (they210 <= ymax_45_210)) &&
	 ((thex220 >= xmin_45_220) && (thex220 <= xmax_45_220)) &&
	 ((they220 >= ymin_45_220) && (they220 <= ymax_45_220)))
	{
	  pass = true;
	}
    }
  if(arm == 1)
    {
      if(((thex210 >= xmin_56_210) && (thex210 <= xmax_56_210)) &&
         ((they210 >= ymin_56_210) && (they210 <= ymax_56_210)) &&
         ((thex220 >= xmin_56_220) && (thex220 <= xmax_56_220)) &&
         ((they220 >= ymin_56_220) && (they220 <= ymax_56_220)))
        {
          pass = true;
        }
    }

  if(pass == false)
    {
      std::cout << "\t\t\tIn pixel fiducial cuts: x, y, x, y = " << thex210 << ", " << they210 << ", " 
		<< thex220 << ", " << they220 << std::endl;
      if(arm == 0)
	std::cout << "\t\t\t\tArm 45 Limits = " << xmin_45_210 << "-" << xmax_45_210 << ", " 
		  << ymin_45_210 << "-" << ymax_45_210 << ", " 
		  << xmin_45_220 << "-" << xmax_45_220 << ", "
		  << ymin_45_220 << "-" << ymax_45_220 << std::endl;
      if(arm == 1)
	std::cout << "\t\t\t\tArm 56 Limits = " << xmin_56_210 << "-" << xmax_56_210 << ", "
		  << ymin_56_210 << "-" << ymax_56_210 << ", "
		  << xmin_56_220 << "-" << xmax_56_220 << ", "
		  << ymin_56_220 << "-" << ymax_56_220 << std::endl;
      
    }

  
  return pass;
}

float Dimuons2018Macro::MultiRPEffCorr(Float_t trackx210, Float_t tracky210, Float_t trackx220, Float_t tracky220, Int_t arm, Int_t run)
{
  float effcorrpixrad210 = 1.0;
  float effcorrpixrad220 = 1.0;
  float effcorrmultitrk = 1.0;

  Int_t theera = 0;

  if((run >= 315252) && (run <= 316995))
    theera = 1; // 2018A                                                                                                                                   
  if((run >= 316998) && (run <= 317696))
    theera = 2; // 2018B1                                                                                                                                  
  if((run >= 318622) && (run <= 319312))
    theera = 3; // 2018B2                                                                                                                                  
  if((run >= 319313) && (run <= 320393))
    theera = 4; // 2018C                                                                                                                                   
  if((run >= 320394) && (run <= 322633))
    theera = 5; // 2018D1                                                                                                                                   
  if((run >= 323363) && (run <= 325273))
    theera = 6; // 2018D2                                                                                                                                  

  // Rad damage
  if(theera == 1)
    {
      if(arm == 0)
        {
          effcorrpixrad220 = hpixeffA45220->GetBinContent(hpixeffA45220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffA45210->GetBinContent(hpixeffA45210->FindBin(trackx210,tracky210));
        }
      if(arm == 1)
        {
          effcorrpixrad220 = hpixeffA56220->GetBinContent(hpixeffA56220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffA56210->GetBinContent(hpixeffA56210->FindBin(trackx210,tracky210));
        }
    }
  if(theera == 2)
    {
      if(arm == 0)
        {
          effcorrpixrad220 = hpixeffB145220->GetBinContent(hpixeffB145220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffB145210->GetBinContent(hpixeffB145210->FindBin(trackx210,tracky210));
        }
      if(arm == 1)
        {
          effcorrpixrad220 = hpixeffB156220->GetBinContent(hpixeffB156220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffB156210->GetBinContent(hpixeffB156210->FindBin(trackx210,tracky210));
        }
    }
  if(theera == 3)
    {
      if(arm == 0)
        {
          effcorrpixrad220 = hpixeffB245220->GetBinContent(hpixeffB245220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffB245210->GetBinContent(hpixeffB245210->FindBin(trackx210,tracky210));
        }
      if(arm == 1)
        {
          effcorrpixrad220 = hpixeffB256220->GetBinContent(hpixeffB256220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffB256210->GetBinContent(hpixeffB256210->FindBin(trackx210,tracky210));
        }
    }
  if(theera == 4)
    {
      if(arm == 0)
        {
          effcorrpixrad220 = hpixeffC45220->GetBinContent(hpixeffC45220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffC45210->GetBinContent(hpixeffC45210->FindBin(trackx210,tracky210));
        }
      if(arm == 1)
        {
          effcorrpixrad220 = hpixeffC56220->GetBinContent(hpixeffC56220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffC56210->GetBinContent(hpixeffC56210->FindBin(trackx210,tracky210));
        }
    }
  if(theera == 5)
    {
      if(arm == 0)
        {
          effcorrpixrad220 = hpixeffD145220->GetBinContent(hpixeffD145220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffD145210->GetBinContent(hpixeffD145210->FindBin(trackx210,tracky210));
        }
      if(arm == 1)
        {
          effcorrpixrad220 = hpixeffD156220->GetBinContent(hpixeffD156220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffD156210->GetBinContent(hpixeffD156210->FindBin(trackx210,tracky210));
        }
    }
  if(theera == 6)
    {
      if(arm == 0)
        {
          effcorrpixrad220 = hpixeffD245220->GetBinContent(hpixeffD245220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffD245210->GetBinContent(hpixeffD245210->FindBin(trackx210,tracky210));
        }
      if(arm == 1)
        {
          effcorrpixrad220 = hpixeffD256220->GetBinContent(hpixeffD256220->FindBin(trackx220,tracky220));
          effcorrpixrad210 = hpixeffD256210->GetBinContent(hpixeffD256210->FindBin(trackx210,tracky210));
        }
    }

  std::cout << "Efficiencies:" << std::endl
	    << "\t x(210) = " << trackx210 << ", y(210) = " << tracky210 << ", x(220) = " << trackx220 << ", y(220) = " << tracky220 << std::endl 
	    << "\tEra = " << theera << ", arm = " << arm << ", run = " << run << ", eff210 = " << effcorrpixrad210 << ", eff220 = " << effcorrpixrad220 << std::endl;

  float efftotal = (effcorrpixrad210 * effcorrpixrad220);
  
  return efftotal;
}

void Dimuons2018Macro::Loop(Int_t multi, Int_t mc, Int_t sb, Int_t yr, Int_t nearfar)
{
//   In a ROOT session, you can do:
//      root> .L Dimuons2018Macro.C
//      root> Dimuons2018Macro t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TH1F *hacop = new TH1F("hacop","hacop",1000,0,0.01);
   TH2F *hcorr45 = new TH2F("hcorr45","hcorr45",500,0,0.25,500,0,0.25);
   TH2F *hcorr56 = new TH2F("hcorr56","hcorr56",500,0,0.25,500,0,0.25);

   TH2F *hcorrmult45 = new TH2F("hcorrmult45","hcorrmult45",500,0,0.25,500,0,0.25);
   TH2F *hcorrmult56 = new TH2F("hcorrmult56","hcorrmult56",500,0,0.25,500,0,0.25);


   TH1F *hres45 = new TH1F("hres45","hres45",300,-5,1);
   TH1F *hres56 = new TH1F("hres56","hres56",300,-5,1);
   TH1F *hressum = new TH1F("hressum","hressum",300,-5,1);

   TH1F *hbin145 = new TH1F("hbin145","hbin145",300,-5,1);
   TH1F *hbin156 = new TH1F("hbin156","hbin156",300,-5,1);
   TH1F *hbin1sum = new TH1F("hbin1sum","hbin1sum",300,-5,1);
   TH1F *hbin245 = new TH1F("hbin245","hbin245",300,-5,1);
   TH1F *hbin256 = new TH1F("hbin256","hbin256",300,-5,1);
   TH1F *hbin2sum = new TH1F("hbin2sum","hbin2sum",300,-5,1);
   TH1F *hbin345 = new TH1F("hbin345","hbin345",300,-5,1);
   TH1F *hbin356 = new TH1F("hbin356","hbin356",300,-5,1);
   TH1F *hbin3sum = new TH1F("hbin3sum","hbin3sum",300,-5,1);
   TH1F *hbin445 = new TH1F("hbin445","hbin445",300,-5,1);
   TH1F *hbin456 = new TH1F("hbin456","hbin456",300,-5,1);
   TH1F *hbin4sum = new TH1F("hbin4sum","hbin4sum",300,-5,1);

   TH1F *hbin3mult45 = new TH1F("hbin3mult45","hbin3mult45",300,-5,1);
   TH1F *hbin3mult56 = new TH1F("hbin3mult56","hbin3mult56",300,-5,1);
   TH1F *hbin4mult45 = new TH1F("hbin4mult45","hbin4mult45",300,-5,1);
   TH1F *hbin4mult56 = new TH1F("hbin4mult56","hbin4mult56",300,-5,1);


   TH1F *hres45mult = new TH1F("hres45mult","hres45mult",300,-5,1);
   TH1F *hres56mult = new TH1F("hres56mult","hres56mult",300,-5,1);
   TH1F *hressummult = new TH1F("hressummult","hressummult",300,-5,1);

   TH1F *hres45multgenmu = new TH1F("hres45multgenmu","hres45multgenmu",300,-5,1);
   TH1F *hres56multgenmu = new TH1F("hres56multgenmu","hres56multgenmu",300,-5,1);
   TH1F *hressummultgenmu = new TH1F("hressummultgenmu","hressummultgenmu",300,-5,1);


   TH1F *hn45220 = new TH1F("hn45220","hn45220",10,0,10);
   TH1F *hn56220 = new TH1F("hn56220","hn56220",10,0,10);
   TH1F *hn45mult = new TH1F("hn45mult","hn45mult",10,0,10);
   TH1F *hn56mult = new TH1F("hn56mult","hn56mult",10,0,10);

   TH1F *hxi45mult = new TH1F("hxi45mult","hxi45mult",100,0,0.25);
   TH1F *hxi56mult = new TH1F("hxi56mult","hxi56mult",100,0,0.25);
   TH1F *hxangle45mult = new TH1F("hxangle45mult","hxangle45mult",50,120,170);
   TH1F *hxangle56mult = new TH1F("hxangle56mult","hxangle56mult",50,120,170);
   TH1F *hxi45single = new TH1F("hxi45single","hxi45single",100,0,0.25);
   TH1F *hxi56single = new TH1F("hxi56single","hxi56single",100,0,0.25);
   TH1F *hxangle45single = new TH1F("hxangle45single","hxangle45single",50,120,170);
   TH1F *hxangle56single = new TH1F("hxangle56single","hxangle56single",50,120,170);

   TH1F *hmresmulti = new TH1F("hmresmulti","hmresmulti",500,-15,5);
   TH1F *hmressingle = new TH1F("hmressingle","hmressingle",500,-15,5);
   TH1F *hmresmixed = new TH1F("hmresmixed","hmresmixed",500,-15,5);
   TH2F *hmcorrmulti = new TH2F("hmcorrmulti","hmcorrmulti",250,0,2500,250,0,2500);
   TH2F *hmcorrmixed = new TH2F("hmcorrmixed","hmcorrmixed",250,0,2500,250,0,2500);
   TH2F *hmcorrsingle = new TH2F("hmcorrsingle","hmcorrsingle",250,0,2500,250,0,2500);
   TH1F *hmmumu = new TH1F("hmmumu","hmmumu",250,0,2500);
   TH2F *hycorrmulti = new TH2F("hycorrmulti","hycorrmulti",250,-2.5,2.5,250,-2.5,2.5);
   TH2F *hycorrsingle = new TH2F("hycorrsingle","hycorrsingle",250,-2.5,2.5,250,-2.5,2.5);
   TH2F *hycorrmixed = new TH2F("hycorrmixed","hycorrmixed",250,-2.5,2.5,250,-2.5,2.5);
   TH2F *hdmdysingle = new TH2F("hdmdysingle","hdmdysingle",1000,-500,500,200,-2,2);

   TH1F *hmresycutsingle = new TH1F("hmresycutsingle","hmresycutsingle",500,-15,5);
   TH1F *hmresycutmulti = new TH1F("hmresycutmulti","hmresycutmulti",500,-15,5);
   TH1F *hmresycutmixed = new TH1F("hmresycutmixed","hmresycutmixed",500,-15,5);

   TH1F *hpzmumumultmatch45 = new TH1F("hpzmumumultmatch45","hpzmumumultmatch45",500,-2000,2000);
   TH1F *hpzmumusinglematch45 = new TH1F("hpzmumusinglematch45","hpzmumusinglematch45",500,-2000,2000);
   TH1F *hpzmumumultmatch56 = new TH1F("hpzmumumultmatch56","hpzmumumultmatch56",500,-2000,2000);
   TH1F *hpzmumusinglematch56 = new TH1F("hpzmumusinglematch56","hpzmumusinglematch56",500,-2000,2000);
   TH1F *hpzmumumultantimatch45 = new TH1F("hpzmumumultantimatch45","hpzmumumultantimatch45",500,-2000,2000);
   TH1F *hpzmumusingleantimatch45 = new TH1F("hpzmumusingleantimatch45","hpzmumusingleantimatch45",500,-2000,2000);
   TH1F *hpzmumumultantimatch56 = new TH1F("hpzmumumultantimatch56","hpzmumumultantimatch56",500,-2000,2000);
   TH1F *hpzmumusingleantimatch56 = new TH1F("hpzmumusingleantimatch56","hpzmumusingleantimatch56",500,-2000,2000);

   TH1F *hgenpzmumu = new TH1F("hgenpzmumu","hgenpzmumu",500,-2000,2000);

   TH2F *hxivst45 = new TH2F("hxivst45","hxivst45",500,0,0.25,100,0,5);
   TH2F *hxivst56 = new TH2F("hxivst56","hxivst56",500,0,0.25,100,0,5);


   // 0.02-0.05, 0.05-0.075, and 0.075-0.25

   ofstream ofs("TextOutputSingleMultiCorr.txt");
   ofstream ofs2("TextOutputMultiRPwithFidEffAndYstar2018.txt");

   ofs2 << "Run,LS,Event,Arm,xing angle,xi(p),xi(mumu),ximatch,t,theta*x,theta*y,y*,eff" << std::endl;

   Int_t nentries = fChain->GetEntries();

   TLorentzVector mu1,mu2,mumu;

   Int_t usemultitracks = 0;
   Int_t ismc = 1;
   Int_t issideband = 0;
   Int_t year = 2017;
   Int_t usenear = 0;

   usemultitracks = multi;
   ismc = mc;
   issideband = sb;
   year = yr;
   usenear = nearfar;

   Int_t id45 = 23;
   Int_t id56 = 123;
   if(usenear == 1)
     {
       id45 = 3;
       id56 = 103;
     }

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(jentry % 10000 == 0)
	std::cout << "Entry " << jentry << "/" << nentries << std::endl;

      Float_t theacop = 1-fabs(Pair_dphi[0])/3.14159;

      if(Pair_mass[0]<110)
	continue;
      if(MuonCand_pt[0]<50 || MuonCand_pt[1]<50)
	continue;
      if(MuonCand_istight[0]!=1 || MuonCand_istight[1]!=1)
	continue;
      if(MuonCand_charge[0] == MuonCand_charge[1])
	continue;
      if((1-fabs(Pair_dphi[0])/3.14159 > 0.009) && (issideband == 0))
	continue;
      //      if((1-fabs(Pair_dphi[0])/3.14159 <= 0.009) && (issideband == 1))
      if((theacop <= 0.009 || theacop > 0.1) && (issideband == 1)) // JH - testing narrower sideband
	continue;
      if((Pair_extratracks0p5mm[0] > 0) && (issideband == 0))
	continue;
      if((Pair_extratracks0p5mm[0] < 5 || Pair_extratracks0p5mm[0] > 10) && (issideband == 1))
	continue;

      mu1.SetPtEtaPhiE(MuonCand_pt[0],MuonCand_eta[0],MuonCand_phi[0],MuonCand_e[0]);
      mu2.SetPtEtaPhiE(MuonCand_pt[1],MuonCand_eta[1],MuonCand_phi[1],MuonCand_e[1]);
      mumu = mu1+mu2;


      hacop->Fill(1-fabs(Pair_dphi[0])/3.14159);

      Float_t mumuxisol1 = (1.0/13000.0)*((MuonCand_pt[0]*TMath::Exp(MuonCand_eta[0]) + MuonCand_pt[1]*TMath::Exp(MuonCand_eta[1])));
      Float_t mumuxisol2 = (1.0/13000.0)*((MuonCand_pt[0]*TMath::Exp(-1*MuonCand_eta[0]) + MuonCand_pt[1]*TMath::Exp(-1*MuonCand_eta[1])));

      Float_t genmumuxisol1 = (1.0/13000.0)*((GenMuonCand_pt[0]*TMath::Exp(GenMuonCand_eta[0]) + GenMuonCand_pt[1]*TMath::Exp(GenMuonCand_eta[1])));
      Float_t genmumuxisol2 = (1.0/13000.0)*((GenMuonCand_pt[0]*TMath::Exp(-1*GenMuonCand_eta[0]) + GenMuonCand_pt[1]*TMath::Exp(-1*GenMuonCand_eta[1])));


      Float_t protxi45 = 0.0;
      Float_t protxi56 = 0.0;
      Float_t protxi45single = 0.0;
      Float_t protxi56single = 0.0;
      Float_t protxi45multi = 0.0;
      Float_t protxi56multi = 0.0;
      Float_t protx45multi210 = 0.0;
      Float_t proty45multi210 = 0.0;
      Float_t protx45multi220 = 0.0;
      Float_t proty45multi220 = 0.0;
      Float_t protx56multi210 = 0.0;
      Float_t proty56multi210 = 0.0;
      Float_t protx56multi220 = 0.0;
      Float_t proty56multi220 = 0.0;
      Float_t mumuxi45 = 0.0;
      Float_t mumuxi56 = 0.0;
      Float_t prott45 = 0.0;
      Float_t prott56 = 0.0;
      Float_t protthx45 = 0.0;
      Float_t protthy45 = 0.0;
      Float_t protthx56 = 0.0;
      Float_t protthy56 = 0.0;
      Float_t protystar45 = 0.0;
      Float_t protystar56 = 0.0;
      Int_t ntrk45 = 0;
      Int_t ntrk56 = 0;

      // JH - do this to select events with only 1 pixel track!
      Int_t ncountpixel45 = 0; 
      Int_t ncountpixel56 = 0;
      Int_t ncountmulti45 = 0;
      Int_t ncountmulti56 = 0;

      for(Int_t p = 0; p < nRecoProtCand; p++)
        {
          if(ProtCand_rpid[p] == id45 && ProtCand_ismultirp[p]==0)
	    ncountpixel45++;
	  if(ProtCand_rpid[p] == id56 && ProtCand_ismultirp[p]==0)
	    ncountpixel56++;
	  if(ProtCand_ismultirp[p]==1 && ProtCand_arm[p]==0)
	    ncountmulti45++;
	  if(ProtCand_ismultirp[p]==1 && ProtCand_arm[p]==1)
	    ncountmulti56++;
	}

      hn45220->Fill(ncountpixel45);
      hn56220->Fill(ncountpixel56);
      hn45mult->Fill(ncountmulti45);
      hn56mult->Fill(ncountmulti56);


      for(Int_t p = 0; p < nRecoProtCand; p++)
	{
	  if(ProtCand_rpid[p] == id45 && ProtCand_ismultirp[p]==0 && ProtCand_trackpixshift1[p] < 1 && ntrk45 < 1 && (ncountpixel45==1 || usemultitracks==1))
	    {
	      protxi45 = ProtCand_xi[p];
	      protxi45single = protxi45;
	      hcorr45->Fill(protxi45,mumuxisol1);
	      hres45->Fill(1 - (protxi45/mumuxisol1));
	      hressum->Fill(1 - (protxi45/mumuxisol1));

              if(mumuxisol1 >= 0.02 && mumuxisol1 < 0.03)
                hbin145->Fill(1 - (protxi45/mumuxisol1));
              if(mumuxisol1 >= 0.03 && mumuxisol1 < 0.04)
                hbin245->Fill(1 - (protxi45/mumuxisol1));
              if(mumuxisol1 >= 0.04 && mumuxisol1 < 0.06)
                hbin345->Fill(1 - (protxi45/mumuxisol1));
	      if(mumuxisol1 >= 0.06)
		hbin445->Fill(1 - (protxi45/mumuxisol1));

              if(fabs(1 - (protxi45/mumuxisol1)) < 0.25)
		{
		  hxi45single->Fill(protxi45);
		  hpzmumusinglematch45->Fill(mumu.Pz());
		}
	      else
		hpzmumusingleantimatch45->Fill(mumu.Pz());

	      std::cout << "Run, Lumi, Event = " << Run << ", " << LumiSection << ", " << EventNum << std::endl;
	      std::cout << "\t xi(RP=3) = " << protxi45 << std::endl;
	      
	      //	      ofs << Run << " " << LumiSection << " " << EventNum << " Single 23 " << protxi45 << " " << (1 - (protxi45/mumuxisol1)) << std::endl;
	      
	      ntrk45++;
	    }
	  if(ProtCand_rpid[p] == id56 && ProtCand_ismultirp[p]==0 && ProtCand_trackpixshift1[p] < 1 && ntrk56 < 1 && (ncountpixel56==1 || usemultitracks==1))
	    {
	      protxi56 = ProtCand_xi[p];
              protxi56single = protxi56;
	      hcorr56->Fill(protxi56,mumuxisol2);
	      hres56->Fill(1 - (protxi56/mumuxisol2));
	      hressum->Fill(1 - (protxi56/mumuxisol2));

	      if(mumuxisol2 >= 0.02 && mumuxisol2 < 0.03)
		hbin156->Fill(1 - (protxi56/mumuxisol2));
              if(mumuxisol2 >= 0.03 && mumuxisol2 < 0.04)
                hbin256->Fill(1 - (protxi56/mumuxisol2));
              if(mumuxisol2 >= 0.04 && mumuxisol2 < 0.06)
                hbin356->Fill(1 - (protxi56/mumuxisol2));
              if(mumuxisol2 >= 0.06)
                hbin456->Fill(1 - (protxi56/mumuxisol2));

              if(fabs(1 - (protxi56/mumuxisol2)) < 0.25)
		{
		  hxi56single->Fill(protxi56);
		  hpzmumusinglematch56->Fill(mumu.Pz());
		}
              else
                hpzmumusingleantimatch56->Fill(mumu.Pz());

	      std::cout << "Run, Lumi, Event = " << Run << " " << LumiSection << " " << EventNum << std::endl;
	      std::cout << "\t xi(RP=103) = " <<protxi56 << std::endl;

	      //	      ofs << Run << " " << LumiSection << " " << EventNum << " Single 123 " << protxi56 << " " << (1 - (protxi56/mumuxisol2)) << std::endl;

	      ntrk56++;
	    }
          if(ProtCand_arm[p] == 0 && ProtCand_trackpixshift1[p] < 1 && ProtCand_ismultirp[p] == 1 && (ncountmulti45==1 || usemultitracks==1))
            {
	      protxi45 = ProtCand_xi[p];
	      prott45 = -1.0*ProtCand_t[p];
	      protthx45 = ProtCand_ThX[p];
	      protthy45 = ProtCand_ThY[p];
              protystar45 = ProtCand_ystar[p];
              protxi45multi = protxi45;

              // JH: shift up 8%                                                                                                            
	      //              protxi45 = protxi45 - (0.08*protxi45);
              // JH: shift up by Jan's systematics file - old version for DPS                                                              
	      //              grsyst45 = (TGraphErrors *)fsyst->Get("2018_postTS2/multi rp-0/g_xi_unc_vs_xi");
	      grsyst45 = (TGraphErrors *)fsyst->Get("2018_postTS2/multi rp-0/xi/g_systematics_vs_xi");
	      std::cout << "Getting syst. shift for proton with xi = " << protxi45 << std::endl;
              float systshift = grsyst45->Eval(protxi45);
	      std::cout << "xi(orig) = " << protxi45 << " + shift of " << systshift << " (" << systshift/protxi45 << \
		"%)" << std::endl;
              protxi45 = protxi45 - systshift;
	      


              // Tracks for eff. correction                                                                                                 
              protx45multi220 = ProtCand_trackx1[p];
              proty45multi220 = ProtCand_tracky1[p];
              protx45multi210 = ProtCand_trackx2[p];
              proty45multi210 = ProtCand_tracky2[p];

	      //              if((FiducalCuts(protx45multi210, proty45multi210, protx45multi220, proty45multi220, 0, Run) == true))
	      if(1)
		{
		  hres45mult->Fill(1 - (protxi45/mumuxisol1));
		  hressummult->Fill(1 - (protxi45/mumuxisol1));
		  hcorrmult45->Fill(protxi45,mumuxisol1);
		  
		  hressummultgenmu->Fill(1 - (protxi45/genmumuxisol1));
		  hres45multgenmu->Fill(1 - (protxi45/genmumuxisol1));
		  
		  if(mumuxisol1 >= 0.04 && mumuxisol1 < 0.06)
		    hbin3mult45->Fill(1 - (protxi45/mumuxisol1));
		  if(mumuxisol1 >= 0.06)
		    hbin4mult45->Fill(1 - (protxi45/mumuxisol1));
		  
		  std::cout << "Run, Lumi, Event = " << Run << ", " << LumiSection << ", " << EventNum << std::endl;
		  std::cout << "\t xi(multi arm=45) = " << protxi45 << std::endl;
		  
		  if(fabs(1 - (protxi45/mumuxisol1)) < 0.10 && (mumuxisol1 >= 0.04))
		    {
		      hxi45mult->Fill(protxi45);
		      hxangle45mult->Fill(CrossingAngle);
		      hpzmumumultmatch45->Fill(mumu.Pz());
		      hxivst45->Fill(protxi45,prott45);
		      
		      float efficiency = MultiRPEffCorr(protx45multi210, proty45multi210, protx45multi220, proty45multi220, 0, Run);
		      
		      ofs2 << Run << "," << LumiSection << "," << EventNum << ",45," << CrossingAngle << "," << protxi45 << "," << mumuxisol1 << ","
			   << 1 - (protxi45/mumuxisol1) << "," << prott45 
                           << ", " << protthx45 << ", " << protthy45 << ", " << ", " << protystar45 << ", " << efficiency
			   << std::endl;
		      
		      if((FiducalCuts(protx45multi210, proty45multi210, protx45multi220, proty45multi220, 0, Run) == false))
			std::cout << "\t\tFailed fiducial cuts: " << protx45multi210 << ", " << proty45multi210 << ", " 
				  << protx45multi220 << ", " << proty45multi220 << std::endl;
		    }
		  else
		    hpzmumumultantimatch45->Fill(mumu.Pz());
		}
            }

	  if(ProtCand_arm[p] == 1 && ProtCand_trackpixshift1[p] < 1 && ProtCand_ismultirp[p] == 1 && (ncountmulti56==1 || usemultitracks==1))
	    {
	      protxi56 = ProtCand_xi[p];
	      prott56 = -1.0*ProtCand_t[p];
              protthx56 = ProtCand_ThX[p];
              protthy56 = ProtCand_ThY[p];
              protystar56 = ProtCand_ystar[p];
              protxi56multi = protxi56;

              // JH: shift up 8%                                                                                                            
              protxi56 = protxi56 - (0.08*protxi56);
              // JH: shift up by Jan's systematics file - old version for DPS                                                              
	      //              grsyst56 = (TGraphErrors *)fsyst->Get("2018_postTS2/multi rp-1/g_xi_unc_vs_xi");
              grsyst56 = (TGraphErrors *)fsyst->Get("2018_postTS2/multi rp-1/xi/g_systematics_vs_xi");
              float systshift = grsyst56->Eval(protxi56);
	      std::cout << "xi(orig) = " << protxi56 << " + shift of " << systshift << " (" << systshift/protxi56 << "%)" << std::endl;
              protxi56 = protxi56 - systshift;


              // Tracks for eff. correction                                                                                                 
              protx56multi220 = ProtCand_trackx1[p];
              proty56multi220 = ProtCand_tracky1[p];
              protx56multi210 = ProtCand_trackx2[p];
              proty56multi210 = ProtCand_tracky2[p];

	      //              if((FiducalCuts(protx56multi210, proty56multi210, protx56multi220, proty56multi220, 1, Run) == true))
	      if(1)
                {
		  hres56mult->Fill(1 - (protxi56/mumuxisol2));
		  hressummult->Fill(1 - (protxi56/mumuxisol2));
		  hcorrmult56->Fill(protxi56,mumuxisol2);
		  
		  hressummultgenmu->Fill(1 - (protxi56/genmumuxisol2));
		  hres56multgenmu->Fill(1 - (protxi56/genmumuxisol2));
		  
		  
		  if(mumuxisol2 >= 0.04 && mumuxisol2 < 0.06)
		    hbin3mult56->Fill(1 - (protxi56/mumuxisol2));
		  if(mumuxisol2 >= 0.06)
		    hbin4mult56->Fill(1 - (protxi56/mumuxisol2));
		  
		  
		  std::cout << "Run, Lumi, Event = " << Run << ", " << LumiSection << ", " << EventNum << std::endl;
		  std::cout << "\t xi(multi arm=56) = " << protxi56 << std::endl;
		  
		  if((fabs(1 - (protxi56/mumuxisol2)) > -0.05) && (fabs(1 - (protxi56/mumuxisol2)) < 0.20) && (mumuxisol2 >= 0.04))
		    {
		      hxi56mult->Fill(protxi56);
		      hxangle56mult->Fill(CrossingAngle);
		      hpzmumumultmatch56->Fill(mumu.Pz());
		      hxivst56->Fill(protxi56,prott56);
		      
		      float efficiency = MultiRPEffCorr(protx56multi210, proty56multi210, protx56multi220, proty56multi220, 1, Run);
		      
		      ofs2 << Run << "," << LumiSection << "," << EventNum << ",56," << CrossingAngle << "," << protxi56 << "," << mumuxisol2 << ","
			   << 1 - (protxi56/mumuxisol2) << "," << prott56 
                           << ", " << protthx56 << ", " << protthy56 << ", " << protystar56 << ", " << efficiency
			   << std::endl;

                      if((FiducalCuts(protx56multi210, proty56multi210, protx56multi220, proty56multi220, 1, Run) == false))
			std::cout << "\t\tFailed fiducial cuts: " << protx56multi210 << ", " << proty56multi210 << ", "
                                  << protx56multi220 << ", " << proty56multi220 << std::endl;
		      
		    }
		  else
		    hpzmumumultantimatch56->Fill(mumu.Pz());
		}
	    }
	}
      
      ofs << Run << " " << LumiSection << " " << EventNum 
	  << " Single 3 " << protxi45single << " " << (1 - (protxi45single/mumuxisol1)) 
	  << " Single 103 " << protxi56single << " " << (1 - (protxi56single/mumuxisol2))
	  << " Multi 3 " << protxi45multi << " " << (1 - (protxi45multi/mumuxisol1))
	  << " Multi 103 " << protxi56multi << " " << (1 - (protxi56multi/mumuxisol2)) << std::endl;

      Float_t mpp = 0.0;
      Float_t ypp = 0.0;
      Float_t mmumu = mumu.M();
      Float_t ymumu = mumu.Rapidity();

      hmmumu->Fill(mmumu);
      if(protxi45multi > 0 && protxi56multi > 0)
	{
	  mpp = 13000.0 * TMath::Sqrt(protxi45multi*protxi56multi);
	  ypp = 0.5 * TMath::Log(protxi45multi/protxi56multi);
	  std::cout << "mpp(multi) = " << mpp << ", mmumu = " << mmumu << std::endl;  
	  hmresmulti->Fill((1 - (mpp/mmumu)));
	  hmcorrmulti->Fill(mmumu,mpp);
	  hycorrmulti->Fill(ymumu,ypp);
          if(fabs(ypp-ymumu) < 0.3)
            {
	      hmresycutmulti->Fill((1 - (mpp/mmumu)));
	      if(fabs(1 - (mpp/mmumu)) < 0.5)
		std::cout << "HEY!!! Look at multi-multi event " << Run << ":" << LumiSection << ":" << EventNum
			  << " with m(mumu) = " << mmumu << ", m(pp) = " << mpp
			  << ", y(mumu) = " << ymumu << ", y(pp) = " << ypp  
			  << ", xi(45) = " << protxi45multi << ", xi(56) = " << protxi56multi << std::endl;
	    }
	}
      if(protxi45multi > 0 && protxi56single > 0)
	{
	  mpp = 13000.0 * TMath::Sqrt(protxi45multi*protxi56single);
          ypp =0.5 * TMath::Log(protxi45multi/protxi56single);
	  std::cout << "mpp(mixed 45multi 56single) = " << mpp << ", mmumu = " << mmumu << std::endl;
          hmresmixed->Fill((1 - (mpp/mmumu)));
          hmcorrmixed->Fill(mmumu,mpp);
	  hycorrmixed->Fill(ymumu,ypp);
          if(fabs(ypp-ymumu) < 0.3)
	    {
	      hmresycutmixed->Fill((1 - (mpp/mmumu)));
              if(fabs(1 - (mpp/mmumu)) < 0.5)
		std::cout << "HEY!!! Look at multi-sinlge event " << Run << ":" << LumiSection << ":" << EventNum
			  << " with m(mumu) = " << mmumu << ", m(pp) = " << mpp
			  << ", y(mumu) = " << ymumu << ", y(pp) = " << ypp 
                          << ", xi(45) = " << protxi45multi << ", xi(56) = " << protxi56single << std::endl;
	    }
	}
      if(protxi56multi > 0 && protxi45single > 0)
	{
	  mpp = 13000.0 * TMath::Sqrt(protxi45single*protxi56multi);
          ypp =0.5 * TMath::Log(protxi45single/protxi56multi);
	  std::cout << "mpp(mixed 45single 56multi) = " << mpp << ", mmumu = " << mmumu << std::endl;
          hmresmixed->Fill((1 - (mpp/mmumu)));
          hmcorrmixed->Fill(mmumu,mpp);
	  hycorrmixed->Fill(ymumu,ypp);
          if(fabs(ypp-ymumu) < 0.3)
	    {
	      hmresycutmixed->Fill((1 - (mpp/mmumu)));
              if(fabs(1 - (mpp/mmumu)) < 0.5)
		std::cout << "HEY!!! Look at single-multi event " << Run << ":" << LumiSection << ":" << EventNum
			  << " with m(mumu) = " << mmumu << ", m(pp) = " << mpp
			  << ", y(mumu) = " << ymumu << ", y(pp) = " << ypp 
			  << ", xi(45) = " << protxi45single << ", xi(56) = " << protxi56multi << std::endl;
	    }
	}
      if(protxi45single > 0 && protxi56single > 0)
	{
	  mpp = 13000.0 * TMath::Sqrt(protxi45single*protxi56single);
          ypp =0.5 * TMath::Log(protxi45single/protxi56single);
	  std::cout << "mpp(single) = " << mpp << ", mmumu = " << mmumu << std::endl;
	  std::cout << "ypp(single) = " << ypp << ", ymumu = " << ymumu << std::endl;
	  hmressingle->Fill((1 - (mpp/mmumu)));
          hmcorrsingle->Fill(mmumu,mpp);
	  hycorrsingle->Fill(ymumu,ypp);
	  hdmdysingle->Fill(mpp-mmumu,ypp-ymumu);
	  if(fabs(ypp-ymumu) < 0.3)
	    {
	      hmresycutsingle->Fill((1 - (mpp/mmumu)));
	      if(fabs(1 - (mpp/mmumu)) < 0.5)
		std::cout << "HEY!!! Look at single-single event " << Run << ":" << LumiSection << ":" << EventNum 
			  << " with m(mumu) = " << mmumu << ", m(pp) = " << mpp 
			  << ", y(mumu) = " << ymumu << ", y(pp) = " << ypp 
                          << ", xi(45) = " << protxi45single << ", xi(56) = " << protxi56single << std::endl;
	    }
	}
   }

   ofs.close();
   ofs2.close();
          
   //   TFile *fx = new TFile("MoreDimuons2017BCDEFWithMultiLegacyFinalFromDB.root","RECREATE");
   //   TFile *fx = new TFile("MoreDimuons2017BCDEFSingleTrackiNearLegacyFinalFromDB.root","RECREATE");
   //   TFile *fx = new TFile("MoreDimuons2017BCDEFSingleTrackFarLegacyFinalFromDBWithMass.root","RECREATE");                                          

   //   TFile *fx = new TFile("MoreDimuons2017BCDEFWithMultiTrackFarLegacyFinalFromDBWithMass.root","RECREATE");
   //   TFile *fx = new TFile("MoreDimuons2017MCallanglesSingleTrackLegacyFinalFromDB.root","RECREATE");                          
   //   TFile *fx = new TFile("MoreDimuons2017MCallanglesWithMultiLegacyFinalFromDBWithMass.root","RECREATE");
   
   TFile *fx;
   
   if(issideband == 0)
     {
       if(usemultitracks == 1)
	 {
	   if(ismc == 0)
	     {
	       if(year == 2017)
		 fx = new TFile("MoreDimuons2017BCDEFWithMultiTrackFarLegacyFinalFromDBWithMassMoreBins.root","RECREATE");
	       if(year == 2018)
		 {
		   if(usenear == 0){fx = new TFile("MoreDimuons2018BCDWithMultiTrackFarLegacyFinalFromDBWithMassMoreBins.root","RECREATE");}
                   if(usenear == 1){fx = new TFile("MoreDimuons2018BCDWithMultiTrackNearLegacyFinalFromDBWithMassMoreBins.root","RECREATE");}
		 }
	     }
	   if(ismc == 1)
	     {
	       if(year == 2017)
		 fx = new TFile("MoreDimuons2017MCallanglesWithMultiLegacyFinalFromDBWithMassMoreBins.root","RECREATE");
	       if(year == 2018)
		 fx = new TFile("MoreDimuons2018MCallanglesWithMultiLegacyFinalFromDBWithMassMoreBins.root","RECREATE");
	     }
	 }
       if(usemultitracks == 0)
	 {
	   if(ismc == 0)
	     {
	       if(year == 2017)
		 {
		   if(usenear == 0){fx = new TFile("MoreDimuons2017BCDEFSingleTrackPixelsLegacyFinalFromDBWithMass.root","RECREATE");}
		   if(usenear == 1){fx = new TFile("MoreDimuons2017BCDEFSingleTrackStripsLegacyFinalFromDBWithMass.root","RECREATE");}
		   //	     fx = new TFile("MoreDimuons2017BCDEFSingleTrackStripsLegacyFinalFromDBWithMass.root","RECREATE");
		 }
		   if(year == 2018)
		 {
		   if(usenear == 0){fx = new TFile("MoreDimuons2018BCDSingleTrackPixelsFarLegacyFinalFromDBWithMass.root","RECREATE");}
		   if(usenear == 1){fx = new TFile("MoreDimuons2018BCDSingleTrackPixelsNearLegacyFinalFromDBWithMass.root","RECREATE");}
		 }
	     }
	   if(ismc == 1)
	     {
	       if(year == 2017)
		 fx = new TFile("MoreDimuons2017MCallanglesWithMultiLegacyFinalFromDBWithMass.root","RECREATE");
	       if(year == 2018)
		 fx = new TFile("MoreDimuons2018MCallanglesWithMultiLegacyFinalFromDBWithMass.root","RECREATE");
	     }
	 }
     }
   if(issideband == 1)
     {
       if(usemultitracks == 1)
	 if(ismc == 0)
	   {
	     if(year == 2017)
	       fx = new TFile("MoreDimuons2017BCDEFWithMultiTrackFarLegacyFinalFromDBWithMassAntiAcop1to5tracksMoreBins.root","RECREATE");
             if(year == 2018)
	       {
		 if(usenear == 0){fx = new TFile("MoreDimuons2018BCDWithMultiTrackFarLegacyFinalFromDBWithMassAntiAcop1to5tracksMoreBins.root","RECREATE");}
                 if(usenear == 1){fx = new TFile("MoreDimuons2018BCDWithMultiTrackNearLegacyFinalFromDBWithMassAntiAcop1to5tracksMoreBins.root","RECREATE");}
	       }
	   }
       if(usemultitracks == 0)
	 if(ismc == 0)
	   {
	     if(year == 2017)
	       {
		 if(usenear == 0){fx = new TFile("MoreDimuons2017BCDEFSingleTrackFarLegacyFinalFromDBWithMassAntiAcop1to5tracks.root","RECREATE");}
                 if(usenear == 1){fx = new TFile("MoreDimuons2017BCDEFSingleTrackNearLegacyFinalFromDBWithMassAntiAcop1to5tracks.root","RECREATE");}
	       }
             if(year == 2018)
	       {
		 if(usenear == 0){fx = new TFile("MoreDimuons2018BCDSingleTrackFarLegacyFinalFromDBWithMassAntiAcop1to5tracks.root","RECREATE");}
                 if(usenear == 1){fx = new TFile("MoreDimuons2018BCDSingleTrackNearLegacyFinalFromDBWithMassAntiAcop1to5tracks.root","RECREATE");}
	       }
	   }
     }
   
   hacop->Write();
   hcorr45->Write();
   hres45->Write();
   hcorr56->Write();
   hres56->Write();
   hressum->Write();
   hres45mult->Write();
   hres56mult->Write();
   hressummult->Write();
   
   hressummultgenmu->Write();
   hres45multgenmu->Write();
   hres56multgenmu->Write();


   hbin145->Write();
   hbin245->Write();
   hbin345->Write();
   hbin445->Write();

   hbin156->Write();
   hbin256->Write();
   hbin356->Write();
   hbin456->Write();

   hbin3mult45->Write();
   hbin3mult56->Write();
   hbin4mult45->Write();
   hbin4mult56->Write();

   hn45220->Write();
   hn56220->Write();
   hn45mult->Write();
   hn56mult->Write();

   hcorrmult45->Write();
   hcorrmult56->Write();

   hxi45mult->Write();
   hxi56mult->Write();
   hxangle45mult->Write();
   hxangle56mult->Write();
   hxi45single->Write();
   hxi56single->Write();

   hmresmulti->Write();
   hmresmixed->Write();
   hmressingle->Write();
   hmcorrsingle->Write();
   hmcorrmulti->Write();
   hmcorrmixed->Write();

   hmmumu->Write();
   hycorrsingle->Write();
   hycorrmulti->Write();
   hycorrmixed->Write();
   hdmdysingle->Write();
   hmresycutmulti->Write();
   hmresycutmixed->Write();
   hmresycutsingle->Write();
   
   hpzmumusinglematch45->Write();
   hpzmumusinglematch56->Write();
   hpzmumumultmatch45->Write();
   hpzmumumultmatch56->Write();
   hpzmumusingleantimatch45->Write();
   hpzmumusingleantimatch56->Write();
   hpzmumumultantimatch45->Write();
   hpzmumumultantimatch56->Write();

   hxivst45->Write();
   hxivst56->Write();

   fx->Write();
}
