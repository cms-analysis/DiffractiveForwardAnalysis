
double VertexSeparation(int nvtx, int* vtxtrks, double* vtxz, int* ismumuvtx)
{
        double closestvtx = 9999.;
        if(nvtx == 1)
                return 9999.;
        for(Int_t i = 0; i < nvtx; i++)
        {
                if(ismumuvtx[i] == 0)
                        continue;
                for(Int_t j = 0; j < nvtx && (j!=i); j++)
                {
                        if(fabs(vtxz[i]-vtxz[j]) < fabs(closestvtx))
                                closestvtx = vtxz[i]-vtxz[j];
                }
                if(i==0){
                        for(Int_t j = 1; j < nvtx && (j!=i); j++)
                        {
                                if(fabs(vtxz[i]-vtxz[j]) < fabs(closestvtx))
                                        closestvtx = vtxz[i]-vtxz[j];
                        }
                }
        }
        return closestvtx;
}


bool PassesVertexSelection(int vtxNTrack,double vtxChi2,double vtxNdf,double distance_vertex_z,double vtxZ,int isdimuonvtx)
{
        bool pass = false;
        float vtxseparationzthresh = 0.2;
        float PrimVertexZcut = 24.0;

        if((vtxNTrack>=2)  
     //      && (isdimuonvtx == 1)        
           && (TMath::Prob(vtxChi2,vtxNdf+0.5)>0.001) 
           && (fabs(distance_vertex_z) > vtxseparationzthresh)
           && (fabs(vtxZ)<PrimVertexZcut))
                pass = true;
        return pass;
}



void FillBS_XXX_NUM_XXX()
{
gROOT->Reset();

gStyle->SetPalette(1);

#define pi 3.1413565359
//definition des fichiers + Tree
//  TFile *f0 = new TFile("zerobias_merge1.root"); // 
//  TTree *t0 = f0->Get("ntp1");

  TChain *t0 = new TChain("ntp1");
  t0->Add("/storage/data/cms/store/user/schul/384_InclusiveData2/Inclu2_OniaRunA_*");
//  t0->Add("/storage/data/cms/store/user/schul/384_InclusiveData2/Inclu2_OniaRunB_*");
//  t0->Add("/storage/data/cms/store/user/schul/384_InclusiveData2/Inclu2_MuRunB_*");

// definitions des Trees pour la sauvegarde

// definitions des # d'entrée
  const int NUM0 = t0->GetEntries();

  Int_t techBit0[1][128];
  Int_t var_nTrack0[1];
  Int_t hlt_d0[1];
  TString hlttrigger3 = "HLT_DoubleMu3";
   
Int_t var_nvtx0[1];
Double_t var_vtxZ0[20];
Double_t var_vtxX0[20];
Double_t var_vtxY0[20];
Double_t var_vertexChi2_0[20];
Double_t var_vertexNdf0[20];
Int_t var_vtxTrack0[20];
Int_t var_vtxmumu0[20];
Int_t var_run0[1];

  t0->SetBranchAddress("nPrimVertexCand",var_nvtx0);
  t0->SetBranchAddress("PrimVertexCand_z",var_vtxZ0);
  t0->SetBranchAddress("PrimVertexCand_x",var_vtxX0);
  t0->SetBranchAddress("PrimVertexCand_y",var_vtxY0);
  t0->SetBranchAddress("PrimVertexCand_chi2",var_vertexChi2_0);
  t0->SetBranchAddress("PrimVertexCand_ndof",var_vertexNdf0);
  t0->SetBranchAddress("PrimVertexCand_tracks",var_vtxTrack0);
  t0->SetBranchAddress("PrimVertexCand_mumuTwoTracks",var_vtxmumu0); 
  t0->SetBranchAddress(hlttrigger3,hlt_d0);
  t0->SetBranchAddress("Run",var_run0);


cout<<""<<NUM0<<endl;


  TH1F* nVertex = new TH1F("nVertex","",12,-1.,11.);
  TH1F* VertexX = new TH1F("VertexX","",200,-0.1,0.3);
  TH1F* VertexY = new TH1F("VertexY","",200,-0.2,0.2);
  TH1F* VertexZ = new TH1F("VertexZ","",250,-25.,25.);


  for(Int_t i = 0;i < NUM0;i++){
        t0->GetEntry(i);
	if(var_run0[0]==XXX_NUM_XXX){
          int nPrimVtx = var_nvtx0[0];
          int label_vertex[20];
  	  for(int k=0; k<20; k++){label_vertex[k]=99;}
  	  int nValidVtx(0);

          if(nPrimVtx>=1){
            double distance_vertex_z = VertexSeparation(nPrimVtx,var_vtxTrack0,var_vtxZ0,var_vtxmumu0);
  
            for(Int_t j=0; j<nPrimVtx; j++){
		if(var_vtxTrack0[j]>0 && TMath::Prob(var_vertexChi2_0[j],var_vertexNdf0[j]+0.5)>0.001) nValidVtx++;
                if(PassesVertexSelection(var_vtxTrack0[j],var_vertexChi2_0[j],var_vertexNdf0[j],distance_vertex_z,var_vtxZ0[j],var_vtxmumu0[j]) 
//  		   && sqrt(pow(var_vtxX0[j],2)+pow(var_vtxY0[j],2))<0.2
		   && hlt_d0[0]==1)   
                        {label_vertex[j]=j;}
            }
          }
          for(Int_t j=0; j<nPrimVtx; j++){
	     if(label_vertex[j]!=99){
//		cout<<var_vtxX0[label_vertex[j]]<<"\t"<<var_vtxY0[label_vertex[j]]<<"\t"<<var_vtxZ0[label_vertex[j]]<<endl;
		VertexX->Fill(var_vtxX0[label_vertex[j]]);
                VertexY->Fill(var_vtxY0[label_vertex[j]]);
                VertexZ->Fill(var_vtxZ0[label_vertex[j]]);
	     }			
	  }
	  nVertex->Fill(nValidVtx);
        }
   }


cout<<"Mean x="<<VertexX->GetMean()<<endl;
cout<<"Mean y="<<VertexY->GetMean()<<endl;
cout<<"Mean z="<<VertexZ->GetMean()<<"\t sigma="<<VertexZ->GetRMS()<<endl;

 TFile *thefile1 = new TFile("inclusive_results_XXX_NUM_XXX.root","RECREATE");
  thefile1->cd();
  nVertex->Write();
  VertexX->Write();
  VertexY->Write();
  VertexZ->Write();
}
