void setJHStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR"); 

  // For the canvas: 
  tdrStyle->SetCanvasBorderMode(0); 
  tdrStyle->SetCanvasColor(kWhite); 
  tdrStyle->SetCanvasDefH(600); //Height of canvas 
  tdrStyle->SetCanvasDefW(600); //Width of canvas 
  tdrStyle->SetCanvasDefX(0);   //POsition on screen 
  tdrStyle->SetCanvasDefY(0); 
 
  // For the Pad: 
  tdrStyle->SetPadBorderMode(0); 
  // tdrStyle->SetPadBorderSize(Width_t size = 1); 
  tdrStyle->SetPadColor(kWhite); 
  tdrStyle->SetPadGridX(false); 
  tdrStyle->SetPadGridY(false); 
  tdrStyle->SetGridColor(0); 
  tdrStyle->SetGridStyle(3); 
  tdrStyle->SetGridWidth(1); 

  // For the frame: 
  tdrStyle->SetFrameBorderMode(0); 
  tdrStyle->SetFrameBorderSize(1); 
  tdrStyle->SetFrameFillColor(0); 
  tdrStyle->SetFrameFillStyle(0); 
  tdrStyle->SetFrameLineColor(1); 
  tdrStyle->SetFrameLineStyle(1); 
  tdrStyle->SetFrameLineWidth(1); 

  // For the histo: 
  tdrStyle->SetHistLineColor(1); 
  //  tdrStyle->SetHistLineStyle(0); 
  tdrStyle->SetHistLineWidth(1); 
 
  //  tdrStyle->SetEndErrorSize(2); 
  //  tdrStyle->SetErrorX(0.); 
   
  tdrStyle->SetMarkerStyle(20); 


  //For the fit/function: 
  tdrStyle->SetOptFit(1); 
  tdrStyle->SetFitFormat("5.4g"); 
  tdrStyle->SetFuncColor(2); 
  tdrStyle->SetFuncStyle(1); 
  tdrStyle->SetFuncWidth(1); 
 
  //For the date: 
  tdrStyle->SetOptDate(0); 
  // tdrStyle->SetDateX(Float_t x = 0.01); 
  // tdrStyle->SetDateY(Float_t y = 0.01); 

  // For the statistics box: 
  tdrStyle->SetOptFile(0); 
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr"); 
  tdrStyle->SetStatColor(kWhite); 
  tdrStyle->SetStatFont(42); 
  tdrStyle->SetStatFontSize(0.025); 
  tdrStyle->SetStatTextColor(1); 
  tdrStyle->SetStatFormat("6.4g"); 
  tdrStyle->SetStatBorderSize(1); 
  tdrStyle->SetStatH(0.1); 
  tdrStyle->SetStatW(0.15); 
  // tdrStyle->SetStatStyle(Style_t style = 1001); 
  // tdrStyle->SetStatX(Float_t x = 0); 
  // tdrStyle->SetStatY(Float_t y = 0); 
 
  // Margins: 
  tdrStyle->SetPadTopMargin(0.05); 
  tdrStyle->SetPadBottomMargin(0.13); 
  tdrStyle->SetPadLeftMargin(0.13); 
  tdrStyle->SetPadRightMargin(0.05); 

  // For the Global title: 
 
  //  tdrStyle->SetOptTitle(0); 
  tdrStyle->SetTitleFont(42); 
  tdrStyle->SetTitleColor(1); 
  tdrStyle->SetTitleTextColor(1); 
  tdrStyle->SetTitleFillColor(10); 
  tdrStyle->SetTitleFontSize(0.05); 
  // tdrStyle->SetTitleH(0); // Set the height of the title box 
  // tdrStyle->SetTitleW(0); // Set the width of the title box 
  // tdrStyle->SetTitleX(0); // Set the position of the title box 
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box 
  // tdrStyle->SetTitleStyle(Style_t style = 1001); 
  // tdrStyle->SetTitleBorderSize(2); 
 
  // For the axis titles: 
 
  tdrStyle->SetTitleColor(1, "XYZ"); 
  tdrStyle->SetTitleFont(42, "XYZ"); 
  tdrStyle->SetTitleSize(0.06, "XYZ"); 
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size? 
  // tdrStyle->SetTitleYSize(Float_t size = 0.02); 
  tdrStyle->SetTitleXOffset(0.9); 
  tdrStyle->SetTitleYOffset(1.05); 
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset 
 
  // For the axis labels: 
 
  tdrStyle->SetLabelColor(1, "XYZ"); 
  tdrStyle->SetLabelFont(42, "XYZ"); 
  tdrStyle->SetLabelOffset(0.007, "XYZ"); 
  tdrStyle->SetLabelSize(0.05, "XYZ"); 

  // For the axis: 
 
  tdrStyle->SetAxisColor(1, "XYZ"); 
  tdrStyle->SetStripDecimals(kTRUE); 
  tdrStyle->SetTickLength(0.03, "XYZ"); 
  tdrStyle->SetNdivisions(510, "XYZ"); 
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame 
  tdrStyle->SetPadTickY(1); 
 
  // Change for log plots: 
  tdrStyle->SetOptLogx(0); 
  tdrStyle->SetOptLogy(0); 
  tdrStyle->SetOptLogz(0); 
 
  tdrStyle->SetPalette(51);

  // Postscript options: 
  // tdrStyle->SetPaperSize(15.,15.); 
  // tdrStyle->SetLineScalePS(Float_t scale = 3); 
  // tdrStyle->SetLineStyleString(Int_t i, const char* text); 
  // tdrStyle->SetHeaderPS(const char* header); 
  // tdrStyle->SetTitlePS(const char* pstitle); 
 
  // tdrStyle->SetBarOffset(Float_t baroff = 0.5); 
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5); 
  // tdrStyle->SetPaintTextFormat(const char* format = "g"); 
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0); 
  // tdrStyle->SetTimeOffset(Double_t toffset); 
  // tdrStyle->SetHistMinimumZero(kTRUE); 
 
  // personnal additions 
  tdrStyle->cd(); 
 
} 

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.13);
  tdrStyle->SetPadRightMargin(0.05);

// For the Global title:

//  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

// personnal additions
  tdrStyle->cd();

}

void PlotSubset()
{
  PlotSigVsBkgWW(1,3,1,1,"nextratracks_pt30_mcatnlo");  
  PlotSigVsBkgWW(1,3,2,1,"nextratracks_pt30_mcatnlo_log");   
}

void PlotSigVsBkgWWAll()
{
  PlotSigVsBkgWW(1,1,1,1,"nextratracks");
  PlotSigVsBkgWW(1,1,2,1,"nextratracks_log"); 
  PlotSigVsBkgWW(1,3,1,1,"nextratracks_pt30"); 
  PlotSigVsBkgWW(1,3,2,1,"nextratracks_pt30_log");  
  PlotSigVsBkgWW(1,4,1,1,"nextratracks_ptlessthan30");  
  PlotSigVsBkgWW(1,4,2,1,"nextratracks_ptlessthan30_log");   

  PlotSigVsBkgWW(2,1,1,1,"memu_1to6tracks"); 
  PlotSigVsBkgWW(2,1,2,1,"memu_1to6tracks_log");
  PlotSigVsBkgWW(2,3,1,1,"memu_1to6tracks_pt30");  
  PlotSigVsBkgWW(2,3,2,1,"memu_1to6tracks_pt30_log"); 
  PlotSigVsBkgWW(2,4,1,1,"memu_1to6tracks_ptlessthan30");   
  PlotSigVsBkgWW(2,4,2,1,"memu_1to6tracks_ptlessthan30_log");  
  PlotSigVsBkgWW(2,777,1,1,"memu_0tracks_ptlessthan30");    
  PlotSigVsBkgWW(2,888,1,1,"memu_0tracks_pt30_unblind");

  PlotSigVsBkgWW(3,1,1,1,"dphiemu_1to6tracks");  
  PlotSigVsBkgWW(3,1,2,1,"dphiemu_1to6tracks_log"); 
  PlotSigVsBkgWW(3,3,1,1,"dphiemu_1to6tracks_pt30");   
  PlotSigVsBkgWW(3,3,2,1,"dphiemu_1to6tracks_pt30_log");  
  PlotSigVsBkgWW(3,4,1,1,"dphiemu_1to6tracks_ptlessthan30");    
  PlotSigVsBkgWW(3,4,2,1,"dphiemu_1to6tracks_ptlessthan30_log");   
  PlotSigVsBkgWW(3,777,1,1,"dphiemu_0tracks_ptlessthan30");     
  PlotSigVsBkgWW(3,888,1,1,"dphiemu_0tracks_pt30_unblind"); 

  PlotSigVsBkgWW(4,1,1,1,"dptemu_1to6tracks");   
  PlotSigVsBkgWW(4,1,2,1,"dptemu_1to6tracks_log");  
  PlotSigVsBkgWW(4,3,1,1,"dptemu_1to6tracks_pt30");    
  PlotSigVsBkgWW(4,3,2,1,"dptemu_1to6tracks_pt30_log");   
  PlotSigVsBkgWW(4,4,1,1,"dptemu_1to6tracks_ptlessthan30");     
  PlotSigVsBkgWW(4,4,2,1,"dptemu_1to6tracks_ptlessthan30_log");    
  //  PlotSigVsBkgWW(4,777,1,1,"dptemu_0tracks_ptlessthan30");
  PlotSigVsBkgWW(4,888,1,1,"dptemu_0tracks_pt30_unblind");  


  PlotSigVsBkgWW(5,1,1,1,"ptemu_1to6tracks");    
  PlotSigVsBkgWW(5,1,2,1,"ptemu_1to6tracks_log");   
  PlotSigVsBkgWW(5,999,1,1,"ptemu_0tracks_pt30_unblindAQGC");   


  //  PlotSigVsBkgWW(5,888,1,1,"ptemu_0tracks_pt30_unblind");  
  PlotSigVsBkgWW(6,1,1,1,"pfmet_1to6tracks");    
  PlotSigVsBkgWW(6,1,2,1,"pfmet_1to6tracks_log");   
  PlotSigVsBkgWW(6,3,1,1,"pfmet_1to6tracks_pt30");     
  PlotSigVsBkgWW(6,3,2,1,"pfmet_1to6tracks_pt30_log");    
  PlotSigVsBkgWW(6,4,1,1,"pfmet_1to6tracks_ptlessthan30");      
  PlotSigVsBkgWW(6,4,2,1,"pfmet_1to6tracks_ptlessthan30_log");     
  //  PlotSigVsBkgWW(6,777,1,1,"pfmet_0tracks_ptlessthan30");
  PlotSigVsBkgWW(6,888,1,1,"pfmet_0tracks_pt30_unblind");  


  PlotSigVsBkgWW(7,1,1,1,"npv_1to6tracks");     
  PlotSigVsBkgWW(7,1,2,1,"npv_1to6tracks_log");    
  /*
  PlotSigVsBkgWW(7,3,1,1,"npv_1to6tracks_pt30");      
  PlotSigVsBkgWW(7,3,2,1,"npv_1to6tracks_pt30_log");     
  PlotSigVsBkgWW(7,4,1,1,"npv_1to6tracks_ptlessthan30");       
  PlotSigVsBkgWW(7,4,2,1,"npv_1to6tracks_ptlessthan30_log");      
  */
  PlotSigVsBkgWW(7,888,1,1,"npv_0tracks_pt30_unblind");  

  PlotSigVsBkgWW(12,1,1,1,"ptmu_1to6tracks");      
  PlotSigVsBkgWW(12,1,2,1,"ptmu_1to6tracks_log");     
  PlotSigVsBkgWW(12,3,1,1,"ptmu_1to6tracks_pt30");       
  PlotSigVsBkgWW(12,3,2,1,"ptmu_1to6tracks_pt30_log");      
  /*
  PlotSigVsBkgWW(12,4,1,1,"ptmu_1to6tracks_ptlessthan30");        
  PlotSigVsBkgWW(12,4,2,1,"ptmu_1to6tracks_ptlessthan30_log");       
  */
  PlotSigVsBkgWW(12,888,1,1,"ptmu_0tracks_pt30_unblind");  

  PlotSigVsBkgWW(13,1,1,1,"etamu_1to6tracks");       
  PlotSigVsBkgWW(13,1,2,1,"etamu_1to6tracks_log");      
  PlotSigVsBkgWW(13,3,1,1,"etamu_1to6tracks_pt30");        
  PlotSigVsBkgWW(13,3,2,1,"etamu_1to6tracks_pt30_log");       
  /*
  PlotSigVsBkgWW(13,4,1,1,"etamu_1to6tracks_ptlessthan30");         
  PlotSigVsBkgWW(13,4,2,1,"etamu_1to6tracks_ptlessthan30_log");        
  */
  PlotSigVsBkgWW(13,888,1,1,"etamu_0tracks_pt30_unblind");  

  PlotSigVsBkgWW(14,1,1,1,"etele_1to6tracks");       
  PlotSigVsBkgWW(14,1,2,1,"etele_1to6tracks_log");      
  PlotSigVsBkgWW(14,3,1,1,"etele_1to6tracks_pt30");        
  PlotSigVsBkgWW(14,3,2,1,"etele_1to6tracks_pt30_log");       
  /*
  PlotSigVsBkgWW(14,4,1,1,"etele_1to6tracks_ptlessthan30");         
  PlotSigVsBkgWW(14,4,2,1,"etele_1to6tracks_ptlessthan30_log");        
  */
  PlotSigVsBkgWW(14,888,1,1,"elele_0tracks_pt30_unblind");  

  PlotSigVsBkgWW(15,1,1,1,"etaele_1to6tracks");       
  PlotSigVsBkgWW(15,1,2,1,"etaele_1to6tracks_log");      
  PlotSigVsBkgWW(15,3,1,1,"etaele_1to6tracks_pt30");        
  PlotSigVsBkgWW(15,3,2,1,"etaele_1to6tracks_pt30_log");       
  /*
  PlotSigVsBkgWW(15,4,1,1,"etaele_1to6tracks_ptlessthan30");         
  PlotSigVsBkgWW(15,4,2,1,"etaele_1to6tracks_ptlessthan30_log");        
  */
  PlotSigVsBkgWW(15,888,1,1,"etaele_0tracks_pt30_unblind");  

  PlotSigVsBkgWW(16,1,1,1,"ptextratracks_1to6tracks");        
  PlotSigVsBkgWW(16,1,2,1,"ptextratracks_1to6tracks_log");       
  PlotSigVsBkgWW(16,3,1,1,"ptextratracks_1to6tracks_pt30");         
  PlotSigVsBkgWW(16,3,2,1,"ptextratracks_1to6tracks_pt30_log");        
  /*
  PlotSigVsBkgWW(16,4,1,1,"ptextratracks_1to6tracks_ptlessthan30");          
  PlotSigVsBkgWW(16,4,2,1,"ptextratracks_1to6tracks_ptlessthan30_log");        
  */

  PlotSigVsBkgWW(17,1,1,1,"etaextratracks_1to6tracks");         
  PlotSigVsBkgWW(17,1,2,1,"etaextratracks_1to6tracks_log");        
  PlotSigVsBkgWW(17,3,1,1,"etaextratracks_1to6tracks_pt30");          
  PlotSigVsBkgWW(17,3,2,1,"etaextratracks_1to6tracks_pt30_log");         
  /*
  PlotSigVsBkgWW(17,4,1,1,"etaextratracks_1to6tracks_ptlessthan30");           
  PlotSigVsBkgWW(17,4,2,1,"etaextratracks_1to6tracks_ptlessthan30_log");     
  */

  PlotSigVsBkgWW(18,1,1,1,"sumptextratracks_1to6tracks");         
  PlotSigVsBkgWW(18,1,2,1,"sumptextratracks_1to6tracks_log");        
  PlotSigVsBkgWW(18,3,1,1,"sumptextratracks_1to6tracks_pt30");          
  PlotSigVsBkgWW(18,3,2,1,"sumptextratracks_1to6tracks_pt30_log");         
  /*
  PlotSigVsBkgWW(18,4,1,1,"sumptextratracks_1to6tracks_ptlessthan30");           
  PlotSigVsBkgWW(18,4,2,1,"sumptextratracks_1to6tracks_ptlessthan30_log");     
  */  
}

void PlotSigVsBkgWW(Int_t thevar = 1, Int_t cutset = 1, Int_t mode = 1, Int_t save = 0, TString savename="none")
{
  // setTDRStyle();
  setJHStyle();

 ifstream ifs2011Amuid("TIGHT_nL8_2011A_ratio.txt");  
 ifstream ifs2011Bmuid("TIGHT_nL8_2011B_ratio.txt");  
 ifstream ifs2011mu17trg("Mu17_2011_ratio.txt");
 ifstream ifs2011mu8trg("Mu8_2011_ratio.txt");
 ifstream ifs2011eleid("MEDIUM_Ele_2011_ratio.txt");

 Double_t muoneffs[28][8]; 
 Double_t muonhlt8effs[17][8];
 Double_t muonhlt17effs[17][8];
 Double_t eleeffs[30][8];
 
 Int_t iter = 0;  
 Double_t etalo, etahi, ptlo, pthi, eff, errlo, errhi;  

 while(!ifs2011eleid.eof() && iter<30)
   {
     ifs2011eleid >> etalo >> etahi >> ptlo >> pthi >> eff >> errlo >> errhi;
     eleeffs[iter][0] = etalo;  
     eleeffs[iter][1] = etahi;   
     eleeffs[iter][2] = ptlo;   
     eleeffs[iter][3] = pthi;   
     eleeffs[iter][4] = eff;   
     eleeffs[iter][5] = errlo;   
     eleeffs[iter][6] = errhi;   
     eleeffs[iter][7] = 1;   
     iter++; 
   }  
 ifs2011eleid.close();   

 iter=0;
 while(!ifs2011mu17trg.eof())
   {
     ifs2011mu17trg >> etalo >> etahi >> ptlo >> pthi >> eff >> errlo >> errhi;
     muonhlt17effs[iter][0] = etalo; 
     muonhlt17effs[iter][1] = etahi;  
     muonhlt17effs[iter][2] = ptlo;  
     muonhlt17effs[iter][3] = pthi;   
     muonhlt17effs[iter][4] = eff;    
     muonhlt17effs[iter][5] = errlo;     
     muonhlt17effs[iter][6] = errhi;     
     muonhlt17effs[iter][7] = 1;     
     iter++;
   }
 ifs2011mu17trg.close();

 iter = 0;
 while(!ifs2011mu8trg.eof()) 
   { 
     ifs2011mu8trg >> etalo >> etahi >> ptlo >> pthi >> eff >> errlo >> errhi; 
     muonhlt8effs[iter][0] = etalo;  
     muonhlt8effs[iter][1] = etahi;   
     muonhlt8effs[iter][2] = ptlo;   
     muonhlt8effs[iter][3] = pthi;    
     muonhlt8effs[iter][4] = eff;     
     muonhlt8effs[iter][5] = errlo;      
     muonhlt8effs[iter][6] = errhi;      
     muonhlt8effs[iter][7] = 1;      
     iter++; 
   } 
 ifs2011mu8trg.close(); 

 iter = 0;
 while(!ifs2011Amuid.eof())  
   {  
     ifs2011Amuid >> etalo >> etahi >> ptlo >> pthi >> eff >> errlo >> errhi;  
     muoneffs[iter][0] = etalo; 
     muoneffs[iter][1] = etahi;  
     muoneffs[iter][2] = ptlo;  
     muoneffs[iter][3] = pthi;  
     muoneffs[iter][4] = eff;  
     muoneffs[iter][5] = errlo;  
     muoneffs[iter][6] = errhi;  
     muoneffs[iter][7] = 1;  
     iter++;
   } 
 ifs2011Amuid.close();  

 while(!ifs2011Bmuid.eof() && iter < 28)   
   {   
     ifs2011Bmuid >> etalo >> etahi >> ptlo >> pthi >> eff >> errlo >> errhi;   
     muoneffs[iter][0] = etalo;  
     muoneffs[iter][1] = etahi;   
     muoneffs[iter][2] = ptlo;   
     muoneffs[iter][3] = pthi;   
     muoneffs[iter][4] = eff;   
     muoneffs[iter][5] = errlo;   
     muoneffs[iter][6] = errhi;   
     muoneffs[iter][7] = 2;   
     iter++;  
   }  
 ifs2011Bmuid.close();  

 if(0)
   {
     GetMuEHist(22,3,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs);
     return;
   }

 TCanvas *c1 = new TCanvas("c1","c1");
 if(mode == 2)
   c1->Divide(1,2);

 TH1F *htmp[13]; 
 // TH2F *htmp[13];

 htmp[0] = GetMuEHist(thevar,1,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs);
  htmp[1] = GetMuEHist(thevar,2,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs);
  htmp[2] = GetMuEHist(thevar,3,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs);
  htmp[3] = GetMuEHist(thevar,4,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs);
  htmp[4] = GetMuEHist(thevar,5,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs);
  htmp[5] = GetMuEHist(thevar,6,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs);
  htmp[6] = GetMuEHist(thevar,7,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs);
  if((cutset != 3) && (cutset < 777))
    htmp[7] = GetMuEHist(thevar,8,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs); 
  if(cutset > 999)
    htmp[7] = GetMuEHist(thevar,8,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs);  
  htmp[8] = GetMuEHist(thevar,9,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs);
  htmp[10] = GetMuEHist(thevar,11,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs); 
  htmp[11] = GetMuEHist(thevar,12,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs);  
  htmp[12] = GetMuEHist(thevar,13,cutset,save,muoneffs,muonhlt17effs,muonhlt8effs,eleeffs);   

  // JH - rescaled W+jets from anti-ID sample
  TString wjetstemplatefile;
  TString wjetstemplathist;
  
  if(cutset == 3)
    {
      if(thevar == 1)
	{
	  wjetstemplatehist = "hntrk";
	  wjetstemplatefile = "plotsnew/nextratracks_pt30_invertleptonIDrescaled.root";
	}
      if(thevar == 2)
	{
	  wjetstemplatehist = "hmll"; 
	  wjetstemplatefile = "plotsnew/memu_1to6tracks_pt30_invertleptonIDrescale.root";
	}
      if(thevar == 3)
	{
	  wjetstemplatehist = "hdphi";  
	  wjetstemplatefile = "plotsnew/dphiemu_1to6tracks_pt30_invertleptonIDrescale.root"; 
	}
      if(thevar == 4) 
        { 
          wjetstemplatehist = "hdpt";   
          wjetstemplatefile = "plotsnew/dptemu_1to6tracks_pt30_invertleptonIDrescale.root";  
        } 
      if(thevar == 6) 
        { 
          wjetstemplatehist = "hmet";   
          wjetstemplatefile = "plotsnew/pfmet_1to6tracks_pt30_invertleptonIDrescale.root";  
        } 
      if(thevar == 12)   
        {   
          wjetstemplatehist = "hmupt";     
          wjetstemplatefile = "plotsnew/ptmu_1to6tracks_pt30_invertleptonIDrescale.root";    
        }   
      if(thevar == 13)    
        {    
          wjetstemplatehist = "hept";      
          wjetstemplatefile = "plotsnew/etele_1to6tracks_pt30_invertleptonIDrescale.root";     
        }    
      if(thevar == 14)    
        {    
          wjetstemplatehist = "hmueta";      
          wjetstemplatefile = "plotsnew/etamu_1to6tracks_pt30_invertleptonIDrescale.root";     
        }    
      if(thevar == 15)     
        {     
          wjetstemplatehist = "heeta";       
          wjetstemplatefile = "plotsnew/etaele_1to6tracks_pt30_invertleptonIDrescale.root";      
        }     

      if(thevar == 16)  
        {  
          wjetstemplatehist = "hptextra";    
          wjetstemplatefile = "plotsnew/ptextratracks_1to6tracks_pt30_invertleptonIDrescale.root";   
        }  
      if(thevar == 17)   
        {   
          wjetstemplatehist = "hetaextra";     
          wjetstemplatefile = "plotsnew/etaextratracks_1to6tracks_pt30_invertleptonIDrescale.root";    
        }   
      if(thevar == 18)   
        {   
          wjetstemplatehist = "hsumptextra";     
          wjetstemplatefile = "plotsnew/sumptextratracks_1to6tracks_pt30_invertleptonIDrescale.root";    
        }   
      if(thevar == 24)
	{
	  wjetstemplatehist = "hdzextra";
	  wjetstemplatefile = "plotsnew/dzextra_1to6tracks_pt30_invertleptonIDrescale.root";
	}
      if(thevar == 25) 
        { 
          wjetstemplatehist = "hdxyextra"; 
          wjetstemplatefile = "plotsnew/dxyextra_1to6tracks_pt30_invertleptonIDrescale.root"; 
        } 
      if(thevar == 26)  
        {  
          wjetstemplatehist = "hchi2extra";  
          wjetstemplatefile = "plotsnew/chi2extra_1to6tracks_pt30_invertleptonIDrescale.root";  
        }  
      if(thevar == 27)  
        {  
          wjetstemplatehist = "hhitsextra";  
          wjetstemplatefile = "plotsnew/hitsextra_1to6tracks_pt30_invertleptonIDrescale.root";  
        }  
      if(thevar == 28)  
        {  
          wjetstemplatehist = "hpurityextra";  
          wjetstemplatefile = "plotsnew/purityextra_1to6tracks_pt30_invertleptonIDrescale.root";  
        }  
      if(thevar == 29)
	{
	  wjetstemplatehist = "hchi2pv";   
          wjetstemplatefile = "plotsnew/chi2pv_1to6tracks_pt30_invertleptonIDrescale.root";   

	}

      cout << "Opening " << wjetstemplatefile << endl;
      TFile *fr = TFile::Open(wjetstemplatefile);
      htmp[7] = (TH1F *)fr->Get(wjetstemplatehist); 
      htmp[7]->SetFillColor(11);
      htmp[7]->SetLineWidth(0);
    }

  if(cutset == 777)
    {
      if(thevar == 2)
	{
          wjetstemplatehist = "hmll";     
          wjetstemplatefile = "plotsnew/memu_0tracks_ptlessthan30_invertleptonIDrescale.root"; 
	}
      if(thevar == 3)
	{
	  wjetstemplatehist = "hdphi";    
          wjetstemplatefile = "plotsnew/dphiemu_0tracks_ptlessthan30_invertleptonIDrescale.root";
	}

      TFile *fr = TFile::Open(wjetstemplatefile); 
      htmp[7] = (TH1F *)fr->Get(wjetstemplatehist);  
      htmp[7]->SetFillColor(11); 
      htmp[7]->SetLineWidth(0); 
    }

  if(cutset == 888 || cutset == 999)
    {
      if(thevar == 2) 
        { 
          wjetstemplatehist = "hmll";  
          wjetstemplatefile = "plotsnew/memu_0tracks_pt30_invertleptonIDrescale.root"; 
        } 
      if(thevar == 3) 
        { 
          wjetstemplatehist = "hdphi";   
          wjetstemplatefile = "plotsnew/dphiemu_0tracks_pt30_invertleptonIDrescale.root";  
        } 
      if(thevar == 4)  
        {  
          wjetstemplatehist = "hdpt";    
          wjetstemplatefile = "plotsnew/dptemu_0tracks_pt30_invertleptonIDrescale.root";   
        }  
      if(thevar == 5)
	{
	  wjetstemplatehist = "hpt";
	  wjetstemplatefile = "plotsnew/ptemu_0tracks_unblindAQGC_invertleptonIDrescale.root";
	}
      if(thevar == 6)  
        {  
          wjetstemplatehist = "hmet";    
          wjetstemplatefile = "plotsnew/pfmet_0tracks_pt30_invertleptonIDrescale.root";   
        }  
      if(thevar == 7)
	{
	  wjetstemplatehist = "hnvrt";
	  wjetstemplatefile = "plotsnew/npv_0tracks_pt30_invertleptonIDrescale.root";
	}
      if(thevar == 12)  
        {  
          wjetstemplatehist = "hmupt";   
          wjetstemplatefile = "plotsnew/ptmu_0tracks_pt30_invertleptonIDrescale.root";  
        }  
      if(thevar == 13)  
        {  
          wjetstemplatehist = "hept";    
          wjetstemplatefile = "plotsnew/etele_0tracks_pt30_invertleptonIDrescale.root";   
        }  
      if(thevar == 14)   
        {   
          wjetstemplatehist = "hmueta";     
          wjetstemplatefile = "plotsnew/etamu_0tracks_pt30_invertleptonIDrescale.root";    
        }   
      if(thevar == 15)   
        {   
          wjetstemplatehist = "heeta";     
          wjetstemplatefile = "plotsnew/etaele_0tracks_pt30_invertleptonIDrescale.root";    
        }   
      TFile *fr = TFile::Open(wjetstemplatefile); 
      htmp[7] = (TH1F *)fr->Get(wjetstemplatehist);  
      htmp[7]->SetFillColor(11); 
      htmp[7]->SetLineWidth(0); 
    }

  /*
  htmp[7]->SetBinContent(1,0.2373519);
  htmp[7]->SetBinContent(2,0.2769105);
  htmp[7]->SetBinContent(3,1.068083);
  htmp[7]->SetBinContent(4,1.265877);
  htmp[7]->SetBinContent(5,1.424111);
  htmp[7]->SetBinContent(6,1.384553);
  htmp[7]->SetBinContent(7,1.621904);
  htmp[7]->SetBinContent(8,3.322926);
  htmp[7]->SetBinContent(9,3.20425);
  htmp[7]->SetBinContent(10,5.142624);
  htmp[7]->SetBinContent(11,5.775562);
  htmp[7]->SetBinContent(12,7.990847);
  htmp[7]->SetBinContent(13,9.019371);
  htmp[7]->SetBinContent(14,11.94671);
  htmp[7]->SetBinContent(15,14.91361);
  htmp[7]->SetBinContent(16,16.21905);
  htmp[7]->SetBinError(1,0.09689851);
  htmp[7]->SetBinError(2,0.1046623);
  htmp[7]->SetBinError(3,0.2055527);
  htmp[7]->SetBinError(4,0.2237775);
  htmp[7]->SetBinError(5,0.2373519);
  htmp[7]->SetBinError(6,0.2340321);
  htmp[7]->SetBinError(7,0.2532989);
  htmp[7]->SetBinError(8,0.362561);
  htmp[7]->SetBinError(9,0.3560278);
  htmp[7]->SetBinError(10,0.4510378);
  htmp[7]->SetBinError(11,0.4779891);
  htmp[7]->SetBinError(12,0.562234);
  htmp[7]->SetBinError(13,0.5973225);
  htmp[7]->SetBinError(14,0.6874561);
  htmp[7]->SetBinError(15,0.7680901);
  htmp[7]->SetBinError(16,0.8010017);
  */
  // end JH

  Float_t inclww = htmp[2]->GetSumOfWeights();
  Float_t inclwwstats = htmp[2]->GetEntries();
  Float_t inclwwerr = sqrt(inclwwstats)*inclww/inclwwstats;
  Float_t gamgamtautau = htmp[4]->GetSumOfWeights(); 
  Float_t gamgamtautaustats = htmp[4]->GetEntries(); 
  Float_t gamgamtautauerr = sqrt(gamgamtautaustats)*gamgamtautau/gamgamtautaustats;
  Float_t inelgamgamtautau = htmp[8]->GetSumOfWeights();
  Float_t inelgamgamtautaustats = htmp[8]->GetEntries();
  Float_t inelgamgamtautauerr = sqrt(inelgamgamtautaustats)*inelgamgamtautau/inelgamgamtautaustats;
  Float_t dytautau = htmp[1]->GetSumOfWeights();
  Float_t dytautaustats = htmp[1]->GetEntries();
  Float_t dytautauerr = sqrt(dytautaustats)*dytautau/dytautaustats;
  Float_t diffww = htmp[3]->GetSumOfWeights(); 
  Float_t diffwwstats = htmp[3]->GetEntries(); 
  Float_t diffwwerr = sqrt(diffwwstats)*diffww/diffwwstats; 
  Float_t wjets = htmp[7]->GetSumOfWeights();  
  Float_t wjetsstats = htmp[7]->GetEntries();  
  Float_t wjetserr = sqrt(wjetsstats)*wjets/wjetsstats;  
  Float_t aqgc1stats = htmp[11]->GetEntries();
  Float_t aqgc2stats = htmp[12]->GetEntries();

  //  htmp[2]->Scale(44.0/38.9); 
  //  htmp[3]->Scale(44.0/38.9); 
  //  htmp[2]->Scale(44.0/43.48);
  //  htmp[3]->Scale(44.0/43.48);
  //  htmp[2]->Scale(44.0/117.093);
  //  htmp[3]->Scale(44.0/117.093);

  htmp[1]->Add(htmp[2]);
  htmp[1]->Add(htmp[3]); 
  htmp[1]->Add(htmp[4]); 
  htmp[1]->Add(htmp[5]);
  htmp[1]->Add(htmp[6]); 
  htmp[1]->Add(htmp[7]);  
  htmp[1]->Add(htmp[8]);

  htmp[2]->Add(htmp[3]);  
  htmp[2]->Add(htmp[4]);  
  htmp[2]->Add(htmp[5]);   
  htmp[2]->Add(htmp[6]);    
  htmp[2]->Add(htmp[7]);     
  htmp[2]->Add(htmp[8]);

  htmp[3]->Add(htmp[4]);   
  htmp[3]->Add(htmp[5]);    
  htmp[3]->Add(htmp[6]);     
  htmp[3]->Add(htmp[7]);      
  htmp[3]->Add(htmp[8]);

  htmp[6]->Add(htmp[4]); 
  htmp[6]->Add(htmp[5]);  
  htmp[6]->Add(htmp[7]);  
  htmp[6]->Add(htmp[8]);

  htmp[7]->Add(htmp[4]);   
  htmp[7]->Add(htmp[5]);   
  htmp[7]->Add(htmp[8]);

  htmp[8]->Add(htmp[5]);    
  htmp[8]->Add(htmp[4]);

  htmp[4]->Add(htmp[5]);
  
  TH1F *herror=(TH1F*)htmp[1]->Clone();
  herror->SetName("herror");

  TH1F *hnoerror=(TH1F*)htmp[1]->Clone(); 
  hnoerror->SetName("hnoerror"); 
  for(Int_t i = 1; i < hnoerror->GetNbinsX(); i++)
    hnoerror->SetBinError(i,0);

  TH1F *hsub=(TH1F*)htmp[0]->Clone();
  hsub->SetName("hsub");

  TH1F *hmcsub=(TH1F*)htmp[1]->Clone(); 
  hmcsub->SetName("hmcsub"); 


  if(htmp[1]->GetMaximum() > htmp[0]->GetMaximum())
    htmp[1]->SetMaximum(1.75 * htmp[1]->GetMaximum());
  else
    {
      if(htmp[0]->GetMaximum() > 1)
	htmp[1]->SetMaximum(2.5 * htmp[0]->GetMaximum());
      else
	htmp[1]->SetMaximum(4.0);
      if(cutset>=999)
	{
	  if(htmp[12]->GetMaximum() > htmp[1]->GetMaximum())
	    htmp[1]->SetMaximum(htmp[12]->GetMaximum() * 1.5);
	}
    }

  c1->cd(1);
  if(mode == 2) 
    {
      htmp[1]->SetMinimum(0.01);
      htmp[1]->SetMaximum(10.0* htmp[1]->GetMaximum()); 
    }

  htmp[1]->SetStats(0);
  htmp[1]->SetTitle(0);

  // JH
  htmp[10]->Add(htmp[1]);

  htmp[1]->Draw("hist");
  htmp[2]->Draw("histsame");
  htmp[3]->Draw("histsame"); 
  htmp[6]->Draw("histsame");
  htmp[7]->Draw("histsame");
  htmp[8]->Draw("histsame");
  htmp[4]->Draw("histsame"); 
  htmp[10]->Draw("histsame");  
  if(cutset >= 999)
   {
      htmp[11]->Add(htmp[1]);
      htmp[12]->Add(htmp[1]);

      htmp[11]->Draw("histsame");
      htmp[12]->Draw("histsame"); 
    }
  herror->SetMarkerStyle(0);
  herror->SetFillColor(1);
  herror->SetFillStyle(3004);
  herror->Draw("e2same");
  //  htmp[0]->Draw("e0same");

  // Poisson error bars
  const double alpha = 1 - 0.6827;
  TGraphAsymmErrors * g = new TGraphAsymmErrors(htmp[0]);

  for (int i = 0; i < g->GetN(); ++i) {
    int N = g->GetY()[i];
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  (N==0) ?  ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ): (ROOT::Math::gamma_quantile_c(alpha/2,N+1,1));
    if(N==0)
      U=0;
    g->SetPointEXlow(i, 0);
    g->SetPointEXhigh(i, 0);
    g->SetPointEYlow(i, N-L);
    g->SetPointEYhigh(i, U-N);
  }
  g->Draw("P");

  TLegend *l1 = new TLegend(0.55,0.6,0.65,0.9);
  l1->AddEntry(htmp[0],"Data","lf");
  l1->AddEntry(htmp[1],"Drell-Yan #tau^{+}#tau^{-}","lf"); 
  l1->AddEntry(htmp[2],"Inclusive W^{+}W^{-}","lf"); 
  l1->AddEntry(htmp[3],"Diffractive W^{+}W^{-}","lf");  
  l1->AddEntry(htmp[4],"Elastic #gamma#gamma #rightarrow #tau^{+}#tau^{-}","lf");  
  l1->AddEntry(htmp[8],"Inelastic #gamma#gamma #rightarrow #tau^{+}#tau^{-}","lf");
  l1->AddEntry(htmp[6],"t#bar{t}","lf"); 
  l1->AddEntry(htmp[7],"W+jets","lf");
  l1->AddEntry(htmp[10],"#gamma#gamma #rightarrow W^{+}W^{-} (SM)","lf");   
  //  l1->AddEntry(htmp[1],"POWHEG-PYTHIA Z2 DY #tau^{+}#tau^{-}","lf");
  //  l1->AddEntry(htmp[2],"Madgraph W^{+}W^{-}","lf");
  //  l1->AddEntry(htmp[2],"PYTHIA  W^{+}W^{-}","lf");
  //  l1->AddEntry(htmp[2],"MC@NLO W^{+}W^{-}","lf");
  //  l1->AddEntry(htmp[3],"POMPYT diffractive W^{+}W^{-}","lf"); 
  //  l1->AddEntry(htmp[4],"LPAIR #gamma#gamma #rightarrow #tau^{+}#tau^{-}","lf"); 
  //  l1->AddEntry(htmp[6],"Madgraph ttbar","lf");
  //  l1->AddEntry(htmp[10],"#gamma#gamma #rightarrow W^{+}W^{-} (SM)","lf");  

  if(cutset >= 999)
    {
      l1->AddEntry(htmp[11],"#gamma#gamma #rightarrow W^{+}W^{-} (a0W=2E-4, aCW=0, #Lambda=500GeV)","lf");   
      l1->AddEntry(htmp[12],"#gamma#gamma #rightarrow W^{+}W^{-} (a0W=-2E-4, aCW=-8E-4, #Lambda=500GeV)","lf");    
    }
  
  l1->SetFillColor(0); // l1->SetTextSize(0.03);
  l1->SetTextSize(18); 
  l1->SetTextFont(43); 
  l1->Draw("same");

  stringstream ss;
  ss.str(""); ss << "CMS Preliminary 2011, #sqrt{s}=7 TeV, L=5.05 fb^{-1}";
  //  TPaveText *plotlabel = new TPaveText(0.5,0.96,0.92,0.99,"NDC");
  TPaveText *plotlabel = new TPaveText(0.4,0.96,0.92,1.0,"NDC");
  plotlabel->SetTextColor(kBlack);
  plotlabel->SetFillColor(kWhite);
  plotlabel->SetBorderSize(0);
  plotlabel->SetTextAlign(32);
  plotlabel->SetTextSize(20);
  plotlabel->SetTextFont(43);
  plotlabel->AddText(ss.str().c_str());
  plotlabel->Draw("same");

  if(mode == 2)
    {
      c1_1->SetLogy();
      c1->cd(2); 
      //      hsub->Sumw2();
      hsub->Divide(hnoerror);
      hmcsub->Sumw2();
      hmcsub->Divide(herror);

      //      hsub->Add(herror,-1);
      //      hsub->Divide(herror);
      Float_t themin = (-1.0 * hsub->GetMaximum()) - 2;
      Float_t themax = (-1.0 * themin) + 1;
      //      Float_t themin = (-1.0 * hsub->GetMaximum());
      //      hsub->SetMinimum(themin);
      if(themax > 4) themax = 4.0;
      hsub->SetMinimum(0);
      hmcsub->SetMaximum(themax);
      //      hmcsub->SetMinimum(themin);
      hmcsub->SetMinimum(0);
      //      hsub->SetMaximum(5.0); hsub->SetMinimum(-3.0);
      //      hsub->SetMaximum(5.0); hsub->SetMinimum(-5.0);
      hmcsub->SetTitle(0); hmcsub->SetStats(0); hmcsub->SetMarkerColor(6);
      hmcsub->Draw("e2");
      hsub->SetTitle(0); hsub->SetStats(0);
      hsub->Draw("esame");
      Int_t linemin = hsub->GetXaxis().GetXmin();
      Int_t linemax = hsub->GetXaxis().GetXmax();
      TLine *l2 = new TLine(linemin,1,linemax,1);
      //      TLine *l2 = new TLine(linemin,0,linemax,0);
      TLine *l3 = new TLine(linemin,1,linemax,1);
      TLine *l4 = new TLine(linemin,-1,linemax,-1);
      l2->SetLineColor(2);
      l2->SetLineWidth(3);
      l2->SetLineStyle(2);
      l3->SetLineColor(2);
      l3->SetLineWidth(3);
      l3->SetLineStyle(2);
      l4->SetLineColor(2);
      l4->SetLineWidth(3);
      l4->SetLineStyle(2);
      l2->Draw("same");
      //      l3->Draw("same");
      //      l4->Draw("same");
    }

  // MC stat errors
  Float_t sumofsq = 0.0;
  for(Int_t k = 0; k < htmp[1]->GetNbinsX(); k++)
    {
      sumofsq = sumofsq + ((htmp[1]->GetBinError(k+1))*(htmp[1]->GetBinError(k+1)));
    }

  cout << "Data = " << htmp[0]->GetSumOfWeights() << endl;
  cout << "MC SM Signal + Background = " << htmp[10]->GetSumOfWeights() << " ("
       << htmp[10]->GetEntries() << ")" << endl;
  cout << "MC Sum of bkg = " << htmp[1]->GetSumOfWeights() << " +- " 
       << sqrt(sumofsq) << endl;
  cout << "\tWW = " << inclww << " +- " << inclwwerr << " (" << inclwwstats << ")" << endl;
  cout << "\tDiff WW = " << diffww << " +- " << diffwwerr << " (" << diffwwstats << ")" << endl; 
  cout << "\tW+Jets = " << wjets << " +- " << wjetserr << " (" << wjetsstats << ")" << endl;
  cout << "\tgamgam->tautau = " << gamgamtautau << " +- " << gamgamtautauerr << " (" << gamgamtautaustats << ")" << endl; 
  cout << "\tgamgam->tautau = " << inelgamgamtautau << " +- " << inelgamgamtautauerr << " (" << inelgamgamtautaustats << ")" << endl;
  cout << "\tDY->tautau = " << dytautau << " +- " << dytautauerr << " (" << dytautaustats << ")" << endl; 
  cout << "MC AQGC Signal (Point 1) + Background = " << htmp[11]->GetSumOfWeights() << " (" << aqgc1stats << ")" <<endl;
  cout << "MC AQGC Signal (Point 2) + Background = " << htmp[12]->GetSumOfWeights() << " (" << aqgc2stats << ")" <<endl;

  if(save == 1)
    {
      TString savedir = "plotsnew/";
      TString saveext = ".pdf";
      TString fullsavename = savedir + savename + saveext;
      c1->SaveAs(fullsavename);
      TString saveext = ".png"; 
      fullsavename = savedir + savename + saveext;
      c1->SaveAs(fullsavename); 
    }

  //	gPad->SetLogy();
}

// Return histogrammed quantities for mu+mu- samples
TH1F *GetMuEHist(Int_t plotvar = 1,Int_t physsample = 1, Int_t thecuts, bool save = false, Double_t efftable[28][8], Double_t effhltmu17table[17][8], Double_t effhltmu8table[17][8], Double_t effeletable[30][8])
//TH2F *GetMuEHist(Int_t plotvar = 1,Int_t physsample = 1, Int_t thecuts, bool save = false, Double_t efftable[28][8], Double_t effhltmu17table[17][8], Double_t effhltmu8table[17][8], Double_t effeletable[30][8])
{
  Double_t lumi = 5.05; // Run 2011 A+B
  // lumi = 2.31 // Run 2011 A
  // lumi = 2.74; // Run 2011 B
  Double_t xsec = 1.0;
  TString st = "";
  Int_t linecolor = 1;
  Int_t linewidth = 0;
  Int_t linestyle = 1;
  Int_t fillcolor = kWhite; //15;
  Int_t fillstyle = 1001;
  Int_t nent = 0;

  Double_t pTemumin = 0.0;
  Double_t pTemumax = 99999.0;
  Double_t massemumin = 20.0;
  Int_t extratracksmax = 15; 
  Int_t extratracksmin = 0;
  

  switch(thecuts) {
  case 1:
    // Extra tracks sideband
    extratracksmin = 1;
    extratracksmax = 6;
    cout << "\tm(emu) > 20 GeV" << endl;
    cout << "\tpT(emu) > 0 GeV" << endl;    
    cout << "\t1 < N(extra tracks) < 6" << endl;
    break;
  case 2:
    // High pT, no tracks cut
    pTemumin = 30.0;
    cout << "\tm(emu) > 20 GeV" << endl; 
    cout << "\tpT(emu) > 30 GeV" << endl;  
    break;
  case 3:
    // Extra tracks sideband, pT cut
    pTemumin = 30.0; 
    extratracksmin = 1;  
    extratracksmax = 6;  
    cout << "\tm(emu) > 20 GeV" << endl;   
    cout << "\tpT(emu) > 30 GeV" << endl;    
    cout << "\t1 < N(extra tracks) < 6" << endl;  
    break;
  case 4:
    // Extra tracks sideband, pT anti-cut
    pTemumax = 30.0;  
    extratracksmin = 1;   
    extratracksmax = 6;   
    cout << "\tm(emu) > 20 GeV" << endl;    
    cout << "\tpT(emu) < 30 GeV" << endl;     
    cout << "\t1 < N(extra tracks) < 6" << endl;   
    break; 
  case 5:
    // Extra tracks sideband, no mass cut  
    extratracksmin = 1; 
    extratracksmax = 6; 
    massemumin = 0.0;
    cout << "\tm(emu) > 0 GeV" << endl;  
    cout << "\tpT(emu) > 0 GeV" << endl;   
    cout << "\t1 < N(extra tracks) < 6" << endl; 
    break;
  case 6:
    // No extra tracks cut, no pT cut
    extratracksmin = 0;  
    extratracksmax = 999;  
    massemumin = 0.0; 

    cout << "\tm(emu) > 20 GeV" << endl;   
    cout << "\tpT(emu) > 0 GeV" << endl;    
    break; 
  case 777:
    // Unblind - low pT tautau region
    extratracksmin = 0;
    extratracksmax = 0;
    pTemumax = 30;
    cout << "\tm(emu) > 20 GeV" << endl;
    cout << "\tpT(emu) < 30 GeV" << endl;
    cout << "\tN(extra tracks) = 0" << endl;
    break;
  case 888:
    // Unblind - SM region
    extratracksmin = 0; 
    extratracksmax = 0; 
    pTemumin = 30.0;
    cout << "\tm(emu) > 20 GeV" << endl;   
    cout << "\tpT(emu) > 30 GeV" << endl;    
    cout << "\tN(extra tracks) = 0" << endl;  
    break;
  case 999:
    // Unblind - SM region no pT cut
    extratracksmin = 0;  
    extratracksmax = 0;  
    pTemumin = 0.0; 
    cout << "\tm(emu) > 20 GeV" << endl;    
    cout << "\tpT(emu) > 0 GeV" << endl;     
    cout << "\tN(extra tracks) = 0" << endl;   
    break; 
  case 1111:
    // Unblind - AQGC
    extratracksmin = 0;  
    extratracksmax = 0;  
    pTemumin = 100.0; 
    massemumin = 20.0;
    cout << "\tm(emu) > 20 GeV" << endl;    
    cout << "\tpT(emu) > 100 GeV" << endl;     
    cout << "\tN(extra tracks) = 0" << endl;   
    break;
  case 2222: 
    // Unblind - pT30, no extra tracks cut
    pTemumin = 30.0;  
    massemumin = 20.0; 
    cout << "\tm(emu) > 20 GeV" << endl;     
    cout << "\tpT(emu) > 30 GeV" << endl;      
    break; 
  case 3333:
    // Unblind - pT100, no extra tracks cut
    pTemumin = 100.0;
    massemumin = 20.0;
    cout << "\tm(emu) > 20 GeV" << endl;
    cout << "\tpT(emu) > 100 GeV" << endl;
    break;
  default:
    break;
  }

  switch(physsample) {
  case 1:
    xsec = 1.0;
    st = "ExclusiveMuE2011ABReRecoNov08_Round4/ExclusiveMuE2011ABReRecoNov08_Round4_merge.root";
    //    st = "ExclusiveMuE2011ABReRecoNov08_Round4/SSExclusiveMuE2011ABReRecoNov08_Round4_merge.root";
    linecolor = 1;
    linewidth = 3;
    fillcolor = 1;
    //    fillstyle = 3001;
    fillstyle = 3004;
    break;
  case 2:
    //    xsec = 1300.0 *1000 * lumi / 19937479.0;
    xsec = 1666.0 * 1000 * lumi / 19937479.0; // NNLO
    st = "ExclusiveMuE_DYToTauTau_CT10_M20_Z2_START44_Fall11_Round4/ExclusiveMuE_DYToTauTau_CT10_M20_Z2_START44_Fall11_Round4_merge.root";
    fillcolor = 6;
    break;
  case 3:
    // MadGraph WW - NLO cross-section 
    xsec = 4.79 * 1000 * lumi / 1197558.0;
    st = "ExclusiveMuE_WWJetsTo2L2Nu_Z2_MadgraphTauola_Fall11_Round4/ExclusiveMuE_WWJetsTo2L2Nu_Z2_MadgraphTauola_Fall11_Round4_merge.root";
    // MC@NLO - NLO cross-section
    //    xsec = 4.79 * 1000 * lumi / 539746.0;
    //    st = "ExclusiveMuE_WW_MCATNLO_Fall11_Round4/ExclusiveMuE_WW_MCATNLO_Fall11_Round4_merge.root";
    // Pythia6 - NLO cross-section  
    //    xsec = 4.79 * 1000 * lumi / 210667.0;
    //    st = "ExclusiveMuE_WWTo2L2Nu_Z2_START44_Pythia6Tauola_Fall11_Round4/ExclusiveMuE_WWTo2L2Nu_Z2_START44_Pythia6Tauola_Fall11_Round4_merge.root";
    fillcolor = 7;
    break;
  case 4:
    xsec = 77.0 * lumi / 20000.0; 
    st = "ExclusiveMuE_WW_Pompyt_Fall11_Round4/ExclusiveMuE_WW_Pompyt_Fall11_Round4_merge.root";
    fillcolor = 96;
    break;
  case 5:
    //    xsec = 2.0 * 0.087 * 1000 * lumi / 10000.0; pT>25
    //    xsec = 2.0 * 1.02 * 1000 * lumi / 100000.0; // pT>10
    //    xsec = 2.0 * 0.18 * 1000 * lumi * 0.42542 / 3000.0; // pT>19, GEN-level filter
    //    xsec = 2.53 * 0.18 * 1000 * lumi * (6694.0/3000.0) / 100000.0; // pT>19, 3k events, emu only - default
    xsec = 0.022 * 1000 * lumi / 100000.0; // pT 40GeV 
    //    st = "ExclusiveMuE_GammaGammaTauTau_LPAIR_pT25_Fall11_Round4/GammaGammaTauTau_LPair_Fall11_S6.root";
    //    st = "ExclusiveMuE_GammaGammaTauTau_LPAIR_pT25_Fall11_Round4/GammaGammaTauTau_LPAIR_pT10_Fall11_Round4_merge.root";
    //    st = "ExclusiveMuE_GammaGammaTauTau_LPAIR_pT25_Fall11_Round4/GammaGammaTauTau_LPAIR_pT19_1prong_Fall11_Round4_merge.root";
    // JH - default
    //    st = "ExclusiveMuE_GammaGammaTauTau_LPAIR_pT25_Fall11_Round4/GammaGammaTauTau_LPAIR_pT19_emuonly_Fall11_Round4_merge.root";
    // JH - 40GeV cut
    st = "ExclusiveMuE_GammaGammaTauTau_LPAIR_pT25_Fall11_Round4/GammaGammaTauTau_LPAIR_ElEl_pT40_Fall11_S6_Round4_.root";
    fillcolor = 5;
    break;
  case 6: 
    //    xsec = 2.53 * 0.18 * 1000 * lumi * (35848.0/3000.0) / 100000.0; // pT>19, 3k events, lepton+hadron only, 2.53 scale factor  
    xsec = 0.18 * 1000 * lumi * (35848.0/3000.0) / 100000.0; // pT>19, 3k events, lepton+hadron only, no scale factor - default
    // JH - default
    st = "ExclusiveMuE_GammaGammaTauTau_LPAIR_pT25_Fall11_Round4/GammaGammaTauTau_LPAIR_pT19_1hadrononly_Fall11_Round4_merge.root"; 
    fillcolor = 5;  
    break;  
  case 7:
    xsec = 157.5 * 1000 * lumi / 3701947.0; 
    st = "ExclusiveMuE_TTbarJets_Z2_Madgraph_Fall11_Round4/ExclusiveMuE_TTbarJets_Z2_Madgraph_Fall11_Round4_merge.root";
    fillcolor = 4;
    break;
  case 8:
    xsec = 31314 * 1000 * lumi / 81345381.0;
    st = "ExclusiveMuE_WJetsToLNu_Z2_Madgraph_Tauola_Fall11_Round4/ExclusiveMuE_WJetsToLNu_Z2_Madgraph_Tauola_Fall11_Round4_merge.root";
    fillcolor = 11;
    break;
  case 9:
    //    xsec = 0.8099 * 1000 * lumi / 100000.0; // default
    xsec = 0.59 * 2 * 0.031 * 1000 * lumi / 100000.0; // 40GeV pT cut, 0.59 scale factor derived from mumu 
    // JH - default
    //    st = "ExclusiveMuE_GammaGammaTauTau_LPAIR_pT25_Fall11_Round4/GammaGammaTauTau_LPAIR_InelEl_SM_Fall11_S6_Round4_.root";
    // JH - 40GeV cut 
    st = "ExclusiveMuE_GammaGammaTauTau_LPAIR_pT25_Fall11_Round4/GammaGammaTauTau_LPAIR_InelEl_pT40_Fall11_S6_Round4_.root"; 
    fillcolor = 76;
    break;
  case 11:
    //    xsec = 2.53 * 40.0 * lumi * (4*0.107*0.107) / 10000.0; 
    xsec = 3.95 * 40.0 * lumi * (4*0.107*0.107) / 10000.0;
    st = "ExclusiveMuE_GammaGammaWW_CalcHep_SM_START44_Fall11_Round4/GammaGammaWW_CalcHep_SM_Fall11_S6_Round4_merge.root";
    //    fillcolor = 3;
    fillcolor = 0;
    linecolor = 3; 
    linewidth = 3; 
    break;
  case 12:
    //    xsec = 2.53 * 138.0 * lumi * (4*0.107*0.107) / 5000.0;  
    xsec = 3.95 * 138.0 * lumi * (4*0.107*0.107) / 5000.0;
    st = "ExclusiveMuE_GammaGammaWW_CalcHep_SM_START44_Fall11_Round4/GammaGammaWW_CalcHep_Anomalous1_Gustavo_Fall11_S6_Round4_merge.root"; 
    fillcolor = 0; 
    linecolor = 3;  
    linewidth = 3;  
    linestyle = 2;
    break;
  case 13: 
    //    xsec = 2.0 * 369.6 * lumi * (4*0.107*0.107) / 5000.0;   
    xsec = 3.95 * 369.6 * lumi * (4*0.107*0.107) / 5000.0;
    st = "ExclusiveMuE_GammaGammaWW_CalcHep_SM_START44_Fall11_Round4/GammaGammaWW_CalcHep_Anomalous2_Gustavo_Fall11_S6_Round4_merge.root";  
    fillcolor = 0;  
    linecolor = 3;   
    linewidth = 3;   
    linestyle = 3; 
    break; 
  default:
    break;
  }

  TFile *f1 = new TFile(st);

  TTree *tr1 = f1->Get("ntp1");

  TH1F *hmll;
  TH1F *hdphi;
  TH1F *hdpt;
  TH1F *hpt;
  TH1F *hmet;
  TH1F *hntrk;
  TH1F *hnvrt; 
  TH1F *hptextra;
  TH1F *hetaextra;
  TH1F *hmuemetdphi;
  TH1F *hvtxz;
  TH1F *hvtxT;
  TH1F *hemudist;
  TH1F *hsumptextra;
  TH1F *hsumpzextra;
  TH1F *hdzextra, *hdxyextra;
  TH1F *hchi2extra, *hhitsextra, *hpurityextra, *hchi2pv; 
  TH1F *hmupt, *hmueta, *hept, *heeta;
  TH1F *hleadpt, *htrailpt;
  TH2F *hptntrackcorr;

  if(thecuts < 777) 
    hmll = new TH1F("hmll","hmll",50,0,300);
    //    hmll = new TH1F("hmll","hmll",25,0,300);
  if(thecuts == 777)
    hmll = new TH1F("hmll","hmll",10,0,300); 
  if(thecuts == 3)
    //    hmll = new TH1F("hmll","hmll",25,0,300);
    hmll = new TH1F("hmll","hmll",15,0,300);
  if(thecuts > 777) 
    hmll = new TH1F("hmll","hmll",10,0,1000);
  if(thecuts < 777)
    hdphi = new TH1F("hdphi","hdphi",50,0.0,1.0);
  if(thecuts == 3)
    //    hdphi = new TH1F("hdphi","hdphi",20,0.0,1.0);
    hdphi = new TH1F("hdphi","hdphi",10,0.0,1.0);
  if(thecuts == 777)
    hdphi = new TH1F("hdphi","hdphi",20,0.0,1.0);
  if(thecuts > 777)
    hdphi = new TH1F("hdphi","hdphi",5,0.0,1.0);
  if(thecuts < 777) 
    hdpt = new TH1F("hdpt","hdpt",50,0,100);
    //    hdpt = new TH1F("hdpt","hdpt",10,0,100);
  if(thecuts == 777)
    hdpt = new TH1F("hdpt","hdpt",20,0,100); 
  if(thecuts > 777) 
    hdpt = new TH1F("hdpt","hdpt",10,0,100); 
  if(thecuts < 777)
    hpt = new TH1F("hpt","hpt",50,0,300); 
  if(thecuts >= 777)
    hpt = new TH1F("hpt","hpt",15,0,300);  
  //  hpt = new TH1F("hpt","hpt",12,0,300);
  if(thecuts < 777)
    hmet = new TH1F("hmet","hmet",50,0,300);
  if(thecuts == 777)
    hmet = new TH1F("hmet","hmet",5,0,300);  
  if(thecuts > 777)
    hmet = new TH1F("hmet","hmet",10,0,300); 
  //  hntrk = new TH1F("hntrk","hntrk",15,0,15);
  hntrk = new TH1F("hntrk","hntrk",15,-0.5,14.5);
  hnvrt = new TH1F("hnvrt","hnvrt",30,0,30); 
  hptextra = new TH1F("hptextra","hptextra",50,0,20);
  hetaextra = new TH1F("hetaextra","hetaextra",50,-3,3); 
  hvtxz = new TH1F("hvtxz","hvtxz",30,-30,30);
  hvtxT = new TH1F("hvtxT","hvtxT",500,0,0.6); 
  hmuemetdphi = new TH1F("hmuemetdphi","hmuemetdphi",20,0,1);
  hemudist = new TH1F("hemudist","hemudist",100,-0.05,0.05);
  hsumptextra = new TH1F("hsumptextra","hsumptextra",25,0,100);
  hsumpzextra = new TH1F("hsumpzextra","hsumpzextra",25,0,100); 
  hdzextra = new TH1F("hdzextra","hdzextra",100,0,0.5);
  hdxyextra = new TH1F("hdxyextra","hdxyextra",100,0,0.5);
  hchi2extra = new TH1F("hchi2extra","hchi2extra",50,0,5); 
  hhitsextra = new TH1F("hhitsextra","hhitsextra",35,0,35); 
  hpurityextra = new TH1F("hpurityextra","hpurityextra",2,0,2); 
  hchi2pv = new TH1F("hchi2pv","hchi2pv",50,0,5);

  if(thecuts == 6)
    {
      hmupt = new TH1F("hmupt","hmupt",50,0,250); 
      hmueta = new TH1F("hmueta","hmueta",30,-3.0,3.0);  
      hept = new TH1F("hept","hept",50,0,250);  
      heeta = new TH1F("heeta","heeta",30,-3.0,3.0);   
      hleadpt = new TH1F("hleadpt","hleadpt",50,0,250); 
      htrailpt = new TH1F("htrailpt","htrailpt",50,0,250); 
    }
  if((thecuts < 777) && (thecuts != 6))
    {
      hmupt = new TH1F("hmupt","hmupt",25,0,250);
      hmueta = new TH1F("hmueta","hmueta",15,-3.0,3.0); 
      hept = new TH1F("hept","hept",25,0,250); 
      heeta = new TH1F("heeta","heeta",15,-3.0,3.0);  
      hleadpt = new TH1F("hleadpt","hleadpt",25,0,250);
      htrailpt = new TH1F("htrailpt","htrailpt",25,0,250);
    }
  if(thecuts >= 777)
    {
      hmupt = new TH1F("hmupt","hmupt",5,0,250); 
      hmueta = new TH1F("hmueta","hmueta",6,-3.0,3.0);  
      hept = new TH1F("hept","hept",5,0,250);  
      heeta = new TH1F("heeta","heeta",5,-3.0,3.0);   
      hleadpt = new TH1F("hleadpt","hleadpt",5,0,250); 
      htrailpt = new TH1F("htrailpt","htrailpt",5,0,250); 
    }
  hptntrackcorr = new TH2F("hptntrackcorr","hptntrackcorr",15,0,15,20,0,300);
  //  hleadpt = new TH1F("hleadpt","hleadpt",50,0,200);
  //  htrailpt = new TH1F("htrailpt","htrailpt",50,0,200);

   // Declaration of leaf types
   Int_t           nJetCand;
   Double_t        JetCand_px[30];   //[nJetCand]
   Double_t        JetCand_py[30];   //[nJetCand]
   Double_t        JetCand_pz[30];   //[nJetCand]
   Double_t        JetCand_e[30];   //[nJetCand]
   Double_t        JetCand_eta[30];   //[nJetCand]
   Double_t        JetCand_phi[30];   //[nJetCand]
   Double_t        HighestJet_e;
   Double_t        HighestJet_eta;
   Double_t        HighestJet_phi;
   Double_t        SumJet_e;
   Int_t           nMuonCand;
   Double_t        MuonCand_px[10];   //[nMuonCand]
   Double_t        MuonCand_py[10];   //[nMuonCand]
   Double_t        MuonCand_pz[10];   //[nMuonCand]
   Double_t        MuonCand_vtxx[10];   //[nMuonCand]
   Double_t        MuonCand_vtxy[10];   //[nMuonCand]
   Double_t        MuonCand_vtxz[10];   //[nMuonCand]
   Double_t        MuonCand_p[10];   //[nMuonCand]
   Double_t        MuonCand_pt[10];   //[nMuonCand]
   Double_t        MuonCand_eta[10];   //[nMuonCand]
   Double_t        MuonCand_phi[10];   //[nMuonCand]
   Int_t           MuonCand_charge[10];   //[nMuonCand]
   Int_t           MuonCand_tmlsloosemuonid[10];   //[nMuonCand]
   Int_t           MuonCand_tmlsOptLowPtloosemuonid[10];   //[nMuonCand]
   Int_t           MuonCand_tm2dloosemuid[10];   //[nMuonCand]
   Int_t           MuonCand_tmlsAngloosemuonid[10];   //[nMuonCand]
   Int_t           MuonCand_tmlsAngtightmuonid[10];   //[nMuonCand]
   Int_t           MuonCand_tmosAngloosemuonid[10];   //[nMuonCand]
   Int_t           MuonCand_tmosAngtightmuonid[10];   //[nMuonCand]
   Int_t           MuonCand_arbmuid[10];   //[nMuonCand]
   Int_t           MuonCand_gmPromptTight[10];   //[nMuonCand]
   Int_t           MuonCand_isglobal[10];   //[nMuonCand]
   Int_t           MuonCand_istracker[10];   //[nMuonCand]
   Int_t           MuonCand_isstandalone[10];   //[nMuonCand]
   Double_t        MuonCand_hcalisor3[10];   //[nMuonCand]
   Double_t        MuonCand_ecalisor3[10];   //[nMuonCand]
   Double_t        MuonCand_hoisor3[10];   //[nMuonCand]
   Double_t        MuonCand_trkisor3[10];   //[nMuonCand]
   Double_t        MuonCand_hcalisor5[10];   //[nMuonCand]
   Double_t        MuonCand_ecalisor5[10];   //[nMuonCand]
   Double_t        MuonCand_hoisor5[10];   //[nMuonCand]
   Double_t        MuonCand_trkisor5[10];   //[nMuonCand]
   Double_t        MuonCand_timein[10];   //[nMuonCand]
   Double_t        MuonCand_timeout[10];   //[nMuonCand]
   Double_t        MuonCand_timeinerr[10];   //[nMuonCand]
   Double_t        MuonCand_timeouterr[10];   //[nMuonCand]
   Double_t        MuonCand_efficiency[10];   //[nMuonCand]
   Int_t           MuonCand_validtrackhits[10];   //[nMuonCand]
   Int_t           MuonCand_validhits[10];   //[nMuonCand]
   Int_t           MuonCand_validpixelhits[10];   //[nMuonCand]
   Int_t           MuonCand_validmuonhits[10];   //[nMuonCand]
   Int_t           MuonCand_matches[10];   //[nMuonCand]
   Int_t           MuonCand_nlayers[10];   //[nMuonCand]
   Double_t        MuonCand_normchi2[10];   //[nMuonCand]
   Double_t        MuonCand_normtrackchi2[10];   //[nMuonCand]
   Double_t        MuonCand_dB[10];   //[nMuonCand]
   Int_t           MuonCand_tightID[10];   //[nMuonCand]
   Int_t           MuonCand_PF[10];   //[nMuonCand]
   Int_t           MuEPairCand[2];
   Int_t           nEleCand;
   Double_t        EleCand_px[10];   //[nEleCand]
   Double_t        EleCand_py[10];   //[nEleCand]
   Double_t        EleCand_pz[10];   //[nEleCand]
   Double_t        EleCand_p[10];   //[nEleCand]
   Double_t        EleCand_e[10];   //[nEleCand]
   Double_t        EleCand_et[10];   //[nEleCand]
   Double_t        EleCand_eta[10];   //[nEleCand]
   Double_t        EleCand_phi[10];   //[nEleCand]
   Int_t           EleCand_charge[10];   //[nEleCand]
   Int_t           EleCand_looseid[10];   //[nEleCand]
   Double_t        EleCand_likelihoodid[10];   //[nEleCand]
   Int_t           EleCand_robustid[10];   //[nEleCand]
   Double_t        EleCandTrack_p[10];   //[nEleCand]
   Double_t        EleCandTrack_pt[10];   //[nEleCand]
   Double_t        EleCandTrack_eta[10];   //[nEleCand]
   Double_t        EleCandTrack_phi[10];   //[nEleCand]
   Double_t        EleCandTrack_vtxz[10];   //[nEleCand]
   Double_t        EleCand_vtxx[10];   //[nEleCand]
   Double_t        EleCand_vtxy[10];   //[nEleCand]
   Double_t        EleCand_vtxz[10];   //[nEleCand]
   Double_t        EleCand_deltaPhi[10];   //[nEleCand]
   Double_t        EleCand_deltaEta[10];   //[nEleCand]
   Double_t        EleCand_HoverE[10];   //[nEleCand]
   Double_t        EleCand_trackiso[10];   //[nEleCand]
   Double_t        EleCand_ecaliso[10];   //[nEleCand]
   Double_t        EleCand_hcaliso[10];   //[nEleCand]
   Double_t        EleCand_sigmaIetaIeta[10];   //[nEleCand]
   Double_t        EleCand_convDist[10];   //[nEleCand]
   Double_t        EleCand_convDcot[10];   //[nEleCand]
   Double_t        EleCand_ecalDriven[10];   //[nEleCand]
   Int_t           EleCand_wp80[10];   //[nEleCand]
   Int_t           EleCand_mediumID[10];   //[nEleCand]
   Int_t           EleCand_looseID[10];   //[nEleCand]
   Double_t        HLT_Mu10Ele10_MuonCand_pt[10];   //[nHLTMu10Ele10MuonCand]
   Double_t        HLT_Mu10Ele10_MuonCand_eta[10];   //[nHLTMu10Ele10MuonCand]
   Double_t        HLT_Mu10Ele10_MuonCand_phi[10];   //[nHLTMu10Ele10MuonCand]
   Int_t           HLT_Mu10Ele10_MuonCand_charge[10];   //[nHLTMu10Ele10MuonCand]
   Double_t        HLT_Mu8Ele17_MuonCand_pt[10];   //[nHLTMu8Ele17MuonCand]
   Double_t        HLT_Mu8Ele17_MuonCand_eta[10];   //[nHLTMu8Ele17MuonCand]
   Double_t        HLT_Mu8Ele17_MuonCand_phi[10];   //[nHLTMu8Ele17MuonCand]
   Int_t           HLT_Mu8Ele17_MuonCand_charge[10];   //[nHLTMu8Ele17MuonCand]
   Double_t        HLT_Mu17Ele8_MuonCand_pt[10];   //[nHLTMu17Ele8MuonCand]
   Double_t        HLT_Mu17Ele8_MuonCand_eta[10];   //[nHLTMu17Ele8MuonCand]
   Double_t        HLT_Mu17Ele8_MuonCand_phi[10];   //[nHLTMu17Ele8MuonCand]
   Int_t           HLT_Mu17Ele8_MuonCand_charge[10];   //[nHLTMu17Ele8MuonCand]
   Double_t        HLT_Mu8Ele17_EleLCand_pt[10];   //[nHLTMu8Ele17EleLCand]
   Double_t        HLT_Mu8Ele17_EleLCand_eta[10];   //[nHLTMu8Ele17EleLCand]
   Double_t        HLT_Mu8Ele17_EleLCand_phi[10];   //[nHLTMu8Ele17EleLCand]
   Int_t           HLT_Mu8Ele17_EleLCand_charge[10];   //[nHLTMu8Ele17EleLCand]
   Double_t        HLT_Mu17Ele8_EleLCand_pt[10];   //[nHLTMu17Ele8EleLCand]
   Double_t        HLT_Mu17Ele8_EleLCand_eta[10];   //[nHLTMu17Ele8EleLCand]
   Double_t        HLT_Mu17Ele8_EleLCand_phi[10];   //[nHLTMu17Ele8EleLCand]
   Int_t           HLT_Mu17Ele8_EleLCand_charge[10];   //[nHLTMu17Ele8EleLCand]
   Double_t        HLT_Mu8Ele17_EleTCand_pt[10];   //[nHLTMu8Ele17EleTCand]
   Double_t        HLT_Mu8Ele17_EleTCand_eta[10];   //[nHLTMu8Ele17EleTCand]
   Double_t        HLT_Mu8Ele17_EleTCand_phi[10];   //[nHLTMu8Ele17EleTCand]
   Int_t           HLT_Mu8Ele17_EleTCand_charge[10];   //[nHLTMu8Ele17EleTCand]
   Double_t        HLT_Mu17Ele8_EleTCand_pt[10];   //[nHLTMu17Ele8EleTCand]
   Double_t        HLT_Mu17Ele8_EleTCand_eta[10];   //[nHLTMu17Ele8EleTCand]
   Double_t        HLT_Mu17Ele8_EleTCand_phi[10];   //[nHLTMu17Ele8EleTCand]
   Int_t           HLT_Mu17Ele8_EleTCand_charge[10];   //[nHLTMu17Ele8EleTCand]
   Int_t           HLT_Mu17Ele8L;
   Int_t           HLT_Mu8Ele17L;
   Int_t           HLT_Mu17Ele8T;
   Int_t           HLT_Mu8Ele17T;
   Int_t           HLT_Mu10Ele10;
   Int_t           HLT_Mu17Ele8L_Prescl;
   Int_t           HLT_Mu8Ele17L_Prescl;
   Int_t           HLT_Mu17Ele8T_Prescl;
   Int_t           HLT_Mu8Ele17T_Prescl;
   Int_t           HLT_Mu10Ele10_Prescl;
   Int_t           nTrackCand;
   Int_t           nQualityTrackCand;
   Int_t           nExtraTrackCand;
   Double_t        TrackCand_px[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_py[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_pz[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_p[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_pt[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_eta[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_phi[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_vtxdxyz[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_charge[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_purity[1000];   //[nExtraTrackCand]
   Int_t           TrackCand_nhits[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_chi2[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_ndof[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_vtxZ[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_vtxT[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_X[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_Y[1000];   //[nExtraTrackCand]
   Double_t        TrackCand_Z[1000];   //[nExtraTrackCand]
   Double_t        ClosestExtraTrack_vtxdxyz;
   Double_t        ClosestHighPurityExtraTrack_vtxdxyz;
   Double_t        MuE_mass;
   Double_t        MuE_dphi;
   Double_t        MuE_dpt;
   Double_t        MuE_pt;
   Double_t        MuE_phi;
   Double_t        MuE_3Dangle;
   Double_t        MuE_Kalmanvtxx;
   Double_t        MuE_Kalmanvtxy;
   Double_t        MuE_Kalmanvtxz;
   Double_t        MuE_KalmanvtxT;
   Double_t        MuE_Kalmanvtxchi2dof;
   Int_t           MuE_Kalmanvtxisvalid;
   Int_t           MuE_extratracks1mm;
   Int_t           MuE_extratracks2mm;
   Int_t           MuE_extratracks3mm;
   Int_t           MuE_extratracks4mm;
   Int_t           MuE_extratracks5mm;
   Int_t           MuE_extratracks1cm;
   Int_t           MuE_extratracks3cm;
   Int_t           MuE_extratracks5cm;
   Int_t           MuE_extratracks10cm;
   Int_t           nGenMuonCand;
   Double_t        GenMuonCand_px[10];   //[nGenMuonCand]
   Double_t        GenMuonCand_py[10];   //[nGenMuonCand]
   Double_t        GenMuonCand_pz[10];   //[nGenMuonCand]
   Double_t        GenMuonCand_pt[10];   //[nGenMuonCand]
   Double_t        GenMuonCand_eta[10];   //[nGenMuonCand]
   Double_t        GenMuonCand_phi[10];   //[nGenMuonCand]
   Int_t           nGenEleCand;
   Double_t        GenEleCand_px[10];   //[nGenEleCand]
   Double_t        GenEleCand_py[10];   //[nGenEleCand]
   Double_t        GenEleCand_pz[10];   //[nGenEleCand]
   Double_t        GenEleCand_pt[10];   //[nGenEleCand]
   Double_t        GenEleCand_eta[10];   //[nGenEleCand]
   Double_t        GenEleCand_phi[10];   //[nGenEleCand]
   Double_t        GenMuE_eta;
   Double_t        GenMuE_pt;
   Double_t        Etmiss;
   Double_t        Etmiss_phi;
   Double_t        Etmiss_x;
   Double_t        Etmiss_y;
   Double_t        Etmiss_z;
   Double_t        Etmiss_significance;
   Int_t           HLT_Mu10Ele10;
   Int_t           HLT_Mu10Ele10_Prescl;
   Int_t           Run;
   Int_t           LumiSection;
   Int_t           BX;
   Double_t        AvgInstDelLumi;
   Double_t        BunchInstLumi[3];
   Int_t           EventNum;
   Int_t           L1TechnicalTriggers[128];
   Int_t           nPrimVertexCand;
   Double_t        PrimVertexCand_x[60];   //[nPrimVertexCand]
   Double_t        PrimVertexCand_y[60];   //[nPrimVertexCand]
   Double_t        PrimVertexCand_z[60];   //[nPrimVertexCand]
   Int_t           PrimVertexCand_tracks[60];   //[nPrimVertexCand]
   Double_t        PrimVertexCand_chi2[60];   //[nPrimVertexCand]
   Double_t        PrimVertexCand_ndof[60];   //[nPrimVertexCand]
   Int_t           PrimVertexCand_mueTwoTracks[60];   //[nPrimVertexCand]
   Int_t           PrimVertexCand_mueExactlyTwoTracks[60];   //[nPrimVertexCand]
   Int_t           PrimVertexCand_mueTwoTracksMuIndex[60];   //[nPrimVertexCand]
   Int_t           PrimVertexCand_mueTwoTracksEleIndex[60];   //[nPrimVertexCand]
   Int_t           nTruePUforPUWeight;
   Double_t        nTruePUafterPUWeight;
   Int_t           nTruePUforPUWeightBXM1;
   Double_t        nTruePUafterPUWeightBXM1;
   Int_t           nTruePUforPUWeightBXP1;
   Double_t        nTruePUafterPUWeightBXP1;
   Int_t           nTruePUforPUWeightBX0;
   Double_t        nTruePUafterPUWeightBX0;
   Double_t        Weight3D;
   Double_t        PUWeightTrue;

   // List of branches
   TBranch        *b_nJetCand;   //!
   TBranch        *b_JetCand_px;   //!
   TBranch        *b_JetCand_py;   //!
   TBranch        *b_JetCand_pz;   //!
   TBranch        *b_JetCand_e;   //!
   TBranch        *b_JetCand_eta;   //!
   TBranch        *b_JetCand_phi;   //!
   TBranch        *b_HighestJet_e;   //!
   TBranch        *b_HighestJet_eta;   //!
   TBranch        *b_HighestJet_phi;   //!
   TBranch        *b_SumJet_e;   //!
   TBranch        *b_nMuonCand;   //!
   TBranch        *b_MuonCand_px;   //!
   TBranch        *b_MuonCand_py;   //!
   TBranch        *b_MuonCand_pz;   //!
   TBranch        *b_MuonCand_vtxx;   //!
   TBranch        *b_MuonCand_vtxy;   //!
   TBranch        *b_MuonCand_vtxz;   //!
   TBranch        *b_MuonCand_p;   //!
   TBranch        *b_MuonCand_pt;   //!
   TBranch        *b_MuonCand_eta;   //!
   TBranch        *b_MuonCand_phi;   //!
   TBranch        *b_MuonCand_charge;   //!
   TBranch        *b_MuonCand_tmlsloosemuonid;   //!
   TBranch        *b_MuonCand_tmlsOptLowPtloosemuonid;   //!
   TBranch        *b_MuonCand_tm2dloosemuid;   //!
   TBranch        *b_MuonCand_tmlsAngloosemuonid;   //!
   TBranch        *b_MuonCand_tmlsAngtightmuonid;   //!
   TBranch        *b_MuonCand_tmosAngloosemuonid;   //!
   TBranch        *b_MuonCand_tmosAngtightmuonid;   //!
   TBranch        *b_MuonCand_arbmuid;   //!
   TBranch        *b_MuonCand_gmPromptTight;   //!
   TBranch        *b_MuonCand_isglobal;   //!
   TBranch        *b_MuonCand_istracker;   //!
   TBranch        *b_MuonCand_isstandalone;   //!
   TBranch        *b_MuonCand_hcalisor3;   //!
   TBranch        *b_MuonCand_ecalisor3;   //!
   TBranch        *b_MuonCand_hoisor3;   //!
   TBranch        *b_MuonCand_trkisor3;   //!
   TBranch        *b_MuonCand_hcalisor5;   //!
   TBranch        *b_MuonCand_ecalisor5;   //!
   TBranch        *b_MuonCand_hoisor5;   //!
   TBranch        *b_MuonCand_trkisor5;   //!
   TBranch        *b_MuonCand_timein;   //!
   TBranch        *b_MuonCand_timeout;   //!
   TBranch        *b_MuonCand_timeinerr;   //!
   TBranch        *b_MuonCand_timeouterr;   //!
   TBranch        *b_MuonCand_efficiency;   //!
   TBranch        *b_MuonCand_validtrackhits;   //!
   TBranch        *b_MuonCand_validhits;   //!
   TBranch        *b_MuonCand_validpixelhits;   //!
   TBranch        *b_MuonCand_validmuonhits;   //!
   TBranch        *b_MuonCand_matches;   //!
   TBranch        *b_MuonCand_nlayers;   //!
   TBranch        *b_MuonCand_normchi2;   //!
   TBranch        *b_MuonCand_normtrackchi2;   //!
   TBranch        *b_MuonCand_dB;   //!
   TBranch        *b_MuonCand_tightID;   //!
   TBranch        *b_MuonCand_PF;   //!
   TBranch        *b_MuEPairCand;   //!
   TBranch        *b_nEleCand;   //!
   TBranch        *b_EleCand_px;   //!
   TBranch        *b_EleCand_py;   //!
   TBranch        *b_EleCand_pz;   //!
   TBranch        *b_EleCand_p;   //!
   TBranch        *b_EleCand_e;   //!
   TBranch        *b_EleCand_et;   //!
   TBranch        *b_EleCand_eta;   //!
   TBranch        *b_EleCand_phi;   //!
   TBranch        *b_EleCand_charge;   //!
   TBranch        *b_EleCand_looseid;   //!
   TBranch        *b_EleCand_likelihoodid;   //!
   TBranch        *b_EleCand_robustid;   //!
   TBranch        *b_EleCandTrack_p;   //!
   TBranch        *b_EleCandTrack_pt;   //!
   TBranch        *b_EleCandTrack_eta;   //!
   TBranch        *b_EleCandTrack_phi;   //!
   TBranch        *b_EleCandTrack_vtxz;   //!
   TBranch        *b_EleCand_vtxx;   //!
   TBranch        *b_EleCand_vtxy;   //!
   TBranch        *b_EleCand_vtxz;   //!
   TBranch        *b_EleCand_deltaPhi;   //!
   TBranch        *b_EleCand_deltaEta;   //!
   TBranch        *b_EleCand_HoverE;   //!
   TBranch        *b_EleCand_trackiso;   //!
   TBranch        *b_EleCand_ecaliso;   //!
   TBranch        *b_EleCand_hcaliso;   //!
   TBranch        *b_EleCand_sigmaIetaIeta;   //!
   TBranch        *b_EleCand_convDist;   //!
   TBranch        *b_EleCand_convDcot;   //!
   TBranch        *b_EleCand_ecalDriven;   //!
   TBranch        *b_EleCand_wp80;   //!
   TBranch        *b_EleCand_mediumID;   //!
   TBranch        *b_EleCand_looseID;   //! 
   TBranch        *b_HLT_Mu10Ele10_MuonCand_pt;   //!
   TBranch        *b_HLT_Mu10Ele10_MuonCand_eta;   //!
   TBranch        *b_HLT_Mu10Ele10_MuonCand_phi;   //!
   TBranch        *b_HLT_Mu10Ele10_MuonCand_charge;   //!
   TBranch        *b_HLT_Mu8Ele17_MuonCand_pt;   //!
   TBranch        *b_HLT_Mu8Ele17_MuonCand_eta;   //!
   TBranch        *b_HLT_Mu8Ele17_MuonCand_phi;   //!
   TBranch        *b_HLT_Mu8Ele17_MuonCand_charge;   //!
   TBranch        *b_HLT_Mu17Ele8_MuonCand_pt;   //!
   TBranch        *b_HLT_Mu17Ele8_MuonCand_eta;   //!
   TBranch        *b_HLT_Mu17Ele8_MuonCand_phi;   //!
   TBranch        *b_HLT_Mu17Ele8_MuonCand_charge;   //!
   TBranch        *b_HLT_Mu8Ele17_EleLCand_pt;   //!
   TBranch        *b_HLT_Mu8Ele17_EleLCand_eta;   //!
   TBranch        *b_HLT_Mu8Ele17_EleLCand_phi;   //!
   TBranch        *b_HLT_Mu8Ele17_EleLCand_charge;   //!
   TBranch        *b_HLT_Mu17Ele8_EleLCand_pt;   //!
   TBranch        *b_HLT_Mu17Ele8_EleLCand_eta;   //!
   TBranch        *b_HLT_Mu17Ele8_EleLCand_phi;   //!
   TBranch        *b_HLT_Mu17Ele8_EleLCand_charge;   //!
   TBranch        *b_HLT_Mu8Ele17_EleTCand_pt;   //!
   TBranch        *b_HLT_Mu8Ele17_EleTCand_eta;   //!
   TBranch        *b_HLT_Mu8Ele17_EleTCand_phi;   //!
   TBranch        *b_HLT_Mu8Ele17_EleTCand_charge;   //!
   TBranch        *b_HLT_Mu17Ele8_EleTCand_pt;   //!
   TBranch        *b_HLT_Mu17Ele8_EleTCand_eta;   //!
   TBranch        *b_HLT_Mu17Ele8_EleTCand_phi;   //!
   TBranch        *b_HLT_Mu17Ele8_EleTCand_charge;   //!
   TBranch        *b_HLT_Mu17Ele8L;   //!
   TBranch        *b_HLT_Mu8Ele17L;   //!
   TBranch        *b_HLT_Mu17Ele8T;   //!
   TBranch        *b_HLT_Mu8Ele17T;   //!
   TBranch        *b_HLT_Mu10Ele10;   //!
   TBranch        *b_HLT_Mu17Ele8L_Prescl;   //!
   TBranch        *b_HLT_Mu8Ele17L_Prescl;   //!
   TBranch        *b_HLT_Mu17Ele8T_Prescl;   //!
   TBranch        *b_HLT_Mu8Ele17T_Prescl;   //!
   TBranch        *b_HLT_Mu10Ele10_Prescl;   //!
   TBranch        *b_nTrackCand;   //!
   TBranch        *b_nQualityTrackCand;   //!
   TBranch        *b_nExtraTrackCand;   //!
   TBranch        *b_TrackCand_px;   //!
   TBranch        *b_TrackCand_py;   //!
   TBranch        *b_TrackCand_pz;   //!
   TBranch        *b_TrackCand_p;   //!
   TBranch        *b_TrackCand_pt;   //!
   TBranch        *b_TrackCand_eta;   //!
   TBranch        *b_TrackCand_phi;   //!
   TBranch        *b_TrackCand_vtxdxyz;   //!
   TBranch        *b_TrackCand_charge;   //!
   TBranch        *b_TrackCand_purity;   //!
   TBranch        *b_TrackCand_nhits;   //!
   TBranch        *b_TrackCand_chi2;   //!
   TBranch        *b_TrackCand_ndof;   //!
   TBranch        *b_TrackCand_vtxZ;   //!
   TBranch        *b_TrackCand_vtxT;   //!
   TBranch        *b_TrackCand_X;   //!
   TBranch        *b_TrackCand_Y;   //!
   TBranch        *b_TrackCand_Z;   //!
   TBranch        *b_ClosestExtraTrack_vtxdxyz;   //!
   TBranch        *b_ClosestHighPurityExtraTrack_vtxdxyz;   //!
   TBranch        *b_MuE_mass;   //!
   TBranch        *b_MuE_dphi;   //!
   TBranch        *b_MuE_dpt;   //!
   TBranch        *b_MuE_pt;   //!
   TBranch        *b_MuE_phi;   //!
   TBranch        *b_MuE_3Dangle;   //!
   TBranch        *b_MuE_Kalmanvtxx;   //!
   TBranch        *b_MuE_Kalmanvtxy;   //!
   TBranch        *b_MuE_Kalmanvtxz;   //!
   TBranch        *b_MuE_KalmanvtxT;   //!
   TBranch        *b_MuE_Kalmanvtxchi2dof;   //!
   TBranch        *b_MuE_Kalmanvtxisvalid;   //!
   TBranch        *b_MuE_extratracks1mm;   //!
   TBranch        *b_MuE_extratracks2mm;   //!
   TBranch        *b_MuE_extratracks3mm;   //!
   TBranch        *b_MuE_extratracks4mm;   //!
   TBranch        *b_MuE_extratracks5mm;   //!
   TBranch        *b_MuE_extratracks1cm;   //!
   TBranch        *b_MuE_extratracks3cm;   //!
   TBranch        *b_MuE_extratracks5cm;   //!
   TBranch        *b_MuE_extratracks10cm;   //!
   TBranch        *b_nGenMuonCand;   //!
   TBranch        *b_GenMuonCand_px;   //!
   TBranch        *b_GenMuonCand_py;   //!
   TBranch        *b_GenMuonCand_pz;   //!
   TBranch        *b_GenMuonCand_pt;   //!
   TBranch        *b_GenMuonCand_eta;   //!
   TBranch        *b_GenMuonCand_phi;   //!
   TBranch        *b_nGenEleCand;   //!
   TBranch        *b_GenEleCand_px;   //!
   TBranch        *b_GenEleCand_py;   //!
   TBranch        *b_GenEleCand_pz;   //!
   TBranch        *b_GenEleCand_pt;   //!
   TBranch        *b_GenEleCand_eta;   //!
   TBranch        *b_GenEleCand_phi;   //!
   TBranch        *b_GenMuE_eta;   //!
   TBranch        *b_GenMuE_pt;   //!
   TBranch        *b_Etmiss;   //!
   TBranch        *b_Etmiss_phi;   //!
   TBranch        *b_Etmiss_x;   //!
   TBranch        *b_Etmiss_y;   //!
   TBranch        *b_Etmiss_z;   //!
   TBranch        *b_Etmiss_significance;   //!
   TBranch        *b_HLT_Mu10Ele10;   //!
   TBranch        *b_HLT_Mu10Ele10_Prescl;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_LumiSection;   //!
   TBranch        *b_BX;   //!
   TBranch        *b_AvgInstDelLumi;   //!
   TBranch        *b_BunchInstLumi;   //!
   TBranch        *b_EventNum;   //!
   TBranch        *b_L1TechnicalTriggers;   //!
   TBranch        *b_nPrimVertexCand;   //!
   TBranch        *b_PrimVertexCand_x;   //!
   TBranch        *b_PrimVertexCand_y;   //!
   TBranch        *b_PrimVertexCand_z;   //!
   TBranch        *b_PrimVertexCand_tracks;   //!
   TBranch        *b_PrimVertexCand_chi2;   //!
   TBranch        *b_PrimVertexCand_ndof;   //!
   TBranch        *b_PrimVertexCand_mueTwoTracks;   //!
   TBranch        *b_PrimVertexCand_mueExactlyTwoTracks;   //!
   TBranch        *b_PrimVertexCand_mueTwoTracksMuIndex;   //!
   TBranch        *b_PrimVertexCand_mueTwoTracksEleIndex;   //!
   TBranch        *b_nTruePUforPUWeight;   //!
   TBranch        *b_nTruePUafterPUWeight;   //!
   TBranch        *b_nTruePUforPUWeightBXM1;   //!
   TBranch        *b_nTruePUafterPUWeightBXM1;   //!
   TBranch        *b_nTruePUforPUWeightBXP1;   //!
   TBranch        *b_nTruePUafterPUWeightBXP1;   //!
   TBranch        *b_nTruePUforPUWeightBX0;   //!
   TBranch        *b_nTruePUafterPUWeightBX0;   //!
   TBranch        *b_Weight3D;   //!
   TBranch        *b_PUWeightTrue;   //!

   tr1->SetBranchAddress("nMuonCand", &nMuonCand, &b_nMuonCand);
   tr1->SetBranchAddress("MuonCand_px", MuonCand_px, &b_MuonCand_px);
   tr1->SetBranchAddress("MuonCand_py", MuonCand_py, &b_MuonCand_py);
   tr1->SetBranchAddress("MuonCand_pz", MuonCand_pz, &b_MuonCand_pz);
   tr1->SetBranchAddress("MuonCand_vtxx", MuonCand_vtxx, &b_MuonCand_vtxx);
   tr1->SetBranchAddress("MuonCand_vtxy", MuonCand_vtxy, &b_MuonCand_vtxy);
   tr1->SetBranchAddress("MuonCand_vtxz", MuonCand_vtxz, &b_MuonCand_vtxz);
   tr1->SetBranchAddress("MuonCand_p", MuonCand_p, &b_MuonCand_p);
   tr1->SetBranchAddress("MuonCand_pt", MuonCand_pt, &b_MuonCand_pt);
   tr1->SetBranchAddress("MuonCand_eta", MuonCand_eta, &b_MuonCand_eta);
   tr1->SetBranchAddress("MuonCand_phi", MuonCand_phi, &b_MuonCand_phi);
   tr1->SetBranchAddress("MuonCand_charge", MuonCand_charge, &b_MuonCand_charge);
   tr1->SetBranchAddress("MuonCand_tmlsloosemuonid", MuonCand_tmlsloosemuonid, &b_MuonCand_tmlsloosemuonid);
   tr1->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid", MuonCand_tmlsOptLowPtloosemuonid, &b_MuonCand_tmlsOptLowPtloosemuonid);
   tr1->SetBranchAddress("MuonCand_tm2dloosemuid", MuonCand_tm2dloosemuid, &b_MuonCand_tm2dloosemuid);
   tr1->SetBranchAddress("MuonCand_tmlsAngloosemuonid", MuonCand_tmlsAngloosemuonid, &b_MuonCand_tmlsAngloosemuonid);
   tr1->SetBranchAddress("MuonCand_tmlsAngtightmuonid", MuonCand_tmlsAngtightmuonid, &b_MuonCand_tmlsAngtightmuonid);
   tr1->SetBranchAddress("MuonCand_tmosAngloosemuonid", MuonCand_tmosAngloosemuonid, &b_MuonCand_tmosAngloosemuonid);
   tr1->SetBranchAddress("MuonCand_tmosAngtightmuonid", MuonCand_tmosAngtightmuonid, &b_MuonCand_tmosAngtightmuonid);
   tr1->SetBranchAddress("MuonCand_arbmuid", MuonCand_arbmuid, &b_MuonCand_arbmuid);
   tr1->SetBranchAddress("MuonCand_gmPromptTight", MuonCand_gmPromptTight, &b_MuonCand_gmPromptTight);
   tr1->SetBranchAddress("MuonCand_isglobal", MuonCand_isglobal, &b_MuonCand_isglobal);
   tr1->SetBranchAddress("MuonCand_istracker", MuonCand_istracker, &b_MuonCand_istracker);
   tr1->SetBranchAddress("MuonCand_isstandalone", MuonCand_isstandalone, &b_MuonCand_isstandalone);
   tr1->SetBranchAddress("MuonCand_hcalisor3", MuonCand_hcalisor3, &b_MuonCand_hcalisor3);
   tr1->SetBranchAddress("MuonCand_ecalisor3", MuonCand_ecalisor3, &b_MuonCand_ecalisor3);
   tr1->SetBranchAddress("MuonCand_hoisor3", MuonCand_hoisor3, &b_MuonCand_hoisor3);
   tr1->SetBranchAddress("MuonCand_trkisor3", MuonCand_trkisor3, &b_MuonCand_trkisor3);
   tr1->SetBranchAddress("MuonCand_hcalisor5", MuonCand_hcalisor5, &b_MuonCand_hcalisor5);
   tr1->SetBranchAddress("MuonCand_ecalisor5", MuonCand_ecalisor5, &b_MuonCand_ecalisor5);
   tr1->SetBranchAddress("MuonCand_hoisor5", MuonCand_hoisor5, &b_MuonCand_hoisor5);
   tr1->SetBranchAddress("MuonCand_trkisor5", MuonCand_trkisor5, &b_MuonCand_trkisor5);
   tr1->SetBranchAddress("MuonCand_timein", MuonCand_timein, &b_MuonCand_timein);
   tr1->SetBranchAddress("MuonCand_timeout", MuonCand_timeout, &b_MuonCand_timeout);
   tr1->SetBranchAddress("MuonCand_timeinerr", MuonCand_timeinerr, &b_MuonCand_timeinerr);
   tr1->SetBranchAddress("MuonCand_timeouterr", MuonCand_timeouterr, &b_MuonCand_timeouterr);
   tr1->SetBranchAddress("MuonCand_efficiency", MuonCand_efficiency, &b_MuonCand_efficiency);
   tr1->SetBranchAddress("MuonCand_validtrackhits", MuonCand_validtrackhits, &b_MuonCand_validtrackhits);
   tr1->SetBranchAddress("MuonCand_validhits", MuonCand_validhits, &b_MuonCand_validhits);
   tr1->SetBranchAddress("MuonCand_validpixelhits", MuonCand_validpixelhits, &b_MuonCand_validpixelhits);
   tr1->SetBranchAddress("MuonCand_validmuonhits", MuonCand_validmuonhits, &b_MuonCand_validmuonhits);
   tr1->SetBranchAddress("MuonCand_matches", MuonCand_matches, &b_MuonCand_matches);
   tr1->SetBranchAddress("MuonCand_nlayers", MuonCand_nlayers, &b_MuonCand_nlayers);
   tr1->SetBranchAddress("MuonCand_normchi2", MuonCand_normchi2, &b_MuonCand_normchi2);
   tr1->SetBranchAddress("MuonCand_normtrackchi2", MuonCand_normtrackchi2, &b_MuonCand_normtrackchi2);
   tr1->SetBranchAddress("MuonCand_dB", MuonCand_dB, &b_MuonCand_dB);
   tr1->SetBranchAddress("MuonCand_tightID", MuonCand_tightID, &b_MuonCand_tightID);
   tr1->SetBranchAddress("MuonCand_PF", MuonCand_PF, &b_MuonCand_PF);
   tr1->SetBranchAddress("MuEPairCand", MuEPairCand, &b_MuEPairCand);
   tr1->SetBranchAddress("nEleCand", &nEleCand, &b_nEleCand);
   tr1->SetBranchAddress("EleCand_px", EleCand_px, &b_EleCand_px);
   tr1->SetBranchAddress("EleCand_py", EleCand_py, &b_EleCand_py);
   tr1->SetBranchAddress("EleCand_pz", EleCand_pz, &b_EleCand_pz);
   tr1->SetBranchAddress("EleCand_p", EleCand_p, &b_EleCand_p);
   tr1->SetBranchAddress("EleCand_e", EleCand_e, &b_EleCand_e);
   tr1->SetBranchAddress("EleCand_et", EleCand_et, &b_EleCand_et);
   tr1->SetBranchAddress("EleCand_eta", EleCand_eta, &b_EleCand_eta);
   tr1->SetBranchAddress("EleCand_phi", EleCand_phi, &b_EleCand_phi);
   tr1->SetBranchAddress("EleCand_charge", EleCand_charge, &b_EleCand_charge);
   tr1->SetBranchAddress("EleCand_looseid", EleCand_looseid, &b_EleCand_looseid);
   tr1->SetBranchAddress("EleCand_likelihoodid", EleCand_likelihoodid, &b_EleCand_likelihoodid);
   tr1->SetBranchAddress("EleCand_robustid", EleCand_robustid, &b_EleCand_robustid);
   tr1->SetBranchAddress("EleCandTrack_p", EleCandTrack_p, &b_EleCandTrack_p);
   tr1->SetBranchAddress("EleCandTrack_pt", EleCandTrack_pt, &b_EleCandTrack_pt);
   tr1->SetBranchAddress("EleCandTrack_eta", EleCandTrack_eta, &b_EleCandTrack_eta);
   tr1->SetBranchAddress("EleCandTrack_phi", EleCandTrack_phi, &b_EleCandTrack_phi);
   tr1->SetBranchAddress("EleCandTrack_vtxz", EleCandTrack_vtxz, &b_EleCandTrack_vtxz);
   tr1->SetBranchAddress("EleCand_vtxx", EleCand_vtxx, &b_EleCand_vtxx);
   tr1->SetBranchAddress("EleCand_vtxy", EleCand_vtxy, &b_EleCand_vtxy);
   tr1->SetBranchAddress("EleCand_vtxz", EleCand_vtxz, &b_EleCand_vtxz);
   tr1->SetBranchAddress("EleCand_deltaPhi", EleCand_deltaPhi, &b_EleCand_deltaPhi);
   tr1->SetBranchAddress("EleCand_deltaEta", EleCand_deltaEta, &b_EleCand_deltaEta);
   tr1->SetBranchAddress("EleCand_HoverE", EleCand_HoverE, &b_EleCand_HoverE);
   tr1->SetBranchAddress("EleCand_trackiso", EleCand_trackiso, &b_EleCand_trackiso);
   tr1->SetBranchAddress("EleCand_ecaliso", EleCand_ecaliso, &b_EleCand_ecaliso);
   tr1->SetBranchAddress("EleCand_hcaliso", EleCand_hcaliso, &b_EleCand_hcaliso);
   tr1->SetBranchAddress("EleCand_sigmaIetaIeta", EleCand_sigmaIetaIeta, &b_EleCand_sigmaIetaIeta);
   tr1->SetBranchAddress("EleCand_convDist", EleCand_convDist, &b_EleCand_convDist);
   tr1->SetBranchAddress("EleCand_convDcot", EleCand_convDcot, &b_EleCand_convDcot);
   tr1->SetBranchAddress("EleCand_ecalDriven", EleCand_ecalDriven, &b_EleCand_ecalDriven);
   tr1->SetBranchAddress("EleCand_wp80", EleCand_wp80, &b_EleCand_wp80);
   tr1->SetBranchAddress("EleCand_mediumID", EleCand_mediumID, &b_EleCand_mediumID);
   tr1->SetBranchAddress("EleCand_looseID", EleCand_looseID, &b_EleCand_looseID); 
   tr1->SetBranchAddress("HLT_Mu10Ele10_MuonCand_pt", &HLT_Mu10Ele10_MuonCand_pt, &b_HLT_Mu10Ele10_MuonCand_pt);
   tr1->SetBranchAddress("HLT_Mu10Ele10_MuonCand_eta", &HLT_Mu10Ele10_MuonCand_eta, &b_HLT_Mu10Ele10_MuonCand_eta);
   tr1->SetBranchAddress("HLT_Mu10Ele10_MuonCand_phi", &HLT_Mu10Ele10_MuonCand_phi, &b_HLT_Mu10Ele10_MuonCand_phi);
   tr1->SetBranchAddress("HLT_Mu10Ele10_MuonCand_charge", &HLT_Mu10Ele10_MuonCand_charge, &b_HLT_Mu10Ele10_MuonCand_charge);
   tr1->SetBranchAddress("HLT_Mu8Ele17_MuonCand_pt", HLT_Mu8Ele17_MuonCand_pt, &b_HLT_Mu8Ele17_MuonCand_pt);
   tr1->SetBranchAddress("HLT_Mu8Ele17_MuonCand_eta", HLT_Mu8Ele17_MuonCand_eta, &b_HLT_Mu8Ele17_MuonCand_eta);
   tr1->SetBranchAddress("HLT_Mu8Ele17_MuonCand_phi", HLT_Mu8Ele17_MuonCand_phi, &b_HLT_Mu8Ele17_MuonCand_phi);
   tr1->SetBranchAddress("HLT_Mu8Ele17_MuonCand_charge", HLT_Mu8Ele17_MuonCand_charge, &b_HLT_Mu8Ele17_MuonCand_charge);
   tr1->SetBranchAddress("HLT_Mu17Ele8_MuonCand_pt", HLT_Mu17Ele8_MuonCand_pt, &b_HLT_Mu17Ele8_MuonCand_pt);
   tr1->SetBranchAddress("HLT_Mu17Ele8_MuonCand_eta", HLT_Mu17Ele8_MuonCand_eta, &b_HLT_Mu17Ele8_MuonCand_eta);
   tr1->SetBranchAddress("HLT_Mu17Ele8_MuonCand_phi", HLT_Mu17Ele8_MuonCand_phi, &b_HLT_Mu17Ele8_MuonCand_phi);
   tr1->SetBranchAddress("HLT_Mu17Ele8_MuonCand_charge", HLT_Mu17Ele8_MuonCand_charge, &b_HLT_Mu17Ele8_MuonCand_charge);
   tr1->SetBranchAddress("HLT_Mu8Ele17_EleLCand_pt", HLT_Mu8Ele17_EleLCand_pt, &b_HLT_Mu8Ele17_EleLCand_pt);
   tr1->SetBranchAddress("HLT_Mu8Ele17_EleLCand_eta", HLT_Mu8Ele17_EleLCand_eta, &b_HLT_Mu8Ele17_EleLCand_eta);
   tr1->SetBranchAddress("HLT_Mu8Ele17_EleLCand_phi", HLT_Mu8Ele17_EleLCand_phi, &b_HLT_Mu8Ele17_EleLCand_phi);
   tr1->SetBranchAddress("HLT_Mu8Ele17_EleLCand_charge", HLT_Mu8Ele17_EleLCand_charge, &b_HLT_Mu8Ele17_EleLCand_charge);
   tr1->SetBranchAddress("HLT_Mu17Ele8_EleLCand_pt", HLT_Mu17Ele8_EleLCand_pt, &b_HLT_Mu17Ele8_EleLCand_pt);
   tr1->SetBranchAddress("HLT_Mu17Ele8_EleLCand_eta", HLT_Mu17Ele8_EleLCand_eta, &b_HLT_Mu17Ele8_EleLCand_eta);
   tr1->SetBranchAddress("HLT_Mu17Ele8_EleLCand_phi", HLT_Mu17Ele8_EleLCand_phi, &b_HLT_Mu17Ele8_EleLCand_phi);
   tr1->SetBranchAddress("HLT_Mu17Ele8_EleLCand_charge", HLT_Mu17Ele8_EleLCand_charge, &b_HLT_Mu17Ele8_EleLCand_charge);
   tr1->SetBranchAddress("HLT_Mu8Ele17_EleTCand_pt", HLT_Mu8Ele17_EleTCand_pt, &b_HLT_Mu8Ele17_EleTCand_pt);
   tr1->SetBranchAddress("HLT_Mu8Ele17_EleTCand_eta", HLT_Mu8Ele17_EleTCand_eta, &b_HLT_Mu8Ele17_EleTCand_eta);
   tr1->SetBranchAddress("HLT_Mu8Ele17_EleTCand_phi", HLT_Mu8Ele17_EleTCand_phi, &b_HLT_Mu8Ele17_EleTCand_phi);
   tr1->SetBranchAddress("HLT_Mu8Ele17_EleTCand_charge", HLT_Mu8Ele17_EleTCand_charge, &b_HLT_Mu8Ele17_EleTCand_charge);
   tr1->SetBranchAddress("HLT_Mu17Ele8_EleTCand_pt", HLT_Mu17Ele8_EleTCand_pt, &b_HLT_Mu17Ele8_EleTCand_pt);
   tr1->SetBranchAddress("HLT_Mu17Ele8_EleTCand_eta", HLT_Mu17Ele8_EleTCand_eta, &b_HLT_Mu17Ele8_EleTCand_eta);
   tr1->SetBranchAddress("HLT_Mu17Ele8_EleTCand_phi", HLT_Mu17Ele8_EleTCand_phi, &b_HLT_Mu17Ele8_EleTCand_phi);
   tr1->SetBranchAddress("HLT_Mu17Ele8_EleTCand_charge", HLT_Mu17Ele8_EleTCand_charge, &b_HLT_Mu17Ele8_EleTCand_charge);
   tr1->SetBranchAddress("HLT_Mu17Ele8L", &HLT_Mu17Ele8L, &b_HLT_Mu17Ele8L);
   tr1->SetBranchAddress("HLT_Mu8Ele17L", &HLT_Mu8Ele17L, &b_HLT_Mu8Ele17L);
   tr1->SetBranchAddress("HLT_Mu17Ele8T", &HLT_Mu17Ele8T, &b_HLT_Mu17Ele8T);
   tr1->SetBranchAddress("HLT_Mu8Ele17T", &HLT_Mu8Ele17T, &b_HLT_Mu8Ele17T);
   tr1->SetBranchAddress("HLT_Mu10Ele10", &HLT_Mu10Ele10, &b_HLT_Mu10Ele10);
   tr1->SetBranchAddress("HLT_Mu17Ele8L_Prescl", &HLT_Mu17Ele8L_Prescl, &b_HLT_Mu17Ele8L_Prescl);
   tr1->SetBranchAddress("HLT_Mu8Ele17L_Prescl", &HLT_Mu8Ele17L_Prescl, &b_HLT_Mu8Ele17L_Prescl);
   tr1->SetBranchAddress("HLT_Mu17Ele8T_Prescl", &HLT_Mu17Ele8T_Prescl, &b_HLT_Mu17Ele8T_Prescl);
   tr1->SetBranchAddress("HLT_Mu8Ele17T_Prescl", &HLT_Mu8Ele17T_Prescl, &b_HLT_Mu8Ele17T_Prescl);
   tr1->SetBranchAddress("HLT_Mu10Ele10_Prescl", &HLT_Mu10Ele10_Prescl, &b_HLT_Mu10Ele10_Prescl);
   tr1->SetBranchAddress("nTrackCand", &nTrackCand, &b_nTrackCand);
   tr1->SetBranchAddress("nQualityTrackCand", &nQualityTrackCand, &b_nQualityTrackCand);
   tr1->SetBranchAddress("nExtraTrackCand", &nExtraTrackCand, &b_nExtraTrackCand);
   tr1->SetBranchAddress("TrackCand_px", TrackCand_px, &b_TrackCand_px);
   tr1->SetBranchAddress("TrackCand_py", TrackCand_py, &b_TrackCand_py);
   tr1->SetBranchAddress("TrackCand_pz", TrackCand_pz, &b_TrackCand_pz);
   tr1->SetBranchAddress("TrackCand_p", TrackCand_p, &b_TrackCand_p);
   tr1->SetBranchAddress("TrackCand_pt", TrackCand_pt, &b_TrackCand_pt);
   tr1->SetBranchAddress("TrackCand_eta", TrackCand_eta, &b_TrackCand_eta);
   tr1->SetBranchAddress("TrackCand_phi", TrackCand_phi, &b_TrackCand_phi);
   tr1->SetBranchAddress("TrackCand_vtxdxyz", TrackCand_vtxdxyz, &b_TrackCand_vtxdxyz);
   tr1->SetBranchAddress("TrackCand_charge", TrackCand_charge, &b_TrackCand_charge);
   tr1->SetBranchAddress("TrackCand_purity", TrackCand_purity, &b_TrackCand_purity);
   tr1->SetBranchAddress("TrackCand_nhits", TrackCand_nhits, &b_TrackCand_nhits);
   tr1->SetBranchAddress("TrackCand_chi2", TrackCand_chi2, &b_TrackCand_chi2);
   tr1->SetBranchAddress("TrackCand_ndof", TrackCand_ndof, &b_TrackCand_ndof);
   tr1->SetBranchAddress("TrackCand_vtxZ", TrackCand_vtxZ, &b_TrackCand_vtxZ);
   tr1->SetBranchAddress("TrackCand_vtxT", TrackCand_vtxT, &b_TrackCand_vtxT);
   tr1->SetBranchAddress("TrackCand_X", TrackCand_X, &b_TrackCand_X);
   tr1->SetBranchAddress("TrackCand_Y", TrackCand_Y, &b_TrackCand_Y);
   tr1->SetBranchAddress("TrackCand_Z", TrackCand_Z, &b_TrackCand_Z);
   tr1->SetBranchAddress("MuE_mass", &MuE_mass, &b_MuE_mass);
   tr1->SetBranchAddress("MuE_dphi", &MuE_dphi, &b_MuE_dphi);
   tr1->SetBranchAddress("MuE_dpt", &MuE_dpt, &b_MuE_dpt);
   tr1->SetBranchAddress("MuE_pt", &MuE_pt, &b_MuE_pt);
   tr1->SetBranchAddress("MuE_phi", &MuE_phi, &b_MuE_phi);
   tr1->SetBranchAddress("MuE_3Dangle", &MuE_3Dangle, &b_MuE_3Dangle);
   tr1->SetBranchAddress("MuE_Kalmanvtxx", &MuE_Kalmanvtxx, &b_MuE_Kalmanvtxx);
   tr1->SetBranchAddress("MuE_Kalmanvtxy", &MuE_Kalmanvtxy, &b_MuE_Kalmanvtxy);
   tr1->SetBranchAddress("MuE_Kalmanvtxz", &MuE_Kalmanvtxz, &b_MuE_Kalmanvtxz);
   tr1->SetBranchAddress("MuE_KalmanvtxT", &MuE_KalmanvtxT, &b_MuE_KalmanvtxT);
   tr1->SetBranchAddress("MuE_Kalmanvtxchi2dof", &MuE_Kalmanvtxchi2dof, &b_MuE_Kalmanvtxchi2dof);
   tr1->SetBranchAddress("MuE_Kalmanvtxisvalid", &MuE_Kalmanvtxisvalid, &b_MuE_Kalmanvtxisvalid);
   tr1->SetBranchAddress("MuE_extratracks1mm", &MuE_extratracks1mm, &b_MuE_extratracks1mm);
   tr1->SetBranchAddress("MuE_extratracks2mm", &MuE_extratracks2mm, &b_MuE_extratracks2mm);
   tr1->SetBranchAddress("MuE_extratracks3mm", &MuE_extratracks3mm, &b_MuE_extratracks3mm);
   tr1->SetBranchAddress("MuE_extratracks4mm", &MuE_extratracks4mm, &b_MuE_extratracks4mm);
   tr1->SetBranchAddress("MuE_extratracks5mm", &MuE_extratracks5mm, &b_MuE_extratracks5mm);
   tr1->SetBranchAddress("MuE_extratracks1cm", &MuE_extratracks1cm, &b_MuE_extratracks1cm);
   tr1->SetBranchAddress("MuE_extratracks3cm", &MuE_extratracks3cm, &b_MuE_extratracks3cm);
   tr1->SetBranchAddress("MuE_extratracks5cm", &MuE_extratracks5cm, &b_MuE_extratracks5cm);
   tr1->SetBranchAddress("MuE_extratracks10cm", &MuE_extratracks10cm, &b_MuE_extratracks10cm);
   tr1->SetBranchAddress("nGenMuonCand", &nGenMuonCand, &b_nGenMuonCand);
   tr1->SetBranchAddress("GenMuonCand_px", GenMuonCand_px, &b_GenMuonCand_px);
   tr1->SetBranchAddress("GenMuonCand_py", GenMuonCand_py, &b_GenMuonCand_py);
   tr1->SetBranchAddress("GenMuonCand_pz", GenMuonCand_pz, &b_GenMuonCand_pz);
   tr1->SetBranchAddress("GenMuonCand_pt", GenMuonCand_pt, &b_GenMuonCand_pt);
   tr1->SetBranchAddress("GenMuonCand_eta", GenMuonCand_eta, &b_GenMuonCand_eta);
   tr1->SetBranchAddress("GenMuonCand_phi", GenMuonCand_phi, &b_GenMuonCand_phi);
   tr1->SetBranchAddress("nGenEleCand", &nGenEleCand, &b_nGenEleCand);
   tr1->SetBranchAddress("GenEleCand_px", GenEleCand_px, &b_GenEleCand_px);
   tr1->SetBranchAddress("GenEleCand_py", GenEleCand_py, &b_GenEleCand_py);
   tr1->SetBranchAddress("GenEleCand_pz", GenEleCand_pz, &b_GenEleCand_pz);
   tr1->SetBranchAddress("GenEleCand_pt", GenEleCand_pt, &b_GenEleCand_pt);
   tr1->SetBranchAddress("GenEleCand_eta", GenEleCand_eta, &b_GenEleCand_eta);
   tr1->SetBranchAddress("GenEleCand_phi", GenEleCand_phi, &b_GenEleCand_phi);
   tr1->SetBranchAddress("GenMuE_eta", &GenMuE_eta, &b_GenMuE_eta);
   tr1->SetBranchAddress("GenMuE_pt", &GenMuE_pt, &b_GenMuE_pt);
   tr1->SetBranchAddress("Etmiss", &Etmiss, &b_Etmiss);
   tr1->SetBranchAddress("Etmiss_phi", &Etmiss_phi, &b_Etmiss_phi);
   tr1->SetBranchAddress("Etmiss_x", &Etmiss_x, &b_Etmiss_x);
   tr1->SetBranchAddress("Etmiss_y", &Etmiss_y, &b_Etmiss_y);
   tr1->SetBranchAddress("Etmiss_z", &Etmiss_z, &b_Etmiss_z);
   tr1->SetBranchAddress("Etmiss_significance", &Etmiss_significance, &b_Etmiss_significance);
   tr1->SetBranchAddress("HLT_Mu10Ele10", &HLT_Mu10Ele10, &b_HLT_Mu10Ele10);
   tr1->SetBranchAddress("HLT_Mu10Ele10_Prescl", &HLT_Mu10Ele10_Prescl, &b_HLT_Mu10Ele10_Prescl);
   tr1->SetBranchAddress("Run", &Run, &b_Run);
   tr1->SetBranchAddress("LumiSection", &LumiSection, &b_LumiSection);
   tr1->SetBranchAddress("BX", &BX, &b_BX);
   tr1->SetBranchAddress("AvgInstDelLumi", &AvgInstDelLumi, &b_AvgInstDelLumi);
   tr1->SetBranchAddress("BunchInstLumi", BunchInstLumi, &b_BunchInstLumi);
   tr1->SetBranchAddress("EventNum", &EventNum, &b_EventNum);
   tr1->SetBranchAddress("L1TechnicalTriggers", L1TechnicalTriggers, &b_L1TechnicalTriggers);
   tr1->SetBranchAddress("nPrimVertexCand", &nPrimVertexCand, &b_nPrimVertexCand);
   tr1->SetBranchAddress("PrimVertexCand_x", PrimVertexCand_x, &b_PrimVertexCand_x);
   tr1->SetBranchAddress("PrimVertexCand_y", PrimVertexCand_y, &b_PrimVertexCand_y);
   tr1->SetBranchAddress("PrimVertexCand_z", PrimVertexCand_z, &b_PrimVertexCand_z);
   tr1->SetBranchAddress("PrimVertexCand_tracks", PrimVertexCand_tracks, &b_PrimVertexCand_tracks);
   tr1->SetBranchAddress("PrimVertexCand_chi2", PrimVertexCand_chi2, &b_PrimVertexCand_chi2);
   tr1->SetBranchAddress("PrimVertexCand_ndof", PrimVertexCand_ndof, &b_PrimVertexCand_ndof);
   tr1->SetBranchAddress("PrimVertexCand_mueTwoTracks", PrimVertexCand_mueTwoTracks, &b_PrimVertexCand_mueTwoTracks);
   tr1->SetBranchAddress("PrimVertexCand_mueExactlyTwoTracks", PrimVertexCand_mueExactlyTwoTracks, &b_PrimVertexCand_mueExactlyTwoTracks);
   tr1->SetBranchAddress("PrimVertexCand_mueTwoTracksMuIndex", PrimVertexCand_mueTwoTracksMuIndex, &b_PrimVertexCand_mueTwoTracksMuIndex);
   tr1->SetBranchAddress("PrimVertexCand_mueTwoTracksEleIndex", PrimVertexCand_mueTwoTracksEleIndex, &b_PrimVertexCand_mueTwoTracksEleIndex);
   tr1->SetBranchAddress("nTruePUforPUWeight", &nTruePUforPUWeight, &b_nTruePUforPUWeight);
   tr1->SetBranchAddress("nTruePUafterPUWeight", &nTruePUafterPUWeight, &b_nTruePUafterPUWeight);
   tr1->SetBranchAddress("nTruePUforPUWeightBXM1", &nTruePUforPUWeightBXM1, &b_nTruePUforPUWeightBXM1);
   tr1->SetBranchAddress("nTruePUafterPUWeightBXM1", &nTruePUafterPUWeightBXM1, &b_nTruePUafterPUWeightBXM1);
   tr1->SetBranchAddress("nTruePUforPUWeightBXP1", &nTruePUforPUWeightBXP1, &b_nTruePUforPUWeightBXP1);
   tr1->SetBranchAddress("nTruePUafterPUWeightBXP1", &nTruePUafterPUWeightBXP1, &b_nTruePUafterPUWeightBXP1);
   tr1->SetBranchAddress("nTruePUforPUWeightBX0", &nTruePUforPUWeightBX0, &b_nTruePUforPUWeightBX0);
   tr1->SetBranchAddress("nTruePUafterPUWeightBX0", &nTruePUafterPUWeightBX0, &b_nTruePUafterPUWeightBX0);
   tr1->SetBranchAddress("Weight3D", &Weight3D, &b_Weight3D);
   tr1->SetBranchAddress("PUWeightTrue", &PUWeightTrue, &b_PUWeightTrue);

   nent = tr1->GetEntries();

   // Fake run-dependence - split into 2011A and 2011B
   Double_t runAboundary = 0.457*nent;

   for(Int_t i = 0;i < nent;i++)
    {
      tr1->GetEntry(i);
      if(i%10000 == 0)
	cout << "\t" << i << "/" << nent << endl;
      Double_t theweight = 1.0;
      if(physsample>1)
	theweight = Weight3D;

      // Acceptance
      /*
      if(physsample > 1)
      	{
      	  if(GenMuE_pt <= 100)
      	    continue;
      	  if((fabs(GenMuonCand_eta[0])>=2.4) || (fabs(GenEleCand_eta[0])>=2.4))
      	    continue;
      	  if((GenMuonCand_pt <= 20) || (GenEleCand_pt <= 20))
      	    continue;
      	}
      */

      Int_t mu17ele8 = 0;
      Int_t mu8ele17 = 0;
      if((HLT_Mu17Ele8L && HLT_Mu17Ele8L_Prescl == 1) || (HLT_Mu17Ele8T && HLT_Mu17Ele8T_Prescl == 1))
	mu17ele8 = 1;
      if((HLT_Mu8Ele17L && HLT_Mu8Ele17L_Prescl == 1) || (HLT_Mu8Ele17T && HLT_Mu8Ele17T_Prescl == 1)) 
	mu8ele17 = 1;
      
      if((HLT_Mu17Ele8L && HLT_Mu17Ele8L_Prescl == 1) || (HLT_Mu17Ele8T && HLT_Mu17Ele8T_Prescl == 1) ||
      	 (HLT_Mu8Ele17L && HLT_Mu8Ele17L_Prescl == 1) || (HLT_Mu8Ele17T && HLT_Mu8Ele17T_Prescl == 1))
	{
	  for(Int_t j = 0;j< nPrimVertexCand;j++)
	    {
	      if(PrimVertexCand_tracks[j]-2 < 16 && fabs(PrimVertexCand_z[j]) < 24)
		{
		  if(PrimVertexCand_mueTwoTracks[j] == 1)
		    {
		      int muindex = PrimVertexCand_mueTwoTracksMuIndex[j];
		      int eleindex = PrimVertexCand_mueTwoTracksEleIndex[j];

		      if(((MuonCand_pt[muindex]>20 && EleCand_et[eleindex]>20 && mu8ele17==1) || 
			  (MuonCand_pt[muindex]>20 && EleCand_et[eleindex]>20 && mu17ele8==1)) && 
			 (fabs(MuonCand_eta[muindex])<2.4 && fabs(EleCand_eta[eleindex])<2.4))
			{
			  if((MuonCand_tightID[muindex] == 1) && (EleCand_mediumID[eleindex] == 1))
			    //			  if((MuonCand_tightID[muindex] == 0) || (EleCand_mediumID[eleindex] == 0))
			    {
			      if((MuE_mass > massemumin) && (MuE_pt > pTemumin) && (MuE_pt < pTemumax))
				{
				  // Apply muon ID efficiency corrections
				  Double_t muoneffcorr = 1.0;
				  Double_t mutrigeffcorr = 1.0;
				  Double_t eleeffcorr = 1.0;
				  Double_t eletrigeffcorr = 1.0;

				  if(physsample>1) 
				    {
				      // Muon offline
				      Double_t absmueta = fabs(MuonCand_eta[muindex]);
				      Double_t mupt = MuonCand_pt[muindex];
				      Double_t mueta = MuonCand_eta[muindex];
				      Double_t abseleeta = fabs(EleCand_eta[eleindex]);
				      Double_t elept = EleCand_et[eleindex];
				      for(Int_t x = 0; x < 28;x++)
					{
					  if(i < runAboundary) 
					    {
					      if(efftable[x][7] == 1)
						if((mupt > efftable[x][2]) && (mupt < efftable[x][3])) 
						  if((absmueta > efftable[x][0]) && (absmueta < efftable[x][1])) 
						    muoneffcorr = efftable[x][4];
					    }
					  else
					    {
                                              if(efftable[x][7] == 2) 
                                                if((mupt > efftable[x][2]) && (mupt < efftable[x][3]))  
                                                  if((absmueta > efftable[x][0]) && (absmueta < efftable[x][1]))  
                                                    muoneffcorr = efftable[x][4];
					    }
					}
				      // Muon HLT
				      for(Int_t x = 0; x < 17; x++)
					{
                                          if(mu17ele8==1) 
                                            { 
                                              if((mupt > effhltmu17table[x][2]) && (mupt < effhltmu17table[x][3]))   
                                                if((mueta > effhltmu17table[x][0]) && (mueta < effhltmu17table[x][1]))  
                                                  mutrigeffcorr = effhltmu17table[x][4];
                                            } 
					  else if(mu8ele17==1)
					    {
					      if((mupt > effhltmu8table[x][2]) && (mupt < effhltmu8table[x][3])) 
						if((mueta > effhltmu8table[x][0]) && (mueta < effhltmu8table[x][1]))
						  mutrigeffcorr = effhltmu8table[x][4];
					    }
					}
				      // Electron offline
				      for(Int_t x = 0; x < 30; x++)
					{
					  if((elept > effeletable[x][2]) && (elept < effeletable[x][3]))
					    if((abseleeta > effeletable[x][0]) && (abseleeta < effeletable[x][1]))
					      eleeffcorr = effeletable[x][4];
					}

				      // Electron HLT
				      eletrigeffcorr *= 1.0;
				      
				      theweight *= muoneffcorr;
				      theweight *= mutrigeffcorr;
				      theweight *= eleeffcorr;
				      theweight *= eletrigeffcorr;
				    }
				  
				  // JH
				  Int_t eleismu = 0;
				  for(Int_t z=0;z<nMuonCand;z++)
				    {
				      if(sqrt((EleCand_phi[eleindex]-MuonCand_phi[z])*(EleCand_phi[eleindex]-MuonCand_phi[z]) 
				  	      + (EleCand_eta[eleindex]-MuonCand_eta[z])*(EleCand_eta[eleindex]-MuonCand_eta[z])) < 0.01)
				  	if(MuonCand_tightID[z] == 1)
				  	  eleismu=1;
				    }
				  
				  if(eleismu == 1)
				    continue;
				  // JH

				  hntrk->Fill(PrimVertexCand_tracks[j]-2-0.5,theweight);
				  hptntrackcorr->Fill(PrimVertexCand_tracks[j]-2,MuE_pt,theweight);

				  if((PrimVertexCand_tracks[j]-2) == 0) 
				    {
                                      cout << "Run:LS:Event = " << Run << ":" << LumiSection << ":" << EventNum << endl; 
                                      cout << "pT(mu) = " << MuonCand_pt[muindex] << ", ET(e) = " << EleCand_et[eleindex]  
                                           << ",eta(mu) = " << MuonCand_eta[muindex] << ", eta(e) = " << EleCand_eta[eleindex]  
                                           << endl << "m(emu) = " << MuE_mass << ", dphi(emu) = " << 1 - fabs(MuE_dphi)/3.14159  
                                           << endl << "MET = " << Etmiss << ", pT(emu) = " << MuE_pt  
                                           << endl; 
				      
				      /*
				      cout << "Run:LS:Event = " << Run << ":" << LumiSection << ":" << EventNum << endl;
				      cout << "pT(mu) = " << MuonCand_pt[muindex] << ", ET(e) = " << EleCand_et[eleindex] 
					   << ",eta(mu) = " << MuonCand_eta[muindex] << ", eta(e) = " << EleCand_eta[eleindex] 
					   << endl << "m(emu) = " << MuE_mass << ", dphi(emu) = " << 1 - fabs(MuE_dphi)/3.14159 
					   << endl << "MET = " << Etmiss << ", pT(emu) = " << MuE_pt 
					   << endl;
				      */
				    }

				  if(((PrimVertexCand_tracks[j]-2) >= extratracksmin) && 
				     ((PrimVertexCand_tracks[j]-2) <= extratracksmax))
				    {

				      hmll->Fill(MuE_mass,theweight);
				      hpt->Fill(MuE_pt,theweight); 
				      hdpt->Fill(MuE_dpt,theweight);

				      hdphi->Fill(1 - fabs(MuE_dphi)/3.14159,theweight); 
				      hnvrt->Fill(nPrimVertexCand,theweight); 
				      hmet->Fill(Etmiss,theweight);
				      hvtxz->Fill(PrimVertexCand_z[j],theweight);
				      hvtxT->Fill(sqrt(PrimVertexCand_x[j]*PrimVertexCand_x[j])+(PrimVertexCand_y[j]*PrimVertexCand_y[j]),theweight);
				      hchi2pv->Fill(PrimVertexCand_chi2[j]/PrimVertexCand_ndof[j],theweight);
                                      Double_t muemetdphi = MuE_phi-Etmiss_phi; 
				      if(muemetdphi > 3.14159)
					muemetdphi = (2.0*3.14159)-muemetdphi;
				      hmuemetdphi->Fill(1-(muemetdphi/3.14159),theweight);
				      hemudist->Fill(MuonCand_vtxx[muindex] - EleCand_vtxx[eleindex],theweight);
				      hmupt->Fill(MuonCand_pt[muindex],theweight);
                                      hept->Fill(EleCand_et[eleindex],theweight); 
                                      hmueta->Fill(MuonCand_eta[muindex],theweight); 
                                      heeta->Fill(EleCand_eta[eleindex],theweight); 
				      //				      hptntrackcorr->Fill(PrimVertexCand_tracks[j]-2,MuE_pt,theweight);

				      if(MuonCand_pt[muindex] > EleCand_et[eleindex])
					{
					  hleadpt->Fill(MuonCand_pt[muindex],theweight);
					  htrailpt->Fill(EleCand_et[eleindex],theweight);
					}
				      else
					{
					  hleadpt->Fill(EleCand_et[eleindex],theweight);
					  htrailpt->Fill(MuonCand_pt[muindex],theweight);
					}

				      Double_t sumpttracks = 0.0;
				      Double_t sumpztracks = 0.0;
				      Double_t sumptracks = 0.0;
				      for(Int_t k = 0; k < nExtraTrackCand; k++)
					{
					  hptextra->Fill(TrackCand_pt[k],theweight);
					  hetaextra->Fill(TrackCand_eta[k],theweight);
                                          hdzextra->Fill(TrackCand_vtxZ[k],theweight); 
                                          hdxyextra->Fill(TrackCand_vtxT[k],theweight); 
                                          hchi2extra->Fill(TrackCand_chi2[k]/TrackCand_ndof[k],theweight); 
                                          hhitsextra->Fill(TrackCand_nhits[k],theweight); 
                                          hpurityextra->Fill(TrackCand_purity[k],theweight); 

					  sumpttracks += TrackCand_pt[k];
					  sumpztracks += TrackCand_pz[k];
					  sumptracks += TrackCand_p[k];
					}
				      hsumptextra->Fill(sumpttracks,theweight);
				      hsumpzextra->Fill(sumptracks-sumpztracks,theweight);
				    }
				}
			    }
			}
		    }
		} 
	    }
	}
    }

   cout << "Returning from " << st << endl;

   hmll->SetFillColor(fillcolor); 
   hdphi->SetFillColor(fillcolor); 
   hdpt->SetFillColor(fillcolor); 
   hpt->SetFillColor(fillcolor); 
   hmet->SetFillColor(fillcolor); 
   hntrk->SetFillColor(fillcolor); 
   hnvrt->SetFillColor(fillcolor);  
   hptextra->SetFillColor(fillcolor); 
   hetaextra->SetFillColor(fillcolor); 
   hvtxz->SetFillColor(fillcolor);
   hvtxT->SetFillColor(fillcolor);
   hmuemetdphi->SetFillColor(fillcolor);
   hemudist->SetFillColor(fillcolor); 
   hsumptextra->SetFillColor(fillcolor);
   hsumpzextra->SetFillColor(fillcolor); 
   hdzextra->SetFillColor(fillcolor);  
   hdxyextra->SetFillColor(fillcolor);  
   hchi2extra->SetFillColor(fillcolor); 
   hhitsextra->SetFillColor(fillcolor); 
   hpurityextra->SetFillColor(fillcolor); 
   hmupt->SetFillColor(fillcolor);
   hmueta->SetFillColor(fillcolor); 
   hept->SetFillColor(fillcolor); 
   heeta->SetFillColor(fillcolor); 
   hleadpt->SetFillColor(fillcolor);
   htrailpt->SetFillColor(fillcolor);
   hchi2pv->SetFillColor(fillcolor);

   hmll->SetLineColor(linecolor);  
   hdphi->SetLineColor(linecolor);  
   hdpt->SetLineColor(linecolor);  
   hpt->SetLineColor(linecolor);  
   hmet->SetLineColor(linecolor);  
   hntrk->SetLineColor(linecolor);  
   hnvrt->SetLineColor(linecolor);   
   hptextra->SetLineColor(linecolor);  
   hetaextra->SetLineColor(linecolor);  
   hvtxz->SetLineColor(linecolor); 
   hvtxT->SetLineColor(linecolor); 
   hmuemetdphi->SetLineColor(linecolor);
   hemudist->SetLineColor(linecolor); 
   hsumptextra->SetLineColor(linecolor);
   hsumpzextra->SetLineColor(linecolor); 
   hdzextra->SetLineColor(linecolor); 
   hdxyextra->SetLineColor(linecolor); 
   hchi2extra->SetLineColor(linecolor); 
   hhitsextra->SetLineColor(linecolor); 
   hpurityextra->SetFillColor(fillcolor); 
   hmupt->SetLineColor(linecolor); 
   hmueta->SetLineColor(linecolor);  
   hept->SetLineColor(linecolor);  
   heeta->SetLineColor(linecolor);  
   hleadpt->SetLineColor(linecolor); 
   htrailpt->SetLineColor(linecolor); 
   hchi2pv->SetLineColor(linecolor);

   hmll->SetLineWidth(linewidth);   
   hdphi->SetLineWidth(linewidth);   
   hdpt->SetLineWidth(linewidth);   
   hpt->SetLineWidth(linewidth);   
   hmet->SetLineWidth(linewidth);   
   hntrk->SetLineWidth(linewidth);   
   hnvrt->SetLineWidth(linewidth);    
   hptextra->SetLineWidth(linewidth);   
   hetaextra->SetLineWidth(linewidth);   
   hvtxz->SetLineWidth(linewidth); 
   hvtxT->SetLineWidth(linewidth); 
   hmuemetdphi->SetLineWidth(linewidth); 
   hemudist->SetLineWidth(linewidth);
   hsumptextra->SetLineWidth(linewidth);
   hsumpzextra->SetLineWidth(linewidth); 
   hdzextra->SetLineWidth(linewidth);
   hdxyextra->SetLineWidth(linewidth); 
   hchi2extra->SetLineColor(linecolor);   
   hhitsextra->SetLineColor(linecolor);   
   hpurityextra->SetLineColor(linecolor);    
   hmupt->SetLineWidth(linewidth); 
   hmueta->SetLineWidth(linewidth);  
   hept->SetLineWidth(linewidth);  
   heeta->SetLineWidth(linewidth);  
   hleadpt->SetLineWidth(linewidth); 
   htrailpt->SetLineWidth(linewidth); 
   hchi2pv->SetLineWidth(linewidth);

   hmll->SetLineStyle(linestyle);    
   hdphi->SetLineStyle(linestyle);    
   hdpt->SetLineStyle(linestyle);    
   hpt->SetLineStyle(linestyle);    
   hmet->SetLineStyle(linestyle);    
   hntrk->SetLineStyle(linestyle);    
   hnvrt->SetLineStyle(linestyle);     
   hptextra->SetLineStyle(linestyle);    
   hetaextra->SetLineStyle(linestyle);    
   hvtxz->SetLineStyle(linestyle);  
   hvtxT->SetLineStyle(linestyle);  
   hmuemetdphi->SetLineStyle(linestyle);  
   hemudist->SetLineStyle(linestyle); 
   hsumptextra->SetLineStyle(linestyle); 
   hsumpzextra->SetLineStyle(linestyle);  
   hdzextra->SetLineStyle(linestyle);   
   hdxyextra->SetLineStyle(linestyle);    
   hchi2extra->SetLineStyle(linestyle);   
   hhitsextra->SetLineStyle(linestyle);   
   hpurityextra->SetLineStyle(linestyle);    
   hmupt->SetLineStyle(linestyle);  
   hmueta->SetLineStyle(linestyle);   
   hept->SetLineStyle(linestyle);   
   heeta->SetLineStyle(linestyle);   
   hleadpt->SetLineStyle(linestyle);  
   htrailpt->SetLineStyle(linestyle);  
   hchi2pv->SetLineStyle(linestyle);

   hmll->Sumw2();   
   hdphi->Sumw2();   
   hdpt->Sumw2();   
   hpt->Sumw2();   
   hmet->Sumw2();   
   hntrk->Sumw2();   
   hnvrt->Sumw2();    
   hptextra->Sumw2();   
   hetaextra->Sumw2();   
   hvtxz->Sumw2();
   hvtxT->Sumw2();
   hmuemetdphi->Sumw2();
   hemudist->Sumw2();
   hsumptextra->Sumw2();
   hsumpzextra->Sumw2(); 
   hdzextra->Sumw2();  
   hdxyextra->Sumw2();   
   hchi2extra->Sumw2(); 
   hhitsextra->Sumw2(); 
   hpurityextra->Sumw2(); 
   hmupt->Sumw2(); 
   hmueta->Sumw2();  
   hept->Sumw2();  
   heeta->Sumw2();  
   hleadpt->Sumw2();
   htrailpt->Sumw2();
   hchi2pv->Sumw2();

   hmll->Scale(xsec);    
   hdphi->Scale(xsec);    
   hdpt->Scale(xsec);    
   hpt->Scale(xsec);    
   hmet->Scale(xsec);    
   hntrk->Scale(xsec);    
   hnvrt->Scale(xsec);     
   hptextra->Scale(xsec);    
   hetaextra->Scale(xsec);    
   hvtxz->Scale(xsec);
   hvtxT->Scale(xsec);
   hmuemetdphi->Scale(xsec);
   hemudist->Scale(xsec);
   hsumptextra->Scale(xsec);
   hsumpzextra->Scale(xsec); 
   hdzextra->Scale(xsec);  
   hdxyextra->Scale(xsec);   
   hchi2extra->Scale(xsec);
   hhitsextra->Scale(xsec);
   hpurityextra->Scale(xsec);
   hmupt->Scale(xsec);  
   hmueta->Scale(xsec);   
   hept->Scale(xsec);   
   heeta->Scale(xsec);   
   hleadpt->Scale(xsec);
   htrailpt->Scale(xsec);
   hchi2pv->Scale(xsec);

  if(plotvar == 1)
    {
      hntrk->SetXTitle("N(extra tracks, e#mu vertex)"); 
      hntrk->SetYTitle("Events");
      return hntrk;  
    }
  if(plotvar == 2) 
    {
      hmll->SetXTitle("m(e#mu) [GeV]"); 
      if(thecuts < 777)
        hmll->SetYTitle("Events/6 GeV"); 
      if(thecuts == 3)
	hmll->SetYTitle("Events/12 GeV");
      if(thecuts == 777)
        hmll->SetYTitle("Events/30 GeV"); 
      if(thecuts > 777)
        hmll->SetYTitle("Events/100 GeV");  
       return hmll; 
    }
  if(plotvar == 3) 
    {
      hdphi->SetXTitle("1 - |#Delta #phi(e#mu)/#pi|"); 
      if(thecuts < 777) 
        hdphi->SetYTitle("Events/0.02");  
      if(thecuts == 3) 
        hdphi->SetYTitle("Events/0.1"); 
      if(thecuts == 777) 
        hdphi->SetYTitle("Events/0.05");  
      if(thecuts > 777) 
        hdphi->SetYTitle("Events/0.2");   
     return hdphi; 
    }
  if(plotvar == 4) 
    {
      hdpt->SetXTitle("#Delta p_{T}(e#mu) [GeV]");  
      if(thecuts < 777) 
	hdpt->SetYTitle("Events/2 GeV");
      if(thecuts == 777)
	hdpt->SetYTitle("Events/5 GeV");
      if(thecuts > 777) 
        hdpt->SetYTitle("Events/10 GeV"); 
      return hdpt; 
    }
  if(plotvar == 5) 
    {
      hpt->SetXTitle("p_{T}(e#mu) [GeV]");  
      if(thecuts < 777)
	hpt->SetYTitle("Events/6 GeV");
      if(thecuts >= 777)
        hpt->SetYTitle("Events/20 GeV"); 
      return hpt; 
    }
  if(plotvar == 6) 
    {
      hmet->SetXTitle("PF MET [GeV]");  
      if(thecuts < 777)
	hmet->SetYTitle("Events/6 GeV");
      if(thecuts == 777)
	hmet->SetYTitle("Events/60 GeV");
      if(thecuts > 777)
	hmet->SetYTitle("Events/30 GeV");
      return hmet; 
    }
  if(plotvar == 7) 
    {
      hnvrt->SetXTitle("N(primary vertices)");  
      hnvrt->SetYTitle("Events");
      return hnvrt; 
    }
  if(plotvar == 8)
    {
      hvtxz->SetXTitle("z position of e#mu vertex [cm]");   
      hvtxz->SetYTitle("Events/2 cm"); 
      return hvtxz;
    }
  if(plotvar == 9)
    {
      hvtxT->SetXTitle("Transverse position of e#mu vertex [cm]");    
      return hvtxT;
    }
  if(plotvar == 10)
    {
      hmuemetdphi->SetXTitle("1 - |#Delta #phi(e#mu - MET)/#pi|");    
      return hmuemetdphi;
    }
  if(plotvar == 11)
    {
      hemudist->SetXTitle("#Delta z(e#mu) [cm]");
      return hemudist;
    }
  if(plotvar == 12) 
    { 
      hmupt->SetXTitle("p_{T} (#mu) [GeV]"); 
      if(thecuts == 6)
	hmupt->SetYTitle("Events/5 GeV"); 
      if((thecuts < 777) && (thecuts != 6))
        hmupt->SetYTitle("Events/10 GeV");  
      if(thecuts >= 777)
        hmupt->SetYTitle("Events/50 GeV");   
      return hmupt; 
    } 
  if(plotvar == 13)  
    {  
      hept->SetXTitle("E_{T} (e) [GeV]");  
      if(thecuts == 6) 
        hept->SetYTitle("Events/5 GeV");  
      if((thecuts < 777) && (thecuts != 6)) 
        hept->SetYTitle("Events/10 GeV");   
      if(thecuts >= 777) 
        hept->SetYTitle("Events/50 GeV");    
      return hept;  
    }  
  if(plotvar == 14)  
    {  
      hmueta->SetXTitle("#eta (#mu)");  
      if(thecuts == 6) 
        hmueta->SetYTitle("Events/0.2");  
      if((thecuts < 777) && (thecuts != 6)) 
        hmueta->SetYTitle("Events/0.4");   
      if(thecuts >= 777) 
        hmueta->SetYTitle("Events/1.0");    
      return hmueta;  
    }  
  if(plotvar == 15)   
    {   
      heeta->SetXTitle("#eta (e)");   
      if(thecuts == 6)  
        heeta->SetYTitle("Events/0.2");   
      if((thecuts < 777) && (thecuts != 6))  
        heeta->SetYTitle("Events/0.4");    
      if(thecuts >= 777)  
        heeta->SetYTitle("Events/1.0");     
      return heeta;   
    }   
  if(plotvar == 16) 
    {
      hptextra->SetXTitle("Extra tracks p_{T} [GeV]");     
      hptextra->SetYTitle("Tracks/0.4 GeV");
      return hptextra; 
    }
  if(plotvar == 17) 
    {
      hetaextra->SetXTitle("Extra tracks #eta");      
      hetaextra->SetYTitle("Tracks/0.12"); 
      return hetaextra; 
    }
  if(plotvar == 18)
    {
      hsumptextra->SetXTitle("#Sigma p_{T} (extra tracks) [GeV]");
      hsumptextra->SetYTitle("Events/4 GeV");
      return hsumptextra;
    }
  if(plotvar == 19) 
    { 
      hsumpzextra->SetXTitle("#Sigma p - p_{Z} (extra tracks) [GeV]"); 
      return hsumpzextra; 
    } 
  if(plotvar == 20)
    {
      hleadpt->SetXTitle("p_{T} (leading lepton) [GeV]");
      return hleadpt;
    }
  if(plotvar == 21) 
    { 
      htrailpt->SetXTitle("p_{T} (trailing lepton) [GeV]"); 
      return htrailpt; 
    } 
  if(plotvar == 22)  
    {  
      hptntrackcorr->SetXTitle("N(extra tracks, e#mu vertex)"); 
      hptntrackcorr->SetYTitle("p_{T}(e#mu) [GeV]");   
      //      hptntrackcorr->Draw("col2z");
      return hptntrackcorr;
    }  
  if(plotvar == 24)
    {
      hdzextra->SetXTitle("#Delta z (extra tracks) [mm]");
      return hdzextra;
    }
  if(plotvar == 25)
    {
      hdxyextra->SetXTitle("#Delta xy (extra tracks) [mm]"); 
      return hdxyextra; 
    }
  if(plotvar == 26)   
    {   
      hchi2extra->SetXTitle("#chi^{2}/ndof (extra tracks)");    
      return hchi2extra;    
    }   
  if(plotvar == 27)   
    {   
      hhitsextra->SetXTitle("Hits (extra tracks)");    
      return hhitsextra;    
    }   
  if(plotvar == 28)   
    {   
      hpurityextra->SetXTitle("High purity (extra tracks)");    
      return hpurityextra;    
    }   
  if(plotvar == 29)
    {
      hchi2pv->SetXTitle("#chi^{2}/ndof (primary vertex)");
      return hchi2pv;
    }
}


