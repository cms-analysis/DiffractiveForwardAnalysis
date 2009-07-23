void makeEventsFile()
{
  TFile *f1 = new TFile("/tmp/jjhollar/lpair-mumu-10tev.root");
  TTree *t1 = f1->Get("h4444");
  Float_t p[4],px[4],py[4],pz[4],pt[4],eta[4],phi[4],en[4],m[4];
  Int_t partid[4];
  Float_t iz[4];

  t1->SetBranchAddress("px",px);
  t1->SetBranchAddress("py",py);
  t1->SetBranchAddress("pz",pz);
  t1->SetBranchAddress("E",en);
  t1->SetBranchAddress("m",m);
  t1->SetBranchAddress("iz",iz);
  t1->SetBranchAddress("icode",partid);

  ofstream output("/tmp/jjhollar/gamgammumu.lpairelastic.10tev.lhe");

  Int_t nevts = t1->GetEntries();

  output << "<LesHouchesEvents version=\"1.0\">"  << endl;
  output << "<header>" << endl;
  output << "This file was created from the output of the LPAIR generator" << endl;
  output << "</header>" << endl;

  output << "<init>" << endl;
  output << "2212  2212  0.50000000000E+04  0.50000000000E+04 0 0 10042 10042 2  1" << endl;
  output << "0.10508723460E+01  0.96530000000E-02  0.26731120000E-03   0" << endl;
  output << "</init>" << endl;

  for(Int_t i = 0;i < 1000;i++)
    {
      t1->GetEntry(i);
      output << "<event>" << endl;
      output << "6   0  0.2983460E-04  0.9118800E+02  0.7546772E-02  0.1300000E+00" << endl;
      // JH - note here we add in two fake photons as the beam particles. The energies don't matter - this is only for 
      // the LHE event record. 
      output << "22   -1    0    0    0    0  0.00000000000E+00  0.00000000000E+00 0.00000000000E+02  0.10000000000E+02  0.00000000000E+00 0.  1." << endl;
      output << "22   -1    0    0    0    0  0.00000000000E+00  0.00000000000E+00 0.00000000000E+00  0.10000000000E+02  0.00000000000E+00 0. -1." << endl; 
      output << partid[0] << " 1 1 2 0 0 " << px[0] << " " << py[0] << " " << pz[0] << " " << en[0] << " " << m[0] << " 0. " << iz[0] << endl; 
      output << partid[1] << " 1 1 2 0 0 " << px[1] << " " << py[1] << " " << pz[1] << " " << en[1] << " " << m[1] << " 0. " << iz[1] << endl;  
      output << partid[2] << " 1 1 2 0 0 " << px[2] << " " << py[2] << " " << pz[2] << " " << en[2] << " " << m[2] << " 0. " << iz[2] << endl;  
      output << partid[3] << " 1 1 2 0 0 " << px[3] << " " << py[3] << " " << pz[3] << " " << en[3] << " " << m[3] << " 0. " << iz[3] << endl;  
      output << "</event>" << endl;
    }
  output << "</LesHouchesEvents>" << endl;

}
