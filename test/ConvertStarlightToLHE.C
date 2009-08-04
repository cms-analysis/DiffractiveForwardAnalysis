/* Routine for the conversion of starlight data format
 * into something readable by CMSSW_1_4_5
 *
 * Modification by X. Rouby on a routine kindly offered by J. Hollar
 * Sept. 28, 2007
 */ 

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
using namespace std;

void makeEventsFile(int num=0)
{
  string filename = "/tmp/jjhollar/starlight_u3s_pp_10TeV_mumu.out"; //input
  ifstream infile(filename.c_str());
  char outfilename[100];
  sprintf(outfilename,"/tmp/jjhollar/starlight10tev_u3s_mumu.lhe",num);
  ofstream output(outfilename);
  if (! infile.is_open()) { cout << "\t ERROR: I can not open \"" << filename << "\"" << endl; return; }

  string temp_string, temp;
  istringstream curstring;
  const unsigned int N = 2; // N_particles
  const unsigned int M = 10000; // N_events
  const unsigned int K = num*M+1;  // first event 
  const double MU = 0.105658369; // muon mass [GeV]
  double charge = 0.0;
  int evt_n=0; // event_number, read from the input-file
  int nn=0; // event_counter, in the output file

  // GENER: STL  1.0            1           1           1           1  200.00     999.999     CMS
  getline(infile,temp_string); // The very first line is useless
  output << "<LesHouchesEvents version=\"1.0\">"  << endl;
  output << "<header>" << endl;
  output << "This file was created from the output of the STARLIGHT generator" << endl;
  output << "</header>" << endl;

  output << "<init>" << endl;
  output << "2212  2212  0.50000000000E+04  0.50000000000E+04 0 0 10042 10042 2  1" << endl;
  output << "0.10508723460E+01  0.96530000000E-02  0.26731120000E-03   0" << endl;
  output << "</init>" << endl;

  while (getline(infile,temp_string)) {
    
    curstring.clear(); // needed when using several tims istringstream::str(string)
    curstring.str(temp_string);

    if(strstr(temp_string.c_str(),"EVENT")){
      curstring >> temp >> evt_n;
      // EVENT:          1       2       1
      if(evt_n >=K && evt_n < K+M) {
	if(nn > 0)
	  output << "</event>" << endl;
	output << "<event>" << endl;
	output << "4   0  0.2983460E-04  0.9118800E+02  0.7546772E-02  0.1300000E+00" << endl;
	// JH - note here we add in two fake photons as the beam particles. The energies don't matter - this is only for   
	// the LHE event record.   
	output << "22   -1    0    0    0    0  0.00000000000E+00  0.00000000000E+00 0.00000000000E+02  0.10000000000E+02  0.00000000000E+00 0.  1." << endl;  
	output << "22   -1    0    0    0    0  0.00000000000E+00  0.00000000000E+00 0.00000000000E+00  0.10000000000E+02  0.00000000000E+00 0. -1." << endl;   

	nn++;
      }

    } else if(strstr(temp_string.c_str(),"VERTEX")){
      float x,y,z,t;
      curstring >> temp >> x >> y >> z >> t; 
      // VERTEX:   0.0000       0.0000       0.0000       0.0000          1      0      0      2
      //      if(evt_n >=K && evt_n < K+M) output << "V " << evt_n << " 0 " << x << " " << y << " " << z << " " << t << " 0 " << N <<  " 0" << endl;

    } else if(strstr(temp_string.c_str(),"TRACK")){
      int useless, part_n, pdg_id;
      float px, py, pz;
      //TRACK:      6   2.9797       3.1399       84.461          1      1      0    -13
      curstring >> temp >> useless >> px >> py >> pz >> part_n >> useless >> useless >> pdg_id;
      //P 5 2212 -2.0 1.2 8.1 5.2 1 0 0 0 0   
      //      if(evt_n >=K && evt_n < K+M) output << "P " << part_n+(evt_n-1)*N << " " << pdg_id << " " << px << " " << py << " " << pz << " " << sqrt(MU*MU + px*px + py*py + pz*pz) << " 1 0 0 0 0\n";
      if(pdg_id == 13)
	charge = -1.;
      if(pdg_id == -13)
	charge = 1.;

      if(evt_n >=K && evt_n < K+M) 
	output << pdg_id << " 1 1 2 0 0 " << px << " " << py << " " << pz << " " << sqrt(MU*MU + px*px + py*py + pz*pz) << " " << MU << " 0. " << charge << endl;
    }

  } // reading loop of the input file
  output << "</LesHouchesEvents>" << endl;
  infile.close();
  output.close();
  cout << nn << " events written in " << outfilename << endl;
  return;
}

void makefiles(int number_of_files=5) {
  for(int i=0; i<number_of_files; i++)
    makeEventsFile(i);
}

/* Explaination of the format :
 * +++ Event +++
 * E 1 -1.0000000000000000e+00 -1.0000000000000000e+00 -1.0000000000000000e+00 20 0 1 0 0
 *   1 : event number  <-------
 *   -1 : event scale
 *   -1 : alpha_QCD
 *   -1 : alpha_QED
 *   20 : signal process ID
 *   0 : signal process vertex barcode
 *   1 : number of vertices <-------
 *   0 : list of vertices
 *   0 : ?
 *
 * +++ Vertex +++
 * V -1 0 0 0 0 0 0 4 0
 *   -1 : vertex barcode (unique)       <-------
 *    0 : vertex id
 *    0 0 0 0 : vertex x,y,z,t
 *    0 : number of orphans
 *    4 : number of out particles       <-------
 *    0 : weights
 *
 * +++ Particle +++
 * P 5 2212 -2.0 1.2 8.1 5.2 1 0 0 0 0   
 *    5 : barcode<-------
 *    0 : pdg_id<-------
 *   -2.0 : px<-------
 *    1.2 : py<-------
 *    8.1 : pz<-------
 *    5.2 : e<-------
 *    1 : status<-------
 *    0 0  : polarization eta , phi
 *    0 0  : vertex and ?
 */
