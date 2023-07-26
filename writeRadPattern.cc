#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <math.h>
#include <iterator>
#include <array>   
#include <cmath> 
#include <algorithm>  

using namespace std;

// Declare functions
vector< vector<double> > readGenes(string inFilePath, int gen, int npop, bool isPhase);
double calcGain(vector<double> genes, int theta, int phi, int freq, bool isDB);

const float PI = 3.14159265359;

int main(int argc, char const *argv[]){

  // Output usage message if wrong args
  if (argc < 3){
    std::cerr << "Usage: " << argv[0] << " <GEN> <NPOP> <SOURCE_DIR>" << endl;;
    return 1;
  }
  
  // Parse args
  int gen = stoi(argv[1]);
  int npop = stoi(argv[2]);
  string source_dir = argv[3];

  // Declare vectors to hold genes from input files
  vector< vector<double> > inputGenesGain(npop, vector<double>(14, 0));
  vector< vector<double> > inputGenesPhase(npop, vector<double>(14, 0));

  // Read gain coefficients from input file
  inputGenesGain = readGenes(source_dir, gen, npop, false);
  inputGenesPhase = readGenes(source_dir, gen, npop, true);

  // Declare and populate theta and phi arrays
  double theta[37];
  double phi[72];
  int freq[60];

  for (int th = 0; th < 37; th++){
    theta[th] = (double) th*5.0f;
  }

  for (int p = 0; p < 72; p++){
    phi[p] = (double) p*5.0f;
  }

  // Loop over all individuals, calculate and write gain, gainDB, phase
  for (int indiv=0; indiv<npop; indiv++){

    // Output status message:
    cout << "Writing rad pattern for indiv " + std::to_string(indiv) + "\n";

    // Set up input and output files for each individual
    string full_file_path = source_dir + "/gen_" + std::to_string(gen) + "/child_" + std::to_string(indiv) + ".txt";
    fstream inFile(full_file_path.c_str());
    stringstream output;
    //ofstream outFile;

    // Declare and populate gain and phase arrays to hold rad pattern values
    //double gain[60][72][37] = {{{0}}};
    //double gainDB[60][72][37] = {{{0}}};
    //double phase[60][72][37] = {{{0}}};

    // Define variables to hold temp values to be overwritten in each iteration
    double tempGain = 0.0;
    double tempGainDB = 0.0;
    double tempPhase = 0.0;

    // Loop over frequency steps
    for (int f=0; f<60; f++){
      // Write "header" for each frequency step
      output << "freq : " << 83.33 + 16.67*f << " MHz" << endl;
      output << " MHz\rSWR : 0.0\n"; 
      output << " Theta \t Phi \t Gain(dB)   \t   Gain   \t   Phase(deg)\n";

      // Loop over phi and theta for each frequency and calc gains and phases
      for (int p = 0; p < 72; p++){
	for (int th = 0; th < 37; th++){

	  // Calculate values at each theta, phi, freq step
	  tempGain = calcGain(inputGenesGain[indiv], theta[th], phi[p], f, false);
	  tempGainDB = calcGain(inputGenesGain[indiv], theta[th], phi[p], f, true);
	  tempPhase = calcGain(inputGenesPhase[indiv], theta[th], phi[p], f, false);

	  // Write values to stringstream output
	  output  << theta[th] << " \t " << phi[p] << " \t " << tempGainDB << "     \t   " << tempGain << "     \t    " << tempPhase << endl;

	}
      }
    }
    
    // Add empty line to match previous formatting
    output << "\n";

    // Add fitness info that was read from input and write everything
    output << inFile.rdbuf();
    inFile.close();

    inFile.open(full_file_path.c_str(), fstream::out | fstream::trunc);
    inFile << output.rdbuf();
    
  }

  return 0;
} 

vector<vector<double>> readGenes(string inFilePath, int gen, int npop, bool isPhase){

  int nColumns = 14;
  int nRows = 1;
  int start_line = 0;
  int line_num = 0;
  int column_count = 0;

  vector<vector<double>> genes(npop, vector<double>(nColumns, 0));

  string full_file_path;

  ifstream inFile;

  // HARD CODE FOR NOW
  // Lines to read for genes
  if (isPhase){
    start_line = 8;
  }
  else{
    start_line = 3;
  }

  for (int i=0; i<npop; i++){

    //Line number starts at 1
    line_num = 1;

    // Open input file
    full_file_path = inFilePath + "/gen_" + std::to_string(gen) + "/child_" + std::to_string(i) + ".txt";
    inFile.open(full_file_path);
    
    if(inFile.is_open()){
      
      // Read until end of file
      while(!(inFile.eof())){
	
	string line;
	getline(inFile, line);
	
	if (line_num == start_line){
	  
	  column_count = 0;
          istringstream iss(line);
          string element;

	  // Data is delimitted by ','
	  while(getline(iss, element, ',')){
	    
	    // Store genes in 2D vector
	    // one entry per individual, each w/ 14 genes
	    string::size_type sz;
	    genes[i][column_count] = stof(element, &sz);
	    
	    column_count++;
	  }

	}
	line_num++;
      }

      inFile.close();
    }

    // If file not open, output error message
    else{
      cout << "ERROR: Problem opening file at : " << full_file_path << endl;
    }
    
  }
  
  return genes;
}

// Function to calculate gain, gainDB, phase at given theta, phi, freq
double calcGain(vector<double> genes, int th, int p, int freq, bool isDB){

  // NOTE: The phase calculation is the same as for gain, so usage is the same

  // NOTE: the gain calculation is gross, but leaving it that way to match the previous version from the GA

  // Declare variables
  double gain = 0.0;
  double gainDB = 0.0;

  gain = /*(1+(genes[13]*(16.67*freq)/83.33))**/(genes[0]*(1/2.0)*(1/sqrt(PI)) +
                                     genes[1]*(1/2.0)*sqrt(3/PI)*cos(th*PI/180.0f) +
                                     genes[2]*(1/4.0)*sqrt(5/PI)*(3*pow(cos(th*PI/180.0f), 2)- 1) +
                                     genes[3]*(1/4.0)*sqrt(7/PI)*(5*pow(cos(th*PI/180.0f),3)- 3*cos(th*PI/180.0f)) +
                                     genes[4]*(3/16.0)*sqrt(1/PI)*(35*pow(cos(th*PI/180.0f),4) - 30*pow(cos(th*PI/180.0f),2)+3) +
                                     genes[5]*(1/16.0)*sqrt(11/PI)*(15*cos(th*PI/180.0f) - 70*pow(cos(th*PI/180.0f),3)+63*pow(cos(th*PI/180.0f),5)) +
                                     genes[6]*(1/32.0)*sqrt(13/PI)*(-5 + 105*pow(cos(th*PI/180.0f),2)-315*pow(cos(th*PI/180.0f),4) + 231*pow(cos(th*PI/180.0f),6)) +
                                     genes[7]*(1/32.0)*sqrt(15/PI)*(-35*cos(th*PI/180.0f)+ 315*pow(cos(th*PI/180.0f),3) -693*pow(cos(th*PI/180.0f),5) + 429*pow(cos(th*PI/180.0f),7)) +
                                     genes[8]*(1/256.0)*sqrt(17/PI)*(35 - 1260*pow(cos(th*PI/180.0f),2) + 6930*pow(cos(th*PI/180.0f),4) - 12012*pow(cos(th*PI/180.0f),6) + 6435*pow((cos(th*PI/180.0f)),8)) +
                                     genes[9]*(1/256.0)*sqrt(19/PI)*(315*cos(th*PI/180.0f)- 4620*pow(cos(th*PI/180.0f),3) + 18018*pow(cos(th*PI/180.0f),5) - 25740*pow(cos(th*PI/180.0f),7) + 12155*pow((cos(th*PI/180.0f)),9)) +
                                     genes[10]*(1/512.0)*sqrt(21/PI)*(-63 +3465*pow(cos(th*PI/180.0f),2) - 30030*pow(cos(th*PI/180.0f),4) + 90090*pow(cos(th*PI/180.0f),6) -109395*pow((cos(th*PI/180.0f)),8)+46189*pow(cos(th*PI/180.0f),10)) +  
                                     genes[11]*(1/512.0)*sqrt(23/PI)*(-693*pow(cos(th*PI/180.0f),1) +15015*pow(cos(th*PI/180.0f),3) - 90090*pow(cos(th*PI/180.0f),5) +218790*pow((cos(th*PI/180.0f)),7)-230945*pow(cos(th*PI/180.0f),9)+88179*pow(cos(th*PI/180.0f),11)) +
                                     genes[12]*(1/2048.0)*sqrt(25/PI)*(231 -18018*pow(cos(th*PI/180.0f),2) +225225*pow(cos(th*PI/180.0f),4) - 1021020*pow(cos(th*PI/180.0f),6) +2078505*pow((cos(th*PI/180.0f)),8)-1939938*pow(cos(th*PI/180.0f),10)+676039*pow(cos(th*PI/180.0f),12)));

  if (isDB){
    if (gain > 0){
      gainDB = 10*log10(gain);
    }
    else{
      gainDB = -80.000;
    }
    
    return gainDB;

  }

  return gain;

}

