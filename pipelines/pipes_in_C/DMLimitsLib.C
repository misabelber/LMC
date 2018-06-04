#include "/home/queenmab/GitHub/Math/matrixes.h"
#include "/home/queenmab/GitHub/Math/linear_fitter.h"
#include "TMatrixD.h"
#include "TFITS.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TMultiGraph.h"
#include "TF2.h"
#include "TRandom.h"
#include "TFile.h"
#include <ctime>


#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

void ReadFits(){
  
  cout << "USAGE: " << endl;
  cout << endl;
  cout << "ReadFits(VM &data, TString filename)" << endl;
  cout << endl;
  cout << "data: Empty Vector of Matrixes to store fits file binned image" << endl;
  cout << "filename: Fits file name" << endl;
  cout << endl;
}
