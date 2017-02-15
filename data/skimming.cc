#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TVector2.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

Float_t dphi_0_phi(Float_t phi1, Float_t phi2){
  Float_t dphi_tmp = fabs(phi1-phi2);
  if (dphi_tmp > 3.141593) dphi_tmp=fabs(dphi_tmp-2*3.141593);
  return dphi_tmp;
}

int main(int argc, char** argv) {

  if( argc<6 ) {
     std::cout << argv[0] << ":  " << std::endl ;
     std::cout << " Functionality: skimming... "  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: " << argv[0] << " inputfile.root outputfile.root Nevts SumWeights cut_string " << std::endl ;
     exit(1) ;
  }

  // input file name
  std::string inputfile((const char*)argv[1]); 
  // output file name
  std::string outputfile((const char*)argv[2]);

  // nevts 
  double SumEvents = atof(argv[3]);  
  // sum weights 
  double SumWeights = atof(argv[4]);  
  // cuts
  std::string Cuts((const char*)argv[5]);

  // initialize
  // root files
  TFile* finput = new TFile(inputfile.c_str());
  TFile* foutput = new TFile(outputfile.c_str(), "recreate");

  // tree
  TTree* tree_in = (TTree*)finput->Get("tree");
  // tree_in->SetBranchStatus("met_*",0); // used to shut down some non-use branches
  std::cout<<"[info]: in tree entries = "<<tree_in->GetEntriesFast()<<std::endl;
  
  // out_tree
  TTree* tree_out = tree_in->CopyTree(Cuts.c_str());
  std::cout<<"[info]: out tree entries = "<<tree_out->GetEntriesFast()<<std::endl;
  
  TBranch *b_SumEvents=tree_out->Branch("SumEvents",&SumEvents,"SumEvents/D");
  TBranch *b_SumWeights=tree_out->Branch("SumWeights",&SumWeights,"SumWeights/D");

  for (int ii=0; ii<(int)tree_out->GetEntries(); ii++){

    b_SumEvents->Fill();
    b_SumWeights->Fill();    
    
    if (ii%1000==0) std::cout<<"Event "<<ii<<std::endl;
    //if (ii==200) break;
  }

  tree_out->Write();
  foutput->Close();
  finput->Close();

  return 0;

}



