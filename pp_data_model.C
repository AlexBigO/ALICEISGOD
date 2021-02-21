#include <string.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <stdbool.h>
#include <iostream>
#include <TColor.h>
#include <TPaveStats.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TArrow.h>
#include <TMath.h>
#include <TFile.h> 
#include <TDirectoryFile.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <fstream>

#if defined(__GLIBCXX__)
#pragma link C++ class MyOtherClass;
#endif


using namespace std;

void data_model()
{

// we import the file (ppbarData_Model.root) containing the data/model graphs for p+pbar produced in pp collisions
	TFile *ppbarFile = new TFile("ppbar_pp_Data_Model.root");
// we import the file (phiData_Model.root) containing the data/model graphs for phi produced in pp collisions
	TFile *phiFile = new TFile("phi_pp_Data_Model.root");
	
// we define a canvas
	TCanvas *T1 = new TCanvas("DATA_MODEL","DATA_MODEL",200,10,600,400);
	T1->SetGrid();
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetOptStat(0);
	float small = 1e-5;
	// we divide the canvas in 2 sub-graphs: 1 at the top and 1 at the bottom 
	T1->Divide(1,2,small,small);

// we define a TGraph for each data/model graph 
	TGraph *ppbar[1];  
	TGraph *phi[1];



//------------------------------------------------------
//----------------- p+pbar -----------------------------	
//------------------------------------------------------
// we explorate the ppbar file
	ppbar[0] = (TGraph*)ppbarFile->Get("ppbar_data_model_1") ; 
	ppbar[0]->GetYaxis()->SetTitle("Data/BW fit");

//------------------------------------------------------
//----------------- phi --------------------------------	
//------------------------------------------------------
// we explorate the ppbar file in order to get the 8 graphs it contains
	phi[0] = (TGraph*)phiFile->Get("phi_data_model_1") 	
	phi[0]->GetYaxis()->SetTitle("Data/BW fit");
	phi[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");



//------------------------------------------------------------
//------------------------------------------------------------
//------------------- PLOT OF THE MULTIGRAPH -----------------
//------------------------------------------------------------
//------------------------------------------------------------


//------------------------------------------------------------
//------------------- TOP GRAPH ------------------------------
//------------------------------------------------------------
	// we create the graph to plot Data/Model fit with rectangle errors for phi
	DATA_MODEL->cd(1);
	gPad->SetTickx(2);
	gPad->SetBottomMargin(small); 
	gPad->SetTickx();
        //gPad->SetTicky(1);
   	//mg_phi->GetYaxis()->SetLabelOffset(0.01);
 	ppbar[0]->Draw("A");

//------------------------------------------------------------
//------------------- BOTTOM GRAPH ------------------------------
//------------------------------------------------------------   	
   	// we create the graph to plot Data/Model fit with rectangle errors for p+pbar
	DATA_MODEL->cd(2); 
	gPad->SetTopMargin(small);	
	phi[0]->Draw("A");
	
	
	
	
	
	
	
// we close the files
	ppbarFile->Close();
	phiFile->Close();	

}

