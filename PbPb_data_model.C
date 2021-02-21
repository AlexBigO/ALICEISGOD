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

// we import the file (ppbarData_Model.root) containing the data/model graphs for p+pbar produced in Pb-Pb collisions
	TFile *ppbarFile = new TFile("ppbar_PbPb_Data_Model.root");
// we import the file (phiData_Model.root) containing the data/model graphs for phi produced in Pb-Pb collisions
	TFile *phiFile = new TFile("phi_PbPb_Data_Model.root");
	
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
	TGraph *ppbar[11];  
	TGraph *phi[9];

// we define the multigraph that will contain the values for p+pbar
	TMultiGraph *mg_ppbar = new TMultiGraph();
	mg_ppbar->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	mg_ppbar->GetYaxis()->SetTitle("Data/BW fit");
// we define the multigraph that will contain the values for phi
	TMultiGraph *mg_phi = new TMultiGraph();
	mg_phi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	mg_phi->GetYaxis()->SetTitle("Data/BW fit");

//------------------------------------------------------
//----------------- p+pbar -----------------------------	
//------------------------------------------------------
// we explorate the ppbar file in order to get the 9 first graphs it contains
	char ppbar_name[] = "ppbar_data_model_X";
	for(int digit=1;digit<=9;digit++)
	{
		ppbar_name[17] = digit+'0';
		ppbar[digit] = (TGraph*)ppbarFile->Get(ppbar_name) ; 
		ppbar[digit]->Draw();   
		mg_ppbar->Add(ppbar[digit]);
	}	
// we add the 10th graph	
	ppbar[10] = (TGraph*)ppbarFile->Get(ppbar_data_model_10) ; 
	ppbar[10]->Draw();  
	mg_ppbar->Add(ppbar[10]);	

	
	

//------------------------------------------------------
//----------------- phi --------------------------------	
//------------------------------------------------------
// we explorate the ppbar file in order to get the 8 graphs it contains
	char phi_name[] = "phi_data_model_X";
	for(int digit=1;digit<=8;digit++)
	{
		phi_name[15] = digit+'0';
		phi[digit] = (TGraph*)ppbarFile->Get(phi_name) ; 
		phi[digit]->Draw();   
		mg_phi->Add(phi[digit]);
	}	





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
 	mg_phi->Draw("A");

//------------------------------------------------------------
//------------------- BOTTOM GRAPH ------------------------------
//------------------------------------------------------------   	
   	// we create the graph to plot Data/Model fit with rectangle errors for p+pbar
	DATA_MODEL->cd(2); 
	gPad->SetTopMargin(small);	
	mg_ppbar->Draw("A");
	
	
	
	
	
	
	
// we close the files
	ppbarFile->Close();
	phiFile->Close();	

}

