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

void contours()
{
// we import the file (ppbarContours.root) with the contours for p+pbar produced in Pb-Pb collisions
	TFile *ppbarFile = new TFile("ppbarContours.root");
// we import the file (phiContours.root) with the contours for phi mesons produced in Pb-Pb collisions
	TFile *phiFile = new TFile("phiContours.root");

// we define a canvas
	TCanvas *T1 = new TCanvas("T1","Contours",200,200);

// il manque des options gstyle ou gpad

// we define the graphs for each contour 
	TGraph *ppbar_1sigma[10];
	TGraph *ppbar_2sigma[10];
	TGraph *phi_1sigma[10];
	TGraph *phi_2sigma[10];

// we define the multigraph that will contain all the contours
	TMultiGraph *mg = new TMultiGraph();
	mg->GetXaxis()->SetTitle("#LT #beta_{T} #GT");
	mg->GetYaxis()->SetTitle("T_{kin}");

//------------------------------------------------------
//----------------- p+pbar -----------------------------	
//------------------------------------------------------
// we explorate the ppbar file in order to get the 9 first pairs of contours it contains
	char ppbar_1sigma_name = "ppbar_1sigma_contour_X";
	char ppbar_2sigma_name = "ppbar_2sigma_contour_X";
	for(int digit=1;digit<=9;digit++)
	{
	// we do it for the 1sigma contours
		ppbar_1sigma_name[21] = digit+'0';
		ppbar_1sigma[digit] = (TGraph*)ppbarFile->Get("ppbar_1sigma_name") ; 
		ppbar_1sigma[digit]->Draw();   // the "A" option has to be there ! 
		mg->Add(ppbar_1sigma[digit]);
	// we do it for the 2sigma contours
		ppbar_2sigma_name[21] = digit+'0';
		ppbar_2sigma[digit] = (TGraph*)ppbarFile->Get("ppbar_2sigma_name") ; 
		ppbar_2sigma[digit]->Draw();   // the "A" option has to be there ! 
		mg->Add(ppbar_2sigma[digit]);
	}
	
	
// we add the 10th pair of contours	
	ppbar_1sigma[10] = (TGraph*)ppbarFile->Get("ppbar_1sigma_contour_10") ; 
	ppbar_1sigma[10]->Draw();   // the "A" option has to be there ! 
	mg->Add(ppbar_1sigma[10]);	
	ppbar_2sigma[10] = (TGraph*)ppbarFile->Get("ppbar_2sigma_contour_10") ; 
	ppbar_1sigma[10]->Draw();   // the "A" option has to be there ! 
	mg->Add(ppbar_2sigma[10]);
	
	

//------------------------------------------------------
//----------------- phi mesons -----------------------------	
//------------------------------------------------------
// we explorate the ppbar file in order to get the 9 first pairs of contours it contains
	char phi_1sigma_name = "phi_1sigma_contour_X";
	char phi_2sigma_name = "phi_2sigma_contour_X";
	for(int digit=1;digit<=10;digit++)
	{
	// we do it for the 1sigma contours
		phi_1sigma_name[21] = digit+'0';
		phi_1sigma[digit] = (TGraph*)phiFile->Get("phi_1sigma_name") ; 
		phi_1sigma[digit]->Draw();   // the "A" option has to be there ! 
		mg->Add(phi_1sigma[digit]);
	// we do it for the 2sigma contours
		phi_2sigma_name[21] = digit+'0';
		phi_2sigma[digit] = (TGraph*)phiFile->Get("phi_2sigma_name") ; 
		phi_2sigma[digit]->Draw();   // the "A" option has to be there ?! 
		mg->Add(phi_2sigma[digit]);
	}
	
	
// we add the 10th pair of contours	
	phi_1sigma[10] = (TGraph*)phiFile->Get("phi_1sigma_contour_10") ; 
	phi_1sigma[10]->Draw();   // the "A" option has to be there ! 
	mg->Add(phi_1sigma[10]);	
	phi_2sigma[10] = (TGraph*)phiFile->Get("phi_2sigma_contour_10") ; 
	phi_1sigma[10]->Draw();   // the "A" option has to be there ?! 
	mg->Add(phi_2sigma[10]);

//------------------------------------------------------
//------------------------------------------------------	
// we draw the multigraph
	mg->Draw("A");	  // the "A" option has to be there ?! 
// change the axis limits
	gpad->Modified();   
	mg->GetXaxis()->SetLimits(0.3,0.9); // we limit the x axis range
	// then we limit the y axis range
	mg->SetMinimum(0.05);
	mg->SetMinimum(0.2);

	
// we add a legend
	TLatex Tl;
	Tl.SetNDC();
	Tl.SetTextSize(0.020);
	Tl.SetTextAlign(13);  //align at top
	Tl.SetTextFont(42);
// we indicate the centrality
	latex.DrawLatex(0.2,0.8, "80-90%");  // top left of the plot
	latex.DrawLatex(0.8,0.2, "0-5%");    // bottom right of the plot
// we indicate the color code
	latex.DrawLatex(0.5,0.65, "Blast-wave fit to");
	latex.DrawLatex(0.5,0.5, "#phi (?-? GeV/c), p+#bar{p} (0.3-3 GeV/c)");
	
// we close the files
	ppbarFile->Close();
	phiFile->Close();
	
}
