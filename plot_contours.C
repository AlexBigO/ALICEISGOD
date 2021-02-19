#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <stdbool.h>
#include <vector>
#include <iostream>
#include <vector>
#include <math.h> 
#include <TRandom.h>
#include <TColor.h>
#include <TPaveStats.h>
#include <TList.h>
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
#include <TMinuit.h>


void plot_contours()
{
// importation of the file containing the p+pbar contours
	TFile *myFile1 = new TFile("contours_proton.root");
// importation of the file containing the phi mesons contours
	//TFile *myFile2 = new TFile("contours_phi.root");

// creation of the canvas that will hold the multigraph
	TCanvas *c1 = new TCanvas("C1" , " C1 " , 200 , 200 );

	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetOptStat(0);

// create a TGraph for each contour 
	TGraph *ppbar_1sigma[10];
	TGraph *ppbar_2sigma[10];
// create a multigraph 
	TMultiGraph  *mg  = new TMultiGraph();
   	mg->GetXaxis()->SetTitle("#LT #beta_{T} #GT");
	mg->GetYaxis()->SetTitle("T_{kin}");	

// we get the 9 first pairs of contours we saved in the TFile
	int digit;
	char sigma_1[] = "Err_def_1_digit_Y";
	char sigma_2[] = "Err_def_9_digit_Y";
	for(int digit=1;digit<=9;digit++)
	{
		sigma_1[16] = digit+'0';
		// we get the 1-sigma contours and add them to the multigraph
		ppbar_1sigma[digit] = (TGraph*)myFile1->Get(sigma_1);
		ppbar_1sigma[digit]->Draw();
		mg->Add(ppbar_1sigma[digit]);
		// we get the 2-sigma contours and add them to the multigraph
		sigma_2[16] = digit+'0';
		ppbar_2sigma[digit] = (TGraph*)myFile1->Get(sigma_2);
		ppbar_2sigma[digit]->Draw();
		mg->Add(ppbar_2sigma[digit]);
	}

/*	
// we get the 10th pair of contours	
	ppbar_1sigma[10] = (TGraph*)myFile1->Get("Err_def_1_digit_10");
	ppbar_1sigma[10]->Draw();
	mg->Add(ppbar_1sigma[10]); 
	ppbar_2sigma[10] = (TGraph*)myFile1->Get("Err_def_4_digit_10");
	ppbar_2sigma[10]->Draw();
	mg->Add(ppbar_1sigma[10]); 
*/

// we draw the multigraph
	mg->Draw("A");


	// Change the axis limits
  	gPad->Modified();
	mg->GetXaxis()->SetLimits(0.8,0.92);  // we limit the x axis range
	// then we limit the y axis range
	mg->SetMinimum(0.05);
	mg->SetMaximum(0.2);	
	
// Ajout des résultats des fits sur le canvas
	TLatex Tl;
		Tl.SetNDC();
		Tl.SetTextSize(0.040);
		Tl.SetTextAlign(13);  //align at top
		Tl.SetTextFont(50);
	   
	 
		Tl.DrawLatex(0.2,0.8, "80-90%" );
		Tl.DrawLatex(0.8,0.5, "0-5% " );
		
		Tl.DrawLatex(0.2,0.2, "Blast-wave fit to #beta" );
		Tl.DrawLatex(0.2,0.18, TString::Format("#phi (0 GeV/#c), p+#bar{p} (0.3-3.0 GeV/#c)"));	
		Tl.DrawText(.1, .5,   "#sqrt{10} #sqrt[3]{10} :"); Tl.DrawLatex(.5, .5, "#sqrt{10} #sqrt[3]{10}");

	myFile1 -> Close();
	//myFile2 -> Close();
}
