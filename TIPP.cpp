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
#include <TGraphErrors.h>
#include <TArrow.h>
#include <TMath.h>
#include <TFile.h> 
#include <TDirectoryFile.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraph.h>

#if defined(__GLIBCXX__)
#pragma link C++ class MyOtherClass;
#endif


using namespace std;

		
// Exponential law
// par[0] = C
// par[1] = m_0
// par[2] = T
double expo_law(double *x , double *par)
{
	double m_T = sqrt(x[0]*x[0] + par[1]*par[1]);
	double res = par[0]*TMath::Power( par[2]*(par[1]+par[2]), -1)*x[0]*exp(- (m_T - par[1]) / par[2]);
	return res;
}


// Boltzmann distribution
// par[0] = C
// par[1] = m_0
// par[2] = T
double boltzmann(double *x , double *par)
{
	double m_T = sqrt(x[0]*x[0] + par[1]*par[1]);
	double res = 2*TMath::Pi()*par[0]*x[0]*m_T*exp(- (m_T - par[1]) / par[2]);
	return res;
}


// Levy-Tsallis distribution
// par[0] = C
// par[1] = m
// par[2] = T
// par[3] = n
double levy(double *x , double *par)
{
	double m_T = sqrt(x[0]*x[0] + par[1]*par[1]);
	double res = par[0]*(par[3]-1)*(par[3]-2)/(par[3]*par[2]*(par[3]*par[2]+par[1]*(par[3]-2))) *x[0]*TMath::Power(1+(m_T - par[1]) / (par[3]*par[2]) , - par[3]);
	return res;
}


// Power law from [5] in biblio
// par[0] = C = dN / dy
// par[1] = p_0
// par[2] = n
double power_law_Five(double *x , double *par)
{
	double res = par[0]*4*(par[2]-1)*(par[2]-2)*TMath::Power(par[2]-3,-2)*TMath::Power(par[1],-2)*x[0]*TMath::Power(1+ x[0] / (par[1]* (par[2]-3)/2 ) , - par[2]);
	return res;
}


// "Power law" distribution  from [7] in biblio
// par[0] = C
// par[1] = p0
// par[2] = n
double power_law_Seven(double *x, double*par) 
{
	double res = par[0]*x[0]*TMath::Power( 1 + TMath::Power(x[0] / par[1] , 2) , - par[2] );
	return res ;
}



// Blast-wave
// par[0] = A, normalization constant
// par[1] = m , the mass of the studied particle (which will be fixed)
// par[2] = T , the kinetic freeze-out temperature
// par[3] = n , the velocity profile
// par[4] = beta_s
double blast_wave(double *x, double *par)
{
// we define the variables and parameters
	//p_T = x[0];
	double a = 0.0;
	double b = TMath::ATanH(par[4]);
	int n = 1000;	
	//double R = 12*TMath::Power(10, -15);
	double R = 1.0;
	double m_T = TMath::Sqrt(par[1]*par[1] + x[0]*x[0]);
// we integrate over rho
	double h = (b-a) / n;
	double z = 0;
	double *rho_0 = new double;  
	double *rho_1 = new double;
	double *rho_2 = new double;
	for(int i = 0; i <= n-1  ; i++)
	{
	    rho_0[0] = a + i*h;   // x_i
	    rho_1[0] = a + (2*i+1)*h/2;  // ( x_{i} + x_{i+1} ) / 2               
	    rho_2[0] = a + (i+1)*h;   // x_{i+1}
	    z = z        +       2*TMath::Pi()*par[0]*x[0]*m_T*R*R*TMath::Power(par[3],-1)*TMath::Power(par[4], -2 / par[3])*TMath::Power( TMath::TanH(rho_0[0] ) , (2-par[3])/par[3])*( 1- TMath::TanH(rho_0[0] )*TMath::TanH(rho_0[0] ) ) * TMath::BesselI0( x[0]*TMath::SinH(rho_0[0] ) / par[2] ) * TMath::BesselK1( m_T*TMath::CosH(rho_0[0] ) / par[2])                +     4* (                 2*TMath::Pi()*par[0]*x[0]*m_T*R*R*TMath::Power(par[3],-1)*TMath::Power(par[4], -2 / par[3])*TMath::Power( TMath::TanH(rho_1[0]) , (2-par[3])/par[3])*( 1- TMath::TanH(rho_1[0] )*TMath::TanH(rho_1[0] ) ) * TMath::BesselI0( x[0]*TMath::SinH(rho_1[0] ) / par[2] ) * TMath::BesselK1( m_T*TMath::CosH(rho_1[0]) / par[2])           )                          +                      2*TMath::Pi()*par[0]*x[0]*m_T*R*R*TMath::Power(par[3],-1)*TMath::Power(par[4], -2 / par[3])*TMath::Power( TMath::TanH(rho_2[0] ) , (2-par[3])/par[3])*( 1- TMath::TanH(rho_2[0] )*TMath::TanH(rho_2[0] ) ) * TMath::BesselI0( x[0]*TMath::SinH(rho_2[0] ) / par[2] ) * TMath::BesselK1( m_T*TMath::CosH(rho_2[0] ) / par[2])   ;
	}
	z =  h*z/6;
	return z;
}




int main()
{
	cout << " Beginning of the main() " << endl ;
	
// Creation of an input file (to retrieve the data) and of an ouput file
	TString outputfilename="result.root" ;
	TFile* OutputHisto = new TFile(outputfilename, "RECREATE");
	TFile *myFile = new TFile("HEPData-1569102768-v1-root.root");

// Retrieval of serveral histograms
	TDirectoryFile* dirFile = (TDirectoryFile*)myFile->Get("Table 5");
		TH1F* H1_pp=(TH1F*)dirFile->Get("Hist1D_y1");
		TH1F* H1_pp_e1=(TH1F*)dirFile->Get("Hist1D_y1_e1");
		TH1F* H1_pp_E1=(TH1F*)dirFile->Get("Hist1D_y1_e2");
// Changing the names of axis and titles
		//H1_pp->SetNameTitle(" Table 4 phi-meson" , "pT-distributions of phi-meson measured in p-pbar collisions at sNN = 5.02 TeV." );
		H1_pp->SetXTitle("p_{T} [GeV/c]");
		H1_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");

// Retrieval of the scale of each bin in order to normalize it for the fit
	double Scale_Histo[11] ;
	Scale_Histo[1]=H1_pp->Integral("width") ;
// Retrieval of the uncertainties for each bin and adding it into the histogram
	int Nbinx = H1_pp->GetNbinsX();
	for(int i = 0; i <= Nbinx  ; i++){
		H1_pp->SetBinError(i, sqrt( pow(H1_pp_e1->GetBinContent(i),2) + pow(H1_pp_E1->GetBinContent(i),2) ) ) ;
	}		
	H1_pp->SetAxisRange(0,20,"X");
	double a= H1_pp->GetBinLowEdge(1);
	double b= H1_pp->GetBinLowEdge(Nbinx);
	b += H1_pp->GetBinWidth(Nbinx);

// Definition of the fitting functions
			
	TF1 *func1 = new TF1("expo_law",expo_law,0,b,3);
		func1->SetParameter(0,Scale_Histo[1]);
		func1->SetParameter(1,0.938);
		func1->SetParLimits(1,0.938,0.938);
		func1->SetParameter(2,0.2);
		func1->SetParNames("C","m","T");
		func1->SetLineColor(kViolet);


	TF1 *func2 = new TF1("boltzmann",boltzmann,0,b,3);
		func2->SetParameter(0,Scale_Histo[1]);
		func2->SetParameter(1,0.938);
		func2->SetParLimits(1,0.938,0.938);
		func2->SetParameter(2,0.2);
		func2->SetParNames("C", "m", "T");	
		func2->SetLineColor(kGreen);


	TF1 *func3 = new TF1("levy",levy,0,b,4);
		func3->SetParameter(0,Scale_Histo[1]);
		func3->SetParameter(1,0.938);
		func3->SetParLimits(1,0.938,0.938);
		func3->SetParameter(2,0.6);
		func3->SetParameter(3,4);
		func3->SetParLimits(3,2,8);
		func3->SetParNames("C","m","T","n");	
		//func3->SetParLimits(3, 1.01 , 10);
		func3->SetLineColor(kCyan);

	TF1 *func4 = new TF1("power_law_Five", power_law_Five,0,b,3);
		func4->SetParameter(0,1);
		func4->SetParameter(1,1);
		func4->SetParameter(2,4);
		func4->SetParNames("C","p0","n");	
		func4->SetLineColor(kMagenta);

	TF1 *func5 = new TF1("power_law_Seven", power_law_Seven ,0,b,3);
		func5->SetParameter(0,1);
		func5->SetParameter(1,1);
		func5->SetParameter(2,1);
		func5->SetParNames("C","p0","n");	
		func5->SetLineColor(kRed);


	TF1 *func6 = new TF1("blast_wave",blast_wave,0,b,5);
		func6->SetParameter(0,100000);
		func6->SetParLimits(0,1,1);
		func6->SetParameter(1,0.938);
		func6->SetParLimits(1, 0.938,0.938);
		func6->SetParameter(2,0.1);
		func6->SetParameter(3,2);
		func6->SetParameter(4,0.5);
		func6->SetParNames("A","m","T","n","beta_s");
		func6->SetLineColor(kBlue);


// Creation of the canvas
	TCanvas *T1 = new TCanvas("Canvas 1" , " CANVAS_1_PP1 " , 200 , 200 );
	T1->SetGrid();
	T1->SetLogy();
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(1111);
	
// Saving the canvas in the ouput-histo-file i.e. result.root
	OutputHisto->cd();

// We fit the data with the different fitting functions
	T1->cd();
	H1_pp->Draw("");
	//func6->Draw();
	H1_pp->Fit("expo_law","+ I  ","SAMES HIST",0,b);  
		double chi2_expo = func1->GetChisquare();
		double ndf_expo = func1->GetNDF();
		double chindf_expo = chi2_expo/ndf_expo ;
		cout <<" Chi square expo law = "<< chi2_expo << endl;
		cout <<" NDF expo law = "<< ndf_expo << endl; 
		cout <<" Chi2 / ndf = " << chindf_expo << endl;


	H1_pp->Fit("boltzmann","+ I ","SAMES HIST",0,b);
		double chi2_bol = func2->GetChisquare();
		double ndf_bol = func2->GetNDF();
		double chindf_bol = chi2_bol/ndf_bol ;
		cout <<" Chi square boltz = "<< chi2_bol << endl;
		cout <<" NDF boltz = "<< ndf_bol << endl; 
		cout <<" Chi2 / ndf = " << chindf_bol << endl;

	H1_pp->Fit("levy","I +","SAMES HIST",0,b);
		double chi2_levy = func3->GetChisquare();
		double ndf_levy = func3->GetNDF();
		double chindf_levy = chi2_levy/ndf_levy ;
		cout <<" Chi square levy = "<< chi2_levy << endl;
		cout <<" NDF levy = "<< ndf_levy << endl; 
		cout <<" Chi2 / ndf = " << chindf_levy << endl;

	//H1_pp->Fit("power_law_Five"," I +","SAMES HIST",0,b);
	//H1_pp->Fit("power_law_Seven","+ I ","SAMES HIST",0,b);

	//H1_pp -> Fit("blast_wave", " I V","SAMES HIST",a,b);

	//H1_pp->GetFunction("expo_law")->SetLineColor(kViolet);
	//H1_pp->GetFunction("boltzmann")->SetLineColor(kGreen);
	//H1_pp->GetFunction("levy")->SetLineColor(kCyan);
	//H1_pp->GetFunction("power_law_Five")->SetLineColor(kMagenta);
	//H1_pp->GetFunction("power_law_Seven")->SetLineColor(kRed);

	//H1_pp->GetFunction("blast_wave")->SetLineColor(kBlue);
	
// Adding the results of the fit on the canvas
	TLatex latex;
		latex.SetNDC();
		latex.SetTextSize(0.020);
		latex.SetTextAlign(13);  //align at top
		latex.SetTextFont(60);
	   
	 
		latex.DrawLatex(0.5,0.65, TString::Format("Parametre Levy C= %g", func3->GetParameter(0)));
		latex.DrawLatex(0.5,0.625, TString::Format("Parametre Levy m= %g GeV/c", func3->GetParameter(1)));
		latex.DrawLatex(0.5,0.60, TString::Format("Parametre Levy T= %g GeV", func3->GetParameter(2)));
		latex.DrawLatex(0.5,0.575, TString::Format("Parametre Levy n= %g ", func3->GetParameter(3)));
		latex.DrawLatex(0.5,0.550, TString::Format("CHI2/NDF Levy = %g ", func3->GetChisquare()/func3->GetNDF()));

	  
		latex.DrawLatex(0.5,0.50, TString::Format("Parametre expo C= %g", func1->GetParameter(0)));
		latex.DrawLatex(0.5,0.475, TString::Format("Parametre expo m_0= %g GeV/c", func1->GetParameter(1)));
		latex.DrawLatex(0.5,0.45, TString::Format("Parametre expo T= %g GeV/c", func1->GetParameter(2)));
		latex.DrawLatex(0.5,0.425, TString::Format("CHI2/NDF expo = %g ", func1->GetChisquare()/func1->GetNDF()));
	 
		latex.DrawLatex(0.5,0.375, TString::Format("Parametre boltzmann C= %g", func2->GetParameter(0)));
		latex.DrawLatex(0.5,0.350, TString::Format("Parametre boltzmann m_0= %g GeV/c", func2->GetParameter(1)));
		latex.DrawLatex(0.5,0.325, TString::Format("Parametre boltzmann T= %g GeV/c", func2->GetParameter(2)));
		latex.DrawLatex(0.5,0.30, TString::Format("CHI2/NDF boltzmann = %g ", func2->GetChisquare()/func2->GetNDF()));
	T1->Update();
// Configuration of the legend
	auto legend = new TLegend(0.4,0.7,0.6,0.9);
		//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
		legend->AddEntry(H1_pp,"Data","lep");
		legend->AddEntry(func1,"Exponential law fit","l");
		legend->AddEntry(func2,"Boltzmann law fit","l");
		legend->AddEntry(func3,"Levy-Tsallis fit","l");
		//legend->AddEntry(func4, "Power law fit from [5]", "l");
		//legend->AddEntry(func5, "Power law fit from [7]", "l");
		//legend->AddEntry(func6,"Blast-wave model","l");
		legend->Draw("SAMES");
	
	T1->SetLogy();	 // log scale for the y axis (better visualization of the fit)		
	T1->Write();
	OutputHisto->Close();
	cout << " End of the main() " << endl ;
}
