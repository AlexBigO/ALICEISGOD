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
#include <TGraphPainter.h>
#include <TGraphErrors.h>
#include <TMinuit.h>

#if defined(__GLIBCXX__)
#pragma link C++ class MyOtherClass;
#endif


using namespace std;



TH1F *hist = 0 ;



double expo( double x , double *par)
{
	double m_T = sqrt(x*x + par[1]*par[1]);
	double res = par[0]*TMath::Power( par[2]*(par[1]+par[2]), -1)*x*exp(- (m_T - par[1]) / par[2]);
	return res;
}

		
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
double blast_wave(double *x , double *par)
{
    // x[0] = p_T  here
// we define the variables and parameters
	double a = 0.0;
	double R = 12 ;
	double b = R;
	int n = 1000;	
	double m_T = TMath::Sqrt(par[1]*par[1] + x[0]*x[0]);
// we integrate over r
	double h = (b-a) / n;
	double z = 0;
	double r_0 , r_1 , r_2;  
	for(int i = 0; i <= n-1  ; i++)
	{
	    r_0 = a + i*h;   // x_i
	    r_1 = a + (2*i+1)*h/2;  // ( x_{i} + x_{i+1} ) / 2               
	    r_2 = a + (i+1)*h;   // x_{i+1}
	    z +=  2*TMath::Pi()*par[0]*x[0]*m_T*
	    (   r_0*TMath::BesselI0( x[0]*TMath::SinH( TMath::ATanH( TMath::Power( r_0 / R, par[3] ) * par[4]  )) / par[2] ) * TMath::BesselK1( m_T*TMath::CosH( TMath::ATanH(  TMath::Power( r_0 / R, par[3] ) * par[4] ) ) / par[2])     + 4*r_1*TMath::BesselI0( x[0]*TMath::SinH( TMath::ATanH( TMath::Power( r_1 / R, par[3] ) * par[4]  )) / par[2] ) * TMath::BesselK1( m_T*TMath::CosH( TMath::ATanH(  TMath::Power( r_1 / R, par[3] ) * par[4] ) ) / par[2])     + r_2*TMath::BesselI0( x[0]*TMath::SinH( TMath::ATanH( TMath::Power( r_2 / R, par[3] ) * par[4]  )) / par[2] ) * TMath::BesselK1( m_T*TMath::CosH( TMath::ATanH(  TMath::Power( r_2 / R, par[3] ) * par[4] ) ) / par[2])  );
	}
	z =  h*z/6;
	return z;
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
	f = 0.;
	for(int i=1;i<=hist->GetNbinsX();i++) 
	{
		double *x = new double ;
		x[0] = hist->GetBinCenter(i);
		double measure = hist->GetBinContent(i);
		double error = hist->GetBinError(i);
		double func = levy(x,par);
		double delta = (func - measure)/error;
		f += delta*delta;
	}
}


// Simpson integration method
// Input parameters: 
// f : the function we want to integrate; this is a user-defined function in our case (expo_law, boltzmann, ...)
// par : the parameters of f; these will be the output parameters of the fit via f
// [a,b] the integration range
// n the number of points we consider such as (b-a)/n is the step of the integration
double simpson(double f(double*,double*), double *par, double a, double b, int n)   
{
    double h = (b-a) / n;    // definition of the integration step
    double z = 0;           // initialize the variable
    // we define different pointers that we will use as input parameters for our f function
    double *x = new double;  
    double *x1 = new double;
    double *x2 = new double;
    for(int i = 0; i <= n-1  ; i++)
    {
        x[0] = a + i*h;   // x_i
        x1[0] = a + (2*i+1)*h/2;  // ( x_{i} + x_{i+1} ) / 2
        x2[0] = a + (i+1)*h;   // x_{i+1}
        z = z + f(x, par) + 4*f(x1,par) + f(x2 , par);
    }
    return h*z/6;
}




// Data / Model fit 
// We compute the ratio data/model for each bin of a given histogram hist
// and store the values in tables
void Data_Model(double X[], double Y[], double Xerr[], double Yerr[], double func(double *,double *), double *par)
{
	int n = 1000;
	double data_value , model_value , model_error;
	double integral, low_edge, up_edge, bin_width ;
	int Nbinx = hist->GetNbinsX();
	for(int i=1;i<=Nbinx;i++) 
	{
	// for each bin, we get its content
		data_value = hist->GetBinContent(i);
	// we compute the mean value of the bin via the model
	// i.e. we comute the integral of the model on the bin's width and divide by the bin's width
		low_edge = hist->GetBinLowEdge(i);
		up_edge = hist->GetBinLowEdge(i);
		bin_width = hist->GetBinWidth(i);
		up_edge += bin_width;
		integral = simpson(func,par,low_edge,up_edge,n);
		model_value = integral / bin_width;
		model_error = TMath::Power((up_edge - low_edge)/n , 4) / bin_width;  // simpson error of order h^4
	// then we compute the ratio data/model fit
		Y[i-1] = data_value / model_value;
	// we compute the uncertainty on this Y value via propagation of uncertainty
	// considering the data and model uncertainties are uncorrelated 
		Yerr[i-1] = TMath::Sqrt(  TMath::Power(hist->GetBinError(i) / model_value,2) + TMath::Power( data_value* model_error /(model_value*model_value) ,2) ) ; 
	// we also get the center value of the bin for the future plot
		X[i-1] = hist->GetBinCenter(i);
		Xerr[i-1] = 0.5*bin_width;  // Xerr will represent the width of the bin (it is not an error strictly speaking)
	}

}

int main()
{
	cout << " début " << endl ;

// Création du fichier root de sortie et récupération du fichier root données
	TString outputfilename="result.root" ;
	TFile* OutputHisto = new TFile(outputfilename, "RECREATE");
	TFile *myFile = new TFile("HEPData-1569102768-v1-root.root");

	TString contours_1_sigma_output="contours_1_sigma.root";
	TFile *Contours_1_sigma_output = new TFile(contours_1_sigma_output, "RECREATE");
	TString contours_2_sigma_output="contours_2_sigma.root";
	TFile *Contours_2_sigma_output = new TFile(contours_2_sigma_output, "RECREATE");

// Récupération des différents histogrammes 
	TDirectoryFile* dirFile = (TDirectoryFile*)myFile->Get("Table 6");
		TH1F* H1_pp=(TH1F*)dirFile->Get("Hist1D_y1");
		TH1F* H1_pp_e1=(TH1F*)dirFile->Get("Hist1D_y1_e1");
		TH1F* H1_pp_E1=(TH1F*)dirFile->Get("Hist1D_y1_e2");
// Changement des noms d'axes et titres
		//H1_pp->SetNameTitle(" Table 4 phi-meson" , "pT-distributions of phi-meson measured in p-pbar collisions at sNN = 5.02 TeV." );
		H1_pp->SetXTitle("p_{T} [GeV/c]");
		H1_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");


//Récupération du scaling en prépartion du fit
	double Scale_Histo[11] ;
	Scale_Histo[1]=H1_pp->Integral("width") ;
// Compilations des différentes erreurs et ajout sur l'histogramme PT
	int Nbinx = H1_pp->GetNbinsX();
	for(int i = 0; i <= Nbinx  ; i++){
		H1_pp->SetBinError(i, sqrt( pow(H1_pp_e1->GetBinContent(i),2) + pow(H1_pp_E1->GetBinContent(i),2) ) ) ;
	}		
	H1_pp->SetAxisRange(0,20,"X");
	double a= H1_pp->GetBinLowEdge(1);
	double b= H1_pp->GetBinLowEdge(Nbinx);
	b += H1_pp->GetBinWidth(Nbinx);




 	hist = (TH1F *)dirFile->Get("Hist1D_y1");
 
	TH1F* hist_e1=(TH1F*)dirFile->Get("Hist1D_y1_e1");
	TH1F* hist_E1=(TH1F*)dirFile->Get("Hist1D_y1_e2");
	hist->SetXTitle("p_{T} [GeV/c]");
	hist->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");

	for(int i = 0; i <= Nbinx  ; i++)
	{
		hist->SetBinError(i, sqrt( pow(hist_e1->GetBinContent(i),2) + pow(hist_E1->GetBinContent(i),2) ) ) ;
	}		
	TMinuit *gMinuit = new TMinuit(4);
	gMinuit->SetFCN(fcn);
	gMinuit->DefineParameter(0, "C", 0.1, 0.01, 0, 100);
	//gMinuit->FixParameter(0);
	gMinuit->DefineParameter(1, "m_0",0.938, 0.01, 0.1,1);
	gMinuit->FixParameter(1);
	gMinuit->DefineParameter(2, "T", 0.05, 0.001, 0, 0.3);
	gMinuit->DefineParameter(3, "n", 7.7938, 0.1, 0, 10);
	gMinuit->FixParameter(3);	
	//gMinuit->DefineParameter(4, "beta_s", 0.5, 0.01, 0, 1);
	gMinuit->Command("MIGRAD");
	gMinuit->Command("MINOS");
	double par[4],err[4];
	for(int i=0;i<4;i++){
		gMinuit->GetParameter(i,par[i],err[i]);
	}
	TH1F* curve = new TH1F("curve","curve",hist->GetNbinsX()*5,0,b);
	for(int i=1;i<=curve->GetNbinsX();i++) 
	{		
		double *x = new double ;
		x[0] = curve->GetBinCenter(i);
		double f = levy(x,par);
		curve->SetBinContent(i,f);
	} 
	curve->SetLineWidth(3);
	hist->Draw();
	curve->Draw("csame");
	
//-------------------------------------------------------------------------------------------------
// DATA / MODEL FIT

	double X[Nbinx];  // x axis of the plot
	double Y[Nbinx];  // y axis of the plot
	double Xerr[Nbinx];  // we take into account the width of a bin (this is not an error strictly speaking!)
	double Yerr[Nbinx];  // error on y axis
	
	Data_Model(X,Y,Xerr,Yerr,levy,par);
	
	//cout << "Yerr after = " <<  Yerr[Nbinx-1] << endl;
	
	TCanvas *c2 = new TCanvas("c2","c2",200,10,600,400);
	c2->SetGrid();
	// we create the graph to plot Data/Model fit with rectangle errors
	TGraphErrors *dataModel = new TGraphErrors(20,X,Y,Xerr,Yerr);
	dataModel->SetLineColor(2);
	dataModel->SetLineWidth(2);
	dataModel->SetMarkerColor(4);
	dataModel->SetMarkerSize(1);
	dataModel->SetMarkerStyle(21);
	dataModel->SetTitle(" ;p_T (GeV/c); Data / model fit");
	dataModel->SetLineWidth(2);
	dataModel->Draw("A5");   // A has to be there (if not, nothing appears on the canvas) & 5 is to plot error rectangles
	dataModel->Draw("PX");   // P is to put the chosen marker & X to remove the error bars (leave only the rectangles)
	
	// we plot a lin at y=1 to show the deviation to Data = Model fit
	TLine *line = new TLine(0,1,b,1);
	line->SetLineColor(30);
	line->SetLineWidth(3);
	line->SetLineStyle(2);
	line->Draw("SAME");
//----------------------------------------------------------------------------------------------------

/*	
   auto c4 = new TCanvas("c4","c4",200,10,600,400);
   double x[] = {0, 1, 2, 3, 4};
   double y[] = {0, 2, 4, 1, 3};
   double ex[] = {0.1, 0.2, 0.3, 0.4, 0.5};
   double ey[] = {1, 0.5, 1, 0.5, 1};
   auto ge = new TGraphErrors(5, x, y, ex, ey);
   ge->Draw("ap"); */	
	
/*
   //Get contour for parameter 2 versus parameter 4  for ERRDEF=2
   	Contours_2_sigma_output -> cd();  // we save the 2-sigma contours 
   	gMinuit->SetErrorDef(4); //note 4 and not 2!
   	TGraph *gr2 = (TGraph*)gMinuit->Contour(100,2,0);
   	//gr2->SetFillColor(42);
   	//gr2->Draw("alf");
   	gr2->Write("Graph");
	Contours_2_sigma_output -> Close();

   
   //Get contour for parameter 2 versus parameter 4 for ERRDEF=1
	Contours_1_sigma_output -> cd();  // we save the 1-sigma contours 
   	gMinuit->SetErrorDef(1);
   	TGraph *gr1 = (TGraph*)gMinuit->Contour(100,2,0);
   	//gr1->SetFillColor(38);
   	//gr1->Draw("lf");
        gr1->Write("Graph");
        Contours_1_sigma_output -> Close();
	 */
	OutputHisto->Close();
	cout << " fin " << endl ;
}
