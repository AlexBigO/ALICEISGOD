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
        //cout << "Simpson iteration " << i << endl;
    }
    return h*z/6;
}



/*
//Simson avec pointeurs
double *simpson(double f(double*,double*), double *par, double a, double b, int n)   
{
    double h = (b-a) / n;    // definition of the integration step
    //double z = 0;           // initialize the variable
    double *z = new double;
    // we define different pointers that we will use as input parameters for our f function
    double *x = new double;  
    double *x1 = new double;
    double *x2 = new double;
    for(int i = 0; i <= n-1  ; i++)
    {
        x[0] = a + i*h;   // x_i
        x1[0] = a + (2*i+1)*h/2;  // ( x_{i} + x_{i+1} ) / 2
        x2[0] = a + (i+1)*h;   // x_{i+1}
        z[0] = z[0] + f(x, par) + 4*f(x1,par) + f(x2 , par);
    }
    z[0] = h*z[0]/6;
    return z;
}
*/


/*
double integrand(double *x , double *par)
{
	double res = par[0]*x[0] + par[1] ;
	return res;
}

double test(double *x, double *par)
{
    double R = 1.0 ;
    par[1] = x[0];  // we explicitely show that p_T is the variable
  // we change the status of p_T (it was a parameter and it becomes the variable of our new function)
    double res = simpson(integrand, par , 0 , R , 1000);
    return res;
}
*/

/*
// par[0] = A, normalization constant
// par[1] = m , the mass of the studied particle (which will be fixed)
// par[2] = T , the kinetic freeze-out temperature
// par[3] = n , the velocity profile
// par[4] = beta_s
// par[5] = p_T in the function "integrand" (we do this because the integration variable in integrand is r, not p_T)
double integrand(double *x, double *par)
{
    // x[0] = rho   here
    double R = 10*TMath::Power(10, -15);
    double m_T = TMath::Sqrt(par[1]*par[1] + par[5]*par[5]);
    double res = 2*TMath::Pi()*par[0]*par[5]*m_T*R*R*TMath::Power(par[3],-1)*TMath::Power(par[4], -2 / par[3])*TMath::Power( TMath::TanH(x[0]) , (2-par[3])/par[3])*( 1- TMath::TanH(x[0])*TMath::TanH(x[0]) ) * TMath::BesselI0( par[5]*TMath::SinH(x[0]) / par[2] ) * TMath::BesselK1( m_T*TMath::CosH(x[0]) / par[2]);
    cout << "par[0] = " << par[0] << endl; 
    cout << "par[1] = " << par[1] << endl; 
    cout << "par[2] = " << par[2] << endl; 
    cout << "par[3] = " << par[3] << endl; 
    cout << "par[4] = " << par[4] << endl;  
    //cout << "par[5] = " << par[5] << endl;
    return res;
}


double integration(double *x,double *par)
{
    double *par1 = new double;
    par1 = par;
    par1[5] = 1.0;  // but in this function p_T is not a parameter, so we have to create a new parameter for the integrand function
    double rho_max = TMath::ATanH(par[4]);  // we define the upper limit of the integration
    double res = simpson(integrand, par1 , 0 , rho_max , 100);  // we integrate the integrand over [0,argtanh beta_s]
    return res;
}



// Then we define our blast-wave function
double blast_wave(double *x, double *par)
{
    //double *ptr = new double;
    //double *res = new double;
    double res = integration(x,par);
    //double *par_new = new double;
    //par_new = ptr;
    par[5] = x[0];
    return res;
    double *par1 = new double;
    par1 = par;
    par1[5] = x[0];  // but in this function p_T is not a parameter, so we have to create a new parameter for the integrand function
    double rho_max = TMath::ATanH(par[4]);  // we define the upper limit of the integration
    double res = simpson(integrand, par1 , 0 , rho_max , 100);  // we integrate the integrand over [0,argtanh beta_s]
    return res;
}

*/





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




int main(){
cout << " début " << endl ;

		// Création du fichier root de sortie et récupération du fichier root données
TString outputfilename="result.root" ;
TFile* OutputHisto = new TFile(outputfilename, "RECREATE");
TFile *myFile = new TFile("HEPData-1569102768-v1-root.root");


		
		// Récupération des différents histogrammes 
TDirectoryFile* dirFile = (TDirectoryFile*)myFile->Get("Table 5");


	TH1F* H1_pp=(TH1F*)dirFile->Get("Hist1D_y1");
	/*
	TH1F* H2_pp=(TH1F*)dirFile->Get("Hist1D_y2");
	TH1F* H3_pp=(TH1F*)dirFile->Get("Hist1D_y3");
	TH1F* H4_pp=(TH1F*)dirFile->Get("Hist1D_y4");
	TH1F* H5_pp=(TH1F*)dirFile->Get("Hist1D_y5");
	TH1F* H6_pp=(TH1F*)dirFile->Get("Hist1D_y6");
	TH1F* H7_pp=(TH1F*)dirFile->Get("Hist1D_y7");
	TH1F* H8_pp=(TH1F*)dirFile->Get("Hist1D_y8");
	TH1F* H9_pp=(TH1F*)dirFile->Get("Hist1D_y9");
	TH1F* H10_pp=(TH1F*)dirFile->Get("Hist1D_y10");
*/
	TH1F* H1_pp_e1=(TH1F*)dirFile->Get("Hist1D_y1_e1");
	/*
	TH1F* H2_pp_e2=(TH1F*)dirFile->Get("Hist1D_y2_e1");
	TH1F* H3_pp_e3=(TH1F*)dirFile->Get("Hist1D_y3_e1");
	TH1F* H4_pp_e4=(TH1F*)dirFile->Get("Hist1D_y4_e1");
	TH1F* H5_pp_e5=(TH1F*)dirFile->Get("Hist1D_y5_e1");
	TH1F* H6_pp_e6=(TH1F*)dirFile->Get("Hist1D_y6_e1");
	TH1F* H7_pp_e7=(TH1F*)dirFile->Get("Hist1D_y7_e1");
	TH1F* H8_pp_e8=(TH1F*)dirFile->Get("Hist1D_y8_e1");
	TH1F* H9_pp_e9=(TH1F*)dirFile->Get("Hist1D_y9_e1");
	TH1F* H10_pp_e10=(TH1F*)dirFile->Get("Hist1D_y10_e1");
	*/
	TH1F* H1_pp_E1=(TH1F*)dirFile->Get("Hist1D_y1_e2");
	/*
	TH1F* H2_pp_E2=(TH1F*)dirFile->Get("Hist1D_y2_e2");
	TH1F* H3_pp_E3=(TH1F*)dirFile->Get("Hist1D_y3_e2");
	TH1F* H4_pp_E4=(TH1F*)dirFile->Get("Hist1D_y4_e2");
	TH1F* H5_pp_E5=(TH1F*)dirFile->Get("Hist1D_y5_e2");
	TH1F* H6_pp_E6=(TH1F*)dirFile->Get("Hist1D_y6_e2");
	TH1F* H7_pp_E7=(TH1F*)dirFile->Get("Hist1D_y7_e2");
	TH1F* H8_pp_E8=(TH1F*)dirFile->Get("Hist1D_y8_e2");
	TH1F* H9_pp_E9=(TH1F*)dirFile->Get("Hist1D_y9_e2");
	TH1F* H10_pp_E10=(TH1F*)dirFile->Get("Hist1D_y10_e2");
	*/
	/*
	TH1F* H1_pp_SU1=(TH1F*)dirFile->Get("Hist1D_y1_e3");
	TH1F* H2_pp_SU2=(TH1F*)dirFile->Get("Hist1D_y2_e3");
	TH1F* H3_pp_SU3=(TH1F*)dirFile->Get("Hist1D_y3_e3");
	TH1F* H4_pp_SU4=(TH1F*)dirFile->Get("Hist1D_y4_e3");
	TH1F* H5_pp_SU5=(TH1F*)dirFile->Get("Hist1D_y5_e3");
	TH1F* H6_pp_SU6=(TH1F*)dirFile->Get("Hist1D_y6_e3");
	TH1F* H7_pp_SU7=(TH1F*)dirFile->Get("Hist1D_y7_e3");
	TH1F* H8_pp_SU8=(TH1F*)dirFile->Get("Hist1D_y8_e3");
	TH1F* H9_pp_SU9=(TH1F*)dirFile->Get("Hist1D_y9_e3");
	TH1F* H10_pp_SU10=(TH1F*)dirFile->Get("Hist1D_y10_e3");
	*/
	
		// Changement des noms d'axes et titres
	//H1_pp->SetNameTitle(" Table 4 phi-meson" , "pT-distributions of phi-meson measured in p-pbar collisions at sNN = 5.02 TeV." );
	H1_pp->SetXTitle("p_{T} [GeV/c]");
	H1_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	/*
	H2_pp->SetNameTitle(" Distribution Pt 2 " , " PT distribution des pions lors de collisions Pb-PB " );
	H2_pp->SetXTitle("p_{T} [GeV/c]");
	H2_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H3_pp->SetNameTitle(" Distribution Pt 3 " , " PT distribution des pions lors de collisions Pb-PB " );
	H3_pp->SetXTitle("p_{T} [GeV/c]");
	H3_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H4_pp->SetNameTitle(" Distribution Pt 4 " , " PT distribution des pions lors de collisions Pb-PB " );
	H4_pp->SetXTitle("p_{T} [GeV/c]");
	H4_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H5_pp->SetNameTitle(" Distribution Pt 5 " , " PT distribution des pions lors de collisions Pb-PB " );
	H5_pp->SetXTitle("p_{T} [GeV/c]");
	H5_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H6_pp->SetNameTitle(" Distribution Pt 6 " , " PT distribution des pions lors de collisions Pb-PB " );
	H6_pp->SetXTitle("p_{T} [GeV/c]");
	H6_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H7_pp->SetNameTitle(" Distribution Pt 7 " , " PT distribution des pions lors de collisions Pb-PB " );
	H7_pp->SetXTitle("p_{T} [GeV/c]");
	H7_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H8_pp->SetNameTitle(" Distribution Pt 8 " , " PT distribution des pions lors de collisions Pb-PB " );
	H8_pp->SetXTitle("p_{T} [GeV/c]");
	H8_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H9_pp->SetNameTitle(" Distribution Pt 9 " , " PT distribution des pions lors de collisions Pb-PB " );
	H9_pp->SetXTitle("p_{T} [GeV/c]");
	H9_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H10_pp->SetNameTitle(" Distribution Pt 10 " , " PT distribution des pions lors de collisions Pb-PB " );
	H10_pp->SetXTitle("p_{T} [GeV/c]");
	H10_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	
	*/



		//Scaling en prépartion du fit

double Scale_Histo[11] ;
Scale_Histo[1]=H1_pp->Integral("width") ;
/*
Scale_Histo[2]=H2_pp->Integral() ;
Scale_Histo[3]=H3_pp->Integral() ;
Scale_Histo[4]=H4_pp->Integral() ;
Scale_Histo[5]=H5_pp->Integral() ;
Scale_Histo[6]=H6_pp->Integral() ;
Scale_Histo[7]=H7_pp->Integral() ;
Scale_Histo[8]=H8_pp->Integral() ;
Scale_Histo[9]=H9_pp->Integral() ;
Scale_Histo[10]=H10_pp->Integral() ;
*/

cout << " Scale Histo 1 = "<< H1_pp->Integral("width") << endl;

H1_pp->Scale(1/H1_pp->Integral("width"));

cout << " Scale after = " << H1_pp->Integral("width") << endl;



		// Compilations des différentes erreurs et ajout sur l'histogramme PT

int Nbinx = H1_pp->GetNbinsX();
cout << " number of bin = " << Nbinx << endl;


for(int i = 0; i <= Nbinx  ; i++){
	H1_pp->SetBinError(i, sqrt( pow(H1_pp_e1->GetBinContent(i),2) + pow(H1_pp_E1->GetBinContent(i),2) ) / Scale_Histo[1] ) ;
	/*
	H2_pp->SetBinError(i, H2_pp_e2->GetBinContent(i) + H2_pp_E2->GetBinContent(i) ) ;
	H3_pp->SetBinError(i, H3_pp_e3->GetBinContent(i) + H3_pp_E3->GetBinContent(i) ) ;
	H4_pp->SetBinError(i, H4_pp_e4->GetBinContent(i) + H4_pp_E4->GetBinContent(i) ) ;
	H5_pp->SetBinError(i, H5_pp_e5->GetBinContent(i) + H5_pp_E5->GetBinContent(i) ) ;
	H6_pp->SetBinError(i, H6_pp_e6->GetBinContent(i) + H6_pp_E6->GetBinContent(i) ) ;
	H7_pp->SetBinError(i, H7_pp_e7->GetBinContent(i) + H7_pp_E7->GetBinContent(i) ) ;
	H8_pp->SetBinError(i, H8_pp_e8->GetBinContent(i) + H8_pp_E8->GetBinContent(i) ) ;
	H9_pp->SetBinError(i, H9_pp_e9->GetBinContent(i) + H9_pp_E9->GetBinContent(i) ) ;
	H10_pp->SetBinError(i, H10_pp_e10->GetBinContent(i) + H10_pp_E10->GetBinContent(i)  ) ;
	*/
	
}



		// Changement de range sur l'axe X
		/*
H2_pp->Scale(1/H2_pp->Integral());
H3_pp->Scale(1/H3_pp->Integral());
H4_pp->Scale(1/H4_pp->Integral());
H5_pp->Scale(1/H5_pp->Integral());
H6_pp->Scale(1/H6_pp->Integral());
H7_pp->Scale(1/H7_pp->Integral());
H8_pp->Scale(1/H8_pp->Integral());
H9_pp->Scale(1/H9_pp->Integral());
H10_pp->Scale(1/H10_pp->Integral());
*/

		
H1_pp->SetAxisRange(0,20,"X");
/*
H2_pp->SetAxisRange(0,5,"X");
H3_pp->SetAxisRange(0,5,"X");
H4_pp->SetAxisRange(0,5,"X");
H5_pp->SetAxisRange(0,5,"X");
H6_pp->SetAxisRange(0,5,"X");
H7_pp->SetAxisRange(0,5,"X");
H8_pp->SetAxisRange(0,5,"X");
H9_pp->SetAxisRange(0,5,"X");
H10_pp->SetAxisRange(0,5,"X");
*/


double a= H1_pp->GetBinLowEdge(1);
double b= H1_pp->GetBinLowEdge(Nbinx);
b += H1_pp->GetBinWidth(Nbinx);
cout << " check bin up edge " << b << endl;
		// Fit par rapport à la distribution de boltzmann
		
TF1 *func1 = new TF1("expo_law",expo_law,0,b,3);
func1->SetParameter(0,1);
func1->SetParameter(1,0.938);
func1->SetParLimits(1, 0.938, 0.938);
func1->SetParameter(2,0.1);
func1->SetParNames("C","m","T");
func1->SetLineColor(kViolet);


TF1 *func2 = new TF1("boltzmann",boltzmann,0,b,3);
func2->SetParameter(0,0.1);
func2->SetParameter(1,0.938);
func2->SetParLimits(1, 0.938,0.938);
func2->SetParameter(2,0.1);
func2->SetParNames("C", "m", "T");	
func2->SetLineColor(kGreen);


TF1 *func3 = new TF1("levy",levy,0,b,4);
func3->SetParameter(0,1);
func3->SetParameter(1,0.938);
func3->SetParLimits(1, 0.938,0.938);
func3->SetParameter(2,0.1);
func3->SetParameter(3,3);
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


		// Création des canvas
TCanvas *T1 = new TCanvas("Canvas 1" , " CANVAS_1_PP1 " , 200 , 200 );
T1->SetGrid();
/*
TCanvas *T2 = new TCanvas("Canvas 2" , " CANVAS_2_PP1 " , 200 , 200 );
T2->SetGrid();
TCanvas *T3 = new TCanvas("Canvas 3" , " CANVAS_3_PP1 " , 200 , 200 );
T3->SetGrid();
TCanvas *T4 = new TCanvas("Canvas 4" , " CANVAS_4_PP1 " , 200 , 200 );
T4->SetGrid();
TCanvas *T5 = new TCanvas("Canvas 5" , " CANVAS_5_PP1 " , 200 , 200 );
T5->SetGrid();
TCanvas *T6 = new TCanvas("Canvas 6" , " CANVAS_6_PP1 " , 200 , 200 );
T6->SetGrid();
TCanvas *T7 = new TCanvas("Canvas 7" , " CANVAS_7_PP1 " , 200 , 200 );
T7->SetGrid();
TCanvas *T8 = new TCanvas("Canvas 8" , " CANVAS_8_PP1 " , 200 , 200 );
T8->SetGrid();
TCanvas *T9 = new TCanvas("Canvas 9" , " CANVAS_9_PP1 " , 200 , 200 );
T9->SetGrid();
TCanvas *T10 = new TCanvas("Canvas 10" , " CANVAS_10_PP1 " , 200 , 200 );
T10->SetGrid();
*/


//TCanvas *D1 = new TCanvas("Canvas Distribution" , " CANVAS_DISTR " , 2000 , 2000 );
//D1->SetGrid();



   
gStyle->SetOptFit(0);
gStyle->SetOptStat(1111);
	

		// Sauvegarde des canvas dans l'output-histo-file càd notre result.root
		// Appels des différents fits
OutputHisto->cd();


T1->cd();
H1_pp->Draw("][ E1");
// We fit the data with the different fitting functions
//func6->Draw();
/*
H1_pp->Fit("expo_law","+ I L ","SAMES HIST",a,b);  
double chi2_expo = func1->GetChisquare();
double ndf_expo = func1->GetNDF();
double chindf_expo = chi2_expo/ndf_expo ;
cout <<" Chi square expo law = "<< chi2_expo << endl;
cout <<" NDF expo law = "<< ndf_expo << endl; 
cout <<" Chi2 / ndf = " << chindf_expo << endl;


H1_pp->Fit("boltzmann","+ I L","SAMES HIST",a,b);
double chi2_bol = func2->GetChisquare();
double ndf_bol = func2->GetNDF();
double chindf_bol = chi2_bol/ndf_bol ;
cout <<" Chi square boltz = "<< chi2_bol << endl;
cout <<" NDF boltz = "<< ndf_bol << endl; 
cout <<" Chi2 / ndf = " << chindf_bol << endl;
*/
/*
H1_pp->Fit("levy","I","SAMES HIST",a,b);
double chi2_levy = func3->GetChisquare();
double ndf_levy = func3->GetNDF();
double chindf_levy = chi2_levy/ndf_levy ;
cout <<" Chi square levy = "<< chi2_levy << endl;
cout <<" NDF levy = "<< ndf_levy << endl; 
cout <<" Chi2 / ndf = " << chindf_levy << endl;
*/
//H1_pp->Fit("power_law_Five"," I L","SAMES HIST",a,b);
//H1_pp->Fit("power_law_Seven","+ R ","SAMES HIST",0,20);

H1_pp -> Fit("blast_wave", " I V","SAMES HIST",a,b);

//H1_pp->GetFunction("expo_law")->SetLineColor(kViolet);
//H1_pp->GetFunction("boltzmann")->SetLineColor(kGreen);
//H1_pp->GetFunction("levy")->SetLineColor(kCyan);
//H1_pp->GetFunction("power_law_Five")->SetLineColor(kMagenta);
//H1_pp->GetFunction("power_law_Seven")->SetLineColor(kRed);

H1_pp->GetFunction("blast_wave")->SetLineColor(kBlue);

   TLatex latex;
   latex.SetNDC();
   latex.SetTextSize(0.020);
   latex.SetTextAlign(13);  //align at top
   latex.SetTextFont(60);
   
/* 
   latex.DrawLatex(0.5,0.65, TString::Format("Parametre Levy C= %g", func3->GetParameter(0)));
   latex.DrawLatex(0.5,0.625, TString::Format("Parametre Levy m= %g GeV/c", func3->GetParameter(1)));
   latex.DrawLatex(0.5,0.60, TString::Format("Parametre Levy T= %g GeV", func3->GetParameter(2)));
   latex.DrawLatex(0.5,0.575, TString::Format("Parametre Levy n= %g ", func3->GetParameter(3)));
*/ 
  
 /*  latex.DrawLatex(0.5,0.65, TString::Format("Parametre expo C= %g", func1->GetParameter(0)));
   latex.DrawLatex(0.5,0.625, TString::Format("Parametre expo m_0= %g GeV/c", func1->GetParameter(1)));
   latex.DrawLatex(0.5,0.60, TString::Format("Parametre expo T= %g GeV/c", func1->GetParameter(2)));
 
   latex.DrawLatex(0.5,0.45, TString::Format("Parametre boltzmann C= %g", func2->GetParameter(0)));
   latex.DrawLatex(0.5,0.425, TString::Format("Parametre boltzmann m_0= %g GeV/c", func2->GetParameter(1)));
   latex.DrawLatex(0.5,0.40, TString::Format("Parametre boltzmann T= %g GeV/c", func2->GetParameter(2)));
 */






T1->Update();
		// Configuration de la légende
auto legend = new TLegend(0.4,0.7,0.6,0.9);
//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
legend->AddEntry(H1_pp,"Data","lep");
//legend->AddEntry(func1,"Exponential law fit","l");
//legend->AddEntry(func2,"Boltzmann law fit","l");
//legend->AddEntry(func3,"Levy-Tsallis fit","l");
//legend->AddEntry(func4, "Power law fit from [5]", "l");
//legend->AddEntry(func5, "Power law fit from [7]", "l");
//legend->AddEntry(func6,"Blast-wave model","l");
legend->Draw("SAMES");



T1->Write();


//cout << "chi2 exponential law " << chi2_expo << endl;
//cout << "chi2 boltzmann law " << chi2_boltzmann << endl;



// TEST SIMPSON
/*
double integration_expo = simpson(expo_law , func1->GetParameters() ,0,b,10000);
cout << "integration expo " << integration_expo << endl;

double integration_boltzmann = simpson(boltzmann , func2->GetParameters() ,0,b,10000);
cout << "integration boltzmann " << integration_boltzmann << endl;
*/

/*
TF1 *func_test = new TF1("test",test,0,10,2);
func_test -> SetParameter(0,3);
func_test -> SetParameter(1,2);
double integration_test = simpson(test, func_test->GetParameters(),0,1,1000);
cout << "intergation test sur x[0] " << integration_test << endl; 
*/
//TF1 *func_test2 = new TF1("blast_wave",blast_wave,0,10,2);



/*
// Tests de la méthode d'intégration pour fonction à 2 variables
// sur une fonction affine "test"
double *par = new double;
par[0] = 1;

double integ_test2 = simpson(test, par , 0 , 1 , 1000);
cout << "integ test 2 " << integ_test2 << endl;


TCanvas *canv1 = new TCanvas("Canvas 1" , " CANVAS_1_PP1 " , 200 , 200 );
canv1->SetGrid();
canv1->cd();
TF1 *ftest = new TF1("test", test, 0 , 4 ,1);
ftest -> SetParameter(0,1);
ftest -> Draw();
double integration_test = simpson(test, ftest->GetParameters(),0,1,1000);
cout << "integ test  " << integration_test << endl;

canv1->Write();
*/

/*
D1->cd();

func3->SetParameter(0,func3->GetParameter(0)) ;
func3->SetParameter(1,func3->GetParameter(1)) ;
func3->SetParameter(2,func3->GetParameter(2)) ;
func3->SetRange(0,20);
func3->SetNpx(10000);
func3->Draw();
cout << " Scale after = " << H1_pp->Integral("width") << endl;
double_t b= func3->Integral(0,20,1E-12);
cout << " integral value = " << b << endl;


D1->Write();
*/

/*
T2->cd();
H2_pp->Draw("][ E1");
H2_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H2_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H2_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H2_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H2_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H2_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();
T2->Write();
T3->cd();
H3_pp->Draw("][ E1");
H3_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H3_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H3_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H3_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H3_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H3_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();
T3->Write();
 
T4->cd();
H4_pp->Draw("][ E1");
H4_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H4_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H4_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H4_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H4_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H4_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();
T4->Write();
T5->cd();
H5_pp->Draw("][ E1");
H5_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H5_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H5_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H5_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H5_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H5_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();
T5->Write();
T6->cd();
H6_pp->Draw("][ E1");
H6_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H6_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H6_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H6_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H6_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H6_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();
T6->Write();
T7->cd();
H7_pp->Draw("][ E1");
H7_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H7_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H7_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H7_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H7_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H7_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();
T7->Write();
T8->cd();
H8_pp->Draw("][ E1");
H8_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H8_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H8_pp->Fit("Plaw","L + I R B","SAME HIST", 0 , 20);
H8_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H8_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H8_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();
T8->Write();
T9->cd();
H9_pp->Draw("][ E1");
H9_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H9_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H9_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H9_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H9_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H9_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();
T9->Write();
T10->cd();
H10_pp->Draw("][ E1");
H10_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H10_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H10_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H10_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H10_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H10_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();
*/

OutputHisto->Close();




cout << " fin " << endl ;
}
